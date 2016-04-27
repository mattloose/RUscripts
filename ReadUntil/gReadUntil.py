from read_until import ReadUntil
import time
import errno
from socket import error as socket_error
import threading
import sys, os, re
from Bio import SeqIO
from StringIO import StringIO
import string
import mlpy
import sklearn.preprocessing
import random
import math
import csv
import numpy as np
import array as ar
import configargparse
import shutil
import pickle
import multiprocessing
import ctypes
import subprocess
import re
import logging
import platform


from ruutils import process_model_file,query_yes_no,send_message,process_ref_fasta,get_seq_len,squiggle_search2,extractsquig,go_or_no,genome_worker

class LockedDict(dict):
    """
    A dict where __setitem__ is synchronised with a new function to
    atomically pop and clear the map.
    """
    def __init__(self, *args, **kwargs):
        self.lock = threading.Lock()
        super(LockedDict, self).__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        with self.lock:
            super(LockedDict, self).__setitem__(key, value)

    def pop_all_and_clear(self):
        with self.lock:
            d=dict(self) # take copy as a normal dict
            super(LockedDict, self).clear()
            return d

def run_analysis(args,analyser):
    #analyser = MyAnalyser(args)
    host = "ws://"+str(args.ip)+":"+str(args.port)+"/"
    setup_conditions = {"ignore_first_events": 75, "padding_length_events": 0,
                        "events_length": 250, "repetitions": 1}

    state=RunningState()
    with ReadUntil(host=host,
                   setup_conditions=setup_conditions,
                   data_received=analyser.apply_async_with_callback,
                   connection_closed=state.closed) as my_client:
        try:
            my_client.start()
        except Exception,err:
            print err
        print "Client connection started. Beginning unblock loop..."
        while state.keep_running:
            print "looping"
            unblock_now = analyser.next_unblock_map()
            if len(unblock_now)>0:
                if args.verbose is True: print "Unblocking channels: ", unblock_now.keys()
                print time.strftime('%Y-%m-%d %H:%M:%S'),
                print "Unblocking ",len(unblock_now.keys())
                my_client.unblock(unblock_now)

            # Throttle rate at which we make unblock controls. Although
            # unblocks should be timely, it is more efficient on the network
            # and on the hardware to unblock a bigger list of channels at once.
            time.sleep(1)
        print "...unblock loop ended. Connection closed."

class MyAnalyser:
    """
    Analyses recent data for each channel and amends the list of the channel
    ids that should be unblocked. ReadUntil should call data_received whenever
    some new data is available. It may call it several times concurrently and
    possibly with updates to the same channel if in the setup conditions we
    have requested more than 1 update for a given read.  We react by keeping up
    to date a current_unblock_set. Then we regularly make a separate unblock
    call using that list in another thread.

    This demonstrates one way of arranging code to analyse and unblock.
    How best to handle concurrent overlapping updates and when to unblock is an
    open problem left to the client script to resolve which will depend on what
    it is trying to achieve.
    """
    def __init__(self,args,p,seqids,threedarray):
        self.current_unblock_map = LockedDict()

    def mycallback(self, actions):
        if actions[0] == "Skip":
            if args.verbose is True: print "Should skip"
            self.current_unblock_map[actions[1]]=actions[2]
            logging.info('%s,%s,%s,%s,%s,%s,%s,%s', actions[1], actions[2], 'REJ',actions[3],actions[4][0],actions[4][1],actions[4][2],actions[4][3])
        elif actions[0] == "timeout":
            if args.verbose is True: print "Read timeout"
            logging.info('%s,%s,%s,%s', actions[1], actions[2], 'TOT',actions[3])
            #logging.info('%s,%s,%s', actions[1],actions[2],"TOT")
        elif actions[0] == "evenskip":
            if args.verbose is True: print "Even Skip"
            logging.info('%s,%s,%s,%s', actions[1], actions[2], 'EVE',actions[3])
        else:
            if args.verbose is True: print "Sequencing from ",action[1],action[2]
            #print actions[1], actions[2], 'SEQ',actions[3],actions[4][0],actions[4][1],actions[4][2],actions[4][3]
            logging.info('%s,%s,%s,%s,%s,%s,%s,%s', actions[1], actions[2], 'SEQ',actions[3],actions[4][0],actions[4][1],actions[4][2],actions[4][3])


    def apply_async_with_callback(self, channels):
        if args.verbose is True: print "Channels length",len(channels)
        d=list()
        if args.verbose is True: print "Checking Channels"
        for channel_id, data in channels.iteritems():
            p.apply_async(genome_worker, args = ((channel_id,data,time.time(),args,seqlen,seqids,threedarray), ), callback = self.mycallback)


    def next_unblock_map(self):
       """
       Returns current map of channel_id to read_id that should be unblocked but
       also clears the map in the assumption that the action to unblock will be
       carried out straightaway.
       """
       return self.current_unblock_map.pop_all_and_clear()

class RunningState:
    def __init__(self):
        self.keep_running=True

    def closed(self, *args):
        self.keep_running=False

if __name__ == "__main__":
    multiprocessing.freeze_support()
    manager=multiprocessing.Manager()
    global oper

    oper = platform.system()
    if oper is 'Windows':  # MS
        oper = 'windows'
    else:
        oper = 'linux'  # MS


    ## linux version
    if (oper is "linux"):
            config_file = os.path.join(os.path.sep, os.path.dirname(os.path.realpath('__file__')), 'amp.config')

    ## linux version
    if (oper is "windows"):
            config_file = os.path.join(os.path.sep, os.path.dirname(os.path.realpath('__file__')), 'ampW.config')

    __version__ = "1.2"
    __date__ = "8th April 2016"

    parser = configargparse.ArgParser(description='gReadUntil.py: A program providing read until for genome sequences with the Oxford Nanopore minION device. This program will ultimately be driven by minoTour to enable selective remote sequencing. This program is partly based on original code generously provided by Oxford Nanopore Technologies.')
    parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file describing the reference sequence for your organism.")
    parser.add('-targets', nargs = '*', dest='targets',required=True, help = 'Positional IDs to enrich for in the form seqid:start-stop . Can be space seperated eg: J02459:10000-15000 J02459:35000-40000')
    parser.add('-procs', '--proc_num', type=int, dest='procs',required=True, help = 'The number of processors to run this on.')
    parser.add('-t', '--time', type=int, dest='time', required=True, default=300, help="This is an error catch for when we cannot keep up with the rate of sequencing on the device. It takes a finite amount of time to process through the all the channels from the sequencer. If we cannot process through the array quickly enough then we will \'fall behind\' and lose the ability to filter sequences. Rather than do that we set a threshold after which we allow the sequencing to complete naturally. The default is 300 seconds which equates to 9kb of sequencing at the standard rate.")
    parser.add('-m', '--model',type=str, required=True, help = 'The appropriate template model file to use', dest='temp_model')
    parser.add('-ip', '--ip-address', type=str ,dest='ip',required=False,default="127.0.0.1", help="The IP address of the machine running minKNOW.")
    parser.add('-p', '--port', type=int, dest='port', default=None,required=True, help='The port that ws_event_sampler is running on.' )
    parser.add('-log', '--log-file', type=str, dest='logfile', default='readuntil.log', help="The name of the log file that data will be written to regarding the decision made by this program to process read until.")
    parser.add('-length', '--library-length', type=int, dest='length', required=False, default=0, help="Provide the average expected length of your library. This offset will be applied to reads that are likely to extend into your region of interest on either strand.")
    parser.add('-skip','--skip_even', action='store_true', help="If set, this will allow all reads from even numbered channels to be sequenced regardless of where they map. This provides an internal control.", default=False,dest='skip')
    parser.add('-v', '--verbose-true', action='store_true', help="Print detailed messages while processing files.", default=False, dest='verbose')
    parser.add_argument('-ver', '--version', action='version',version=('%(prog)s version={version} date={date}').format(version=__version__,date=__date__))
    #global args
    args = parser.parse_args()

    ###Check files

    if not os.path.isfile(args.fasta):
        print "\n**! The fasta file cannot be found - please check your path and file name under the -fasta option.\n"
        sys.exit()

    if not os.path.isfile(args.temp_model):
        print "\n**! The model file cannot be found - please check your path and file name under the -model option.\n"
        sys.exit()


    fasta_file = args.fasta
    #global seqlen
    seqlen = get_seq_len(fasta_file)
    #print type(seqlen)
    print seqlen
    model_file = args.temp_model
    global model_kmer_means
    global kmer_len
    model_kmer_means,kmer_len=process_model_file(model_file)
    seqids,threedarray = process_ref_fasta(fasta_file,model_kmer_means,kmer_len)
    #print "init kmerhash",type(kmerhash)

    print type(threedarray)


    p = multiprocessing.Pool(args.procs)
    analyser = MyAnalyser(args,p,seqids,threedarray)
    messagesend = True
    try:
        if messagesend: send_message("minoTour software is implementing read until on this version of minKNOW and will reject reads being sequenced. If you don\'t want this to happen, you should ensure ws_event_sampler is not running! You proceed at your own risk.",args.ip)
    except:
        messagesend = False
        next
    print "                                    `       "
    print "              ;`               ,;       "
    print "               :;;;;;;;;;;;;;;;,        "
    print "                  , .;;;;;, ,           "
    print "           @@@@@@@@ ;;;;;;; @@@@@@@@    "
    print "         @@@@@@@@@@# ;;;;; +@@@@@@@@@@  "
    print "        #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@ "
    print "        .@@@@@    @@@@@@@@@@@    @@@@@: "
    print "         .@@@@`    @@@@@@@@@     @@@@,  "
    print "           '@@@@@@+ @@@@@@@ '@@@@@@+    "
    print "              ;@@@#  @@@@@. +@@@; inoTour read until routines.       "
    print "                   .;;;;;;;,            "
    print "                  ;;;.   .;;;`          "
    print "                  ;;       ;;`          "
    print "  Welcome to the .;;:     ,;;,          "
    print ""

    print "This script WILL implement read until.\nIf you proceed it is at your own risk.\n\n"
    if not query_yes_no("Seriously - are you happy to proceed? Entering yes will make it your fault..."):
        sys.exit()
    print "***********************************************************************************************"
    if args.skip is True:
        print "**** This version of the code will not process reads derived from even numbered channels.  ****"
    else:
        print "****        This version of the code will process reads regardless of channel.             ****"
    print "***********************************************************************************************"
    logging.basicConfig(format='%(levelname)s:%(message)s',filename=args.logfile, filemode='w', level=logging.INFO    )
    current_time = time.time()
    print current_time

    while 1:
        try:
            print "Running Analysis"
            run_analysis(args,analyser)
        except socket_error as serr:
            if serr.errno != errno.ECONNREFUSED:
                raise serr
        print "Hanging around, waiting for the server..."
        time.sleep(5) # Wait a bit and try again
