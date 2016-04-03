from read_until import ReadUntil
import time
import errno
from socket import error as socket_error
import threading
#import MySQLdb
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
import subprocess
import re
import logging
import platform

from ruutils import process_model_file,query_yes_no,send_message


######################################################
def get_seq_len(ref_fasta):
    seqlens=dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        seq=record.seq
        seqlens[record.id]=len(seq)
    return seqlens

#######################################################################

def process_ref_fasta(ref_fasta,model_kmer_means,kmer_len):
    print "processing the reference fasta."
    kmer_means=dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id]=dict()
        kmer_means[record.id]["F"]=list()
        kmer_means[record.id]["R"]=list()
        kmer_means[record.id]["Fprime"]=list()
        kmer_means[record.id]["Rprime"]=list()
        print "ID", record.id
        print "length", len(record.seq)
        print "FORWARD STRAND"

        seq = record.seq
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))
            #if model_kmer_means[kmer]:
                #print x, kmer, model_kmer_means[kmer]

        print "REVERSE STRAND"
        seq = revcomp = record.seq.reverse_complement()
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["R"].append(float(model_kmer_means[kmer]))

        kmer_means[record.id]["Fprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["F"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=True, copy=True)
    return kmer_means


#######################################################################
def squiggle_search2(squiggle,channel_id,read_id,kmerhash,seqlen,args):
    result=[]

    for ref in kmerhash:
        #print "ss2",ref
        queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)

        dist, cost, path = mlpy.dtw_subsequence(queryarray,kmerhash[ref]['Fprime'])
        result.append((dist,ref,"F",path[1][0],ref))
        dist, cost, path = mlpy.dtw_subsequence(queryarray,kmerhash[ref]['Rprime'])
        result.append((dist,ref,"R",path[1][0],ref))


    return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4]

######################################################################
def extractsquig(events):
    squiggle=list()
    for event in events:
        squiggle.append(event.mean)
    return(squiggle)

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

#######################################################################
def go_or_no(seqid,direction,position,seqlen,args):
    for sequence in args.targets:
        if args.verbose is True: print sequence
        start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
        stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
        length = seqlen[seqid]
        if args.verbose is True:
            print start,stop,length
            print sequence.split(':', 1 )[0]
            print type(seqid)
        #We note that the average template read length is 6kb for the test lambda dataset. Therefore we are interested in reads which start at least 3kb in advance of our position of interest
        balance = args.length/2
        if seqid.find(sequence.split(':', 1 )[0]) >= 0:
            if args.verbose is True: print "Found it"
            if direction == "F":
                if args.verbose is True: print "Forward Strand"
                if position >= ( start - balance ) and position <= stop:
                    return "Sequence"
            elif direction == "R":
                if args.verbose is True: print "Reverse Strand"
                if position >= ( length - stop - balance) and position <= ( length - start ):
                    return "Sequence"
    return "Skip"

###################
def mp_worker((channel_id, data,kmerhash,seqlen,readstarttime,args)):
    if args.verbose is True:
        print "worker running"
        print "channel_id:",channel_id
        print "readstarttime:",readstarttime
    #readnumber = ()
    #print data.read_number
    try:
    #    print "in try"
        if args.verbose is True:
            print "read number", data.read_number
        readnumber = data.read_number
    except Exception, err:
        if args.verbose is True:
            err_string="Read Number Error : %s" % ( err)
        next
    try:
        if args.verbose is True:
            print "read id" ,data.read_id
        readnumber = data.read_id
    except Exception, err:
        if args.verbose is True:
            err_string="Read Number Error : %s" % ( err)
        next
    #print "set readnumber", readnumber
    if args.verbose is True:
        print "read number is", readnumber
    if ((time.time()-readstarttime) > args.time):
        if args.verbose is True:
            print "We have a timeout"
        #logging.info('%s,%s,%s,%s', channel_id, data.read_id, 'TOT',data.events[0].start)
        return 'timeout',channel_id,readnumber,data.events[0].start
    elif (args.skip is True and int(channel_id) %2 == 0) :
        if args.verbose is True:
            print "Even numbered channel so skip"
        return 'evenskip',channel_id,readnumber,data.events[0].start
    else:
        if args.verbose is True:
            print "Odd numbered channel - yay!"
        try:
            squiggle = extractsquig(data.events)
            if args.verbose is True:
                print squiggle
            squiggleres = squiggle_search2(squiggle,channel_id,readnumber,kmerhash,seqlen,args)
            if args.verbose is True: print squiggleres
            result = go_or_no(squiggleres[0],squiggleres[2],squiggleres[3],seqlen,args)
            if args.verbose is True: print "Result ",result
            return result,channel_id,readnumber,data.events[0].start,squiggleres
        except Exception, err:
            err_string="Time Warping Stuff : %s" % ( err)
            print >>sys.stderr, err_string
####################






def run_analysis():
    analyser = MyAnalyser()
    host = "ws://"+str(args.ip)+":"+str(args.port)+"/"
    setup_conditions = {"ignore_first_events": 100, "padding_length_events": 0,
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
            #print "looping"
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

    __version__ = "1.0"
    __date__ = "30th March 2016"

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

        p = multiprocessing.Pool(args.procs)
        fasta_file = args.fasta
        seqlen = get_seq_len(fasta_file)
        #print type(seqlen)
        print seqlen
        model_file = args.temp_model
        model_kmer_means,kmer_len=process_model_file(model_file)
        #model_kmer_means = retrieve_model()
        #global kmerhash
        kmerhash = process_ref_fasta(fasta_file,model_kmer_means,kmer_len)

        print type (kmerhash)


        def __init__(self):
            self.current_unblock_map = LockedDict()
            #self.p = multiprocessing.Pool(48)
            #for id in kmerhash:
            #for ref in kmerhash[id]:
            #    print id,ref

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
            #print "async called"
            if args.verbose is True: print "Channels length",len(channels)
            d=list()
            if args.verbose is True: print "Checking Channels"
            for channel_id, data in channels.iteritems():
                #print channel_id
                MyAnalyser.p.apply_async(mp_worker, args = ((channel_id,data,MyAnalyser.kmerhash,MyAnalyser.seqlen,time.time(),args), ), callback = self.mycallback)

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

    messagesend = True
    try:
        if messagesend: send_message("minoTour software is implementing read until on this version of minKNOW and will reject reads being sequenced. If you don\'t want this to happen, you should ensure ws_event_sampler is not running! You proceed at your own risk.",args.ip)
    except:
        messagesend = False
        next
    print "                                    `       "
#    print "              ,                 ;       "
    print "              ;`               ,;       "
#    print "              :;;:.        `.;;;`       "
    print "               :;;;;;;;;;;;;;;;,        "
#    print "                 `.:;;;;;;;:.`          "
    print "                  , .;;;;;, ,           "
#    print "             #@@@@@ ;;;;;;; @@@@@#      "
    print "           @@@@@@@@ ;;;;;;; @@@@@@@@    "
#    print "          @@@@@@@@@ ,;;;;;: @@@@@@@@@`  "
    print "         @@@@@@@@@@# ;;;;; +@@@@@@@@@@  "
#    print "        ,@@@@@@@@@@@ ,;;;, @@@@@@@@@@@; "
    print "        #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@ "
#    print "        #@@@@@   @@@@@   @@@@@   @@@@@@ "
    print "        .@@@@@    @@@@@@@@@@@    @@@@@: "
#    print "         @@@@@     @@@@@@@@@     @@@@@  "
    print "         .@@@@`    @@@@@@@@@     @@@@,  "
#    print "          +@@@@ @@  @@@@@@@  #@`@@@@#   "
    print "           '@@@@@@+ @@@@@@@ '@@@@@@+    "
#    print "            `@@@@@@ '@@@@@+ @@@@@@`     "
    print "              ;@@@#  @@@@@. +@@@; inoTour read until routines.       "
#    print "                 #   +++++`  @          "
    print "                   .;;;;;;;,            "
#    print "                  .;;;;;;;;;,           "
    print "                  ;;;.   .;;;`          "
#    print "                 ;;;.     `;;;          "
    print "                  ;;       ;;`          "
#    print "                  :;       ;;           "
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
            run_analysis()
        except socket_error as serr:
            if serr.errno != errno.ECONNREFUSED:
                raise serr
        print "Hanging around, waiting for the server..."
        time.sleep(5) # Wait a bit and try again
