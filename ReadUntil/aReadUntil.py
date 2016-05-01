#!/usr/bin/env python
from read_until import ReadUntil
import time
import errno
from socket import error as socket_error
import socket
import threading, thread
import sys, os, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
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
import logging
import urllib2
import json
import hashlib
from watchdog.observers.polling import PollingObserver as Observer
from watchdog.events import FileSystemEventHandler
import datetime
import h5py
import psutil
from Queue import Queue
from cProfile import run
from pstats import Stats
#from blinkstick import blinkstick
import platform

from ruutils import process_model_file,query_yes_no,send_message,execute_command_as_string,get_seq_len,squiggle_search2,get_custom_fasta
from ruutils import file_dict_of_folder,scaleLocally,scale,merge_ranges,correctposition,extractsquig,die_nicely,process_ref_fasta,check_files,checkfasta


######################################################

def get_amplicons(amplicons,seqlengths,args):
    if (args.verbose is True):
        #print "Groking amplicons"
        print "ids is of type", type(amplicons)
    for sequence in amplicons:
        if (args.verbose is True):
            print sequence
        start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
        stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
        if (args.verbose is True):
            print start
            print stop
        REVERSE_stop = seqlengths[sequence.split(':', 1 )[0]]-start
        REVERSE_start = seqlengths[sequence.split(':', 1 )[0]]-stop
        if (args.verbose is True):
            print REVERSE_stop
            print REVERSE_start

######################################################
def process_ref_fasta2(ref_fasta,model_kmer_means):
    if (args.verbose is True):
        print "processing the reference fasta."
    kmer_len=5
    kmer_means=dict()

    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id]=dict()
        kmer_means[record.id]["F"]=list()
        seq = record.seq
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))
    return kmer_means

#######################################################################




def runProcess(exe):
    p=subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while(True):
        retcode= p.poll()
        line=p.stdout.readline()
        yield line
        if(retcode is not None):
            break





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



def channeldictexpire(channeldict):
    '''
    Function loops through a dictionary to test if a read is too old to plausibly exist anymore.
    Removes reads which must either be blocked or failed.
    '''
    for key in channeldict.keys():
        for item,value in channeldict[key].items():
            if time.time()-value > readtimeout:
                del channeldict[key]



def processchanneldict(channeldict):
    dtype = [('channel',int),('amplicon',int),('startsamples',float)]
    arr=[]
    for a,_ in channeldict.items():
        for k,_ in channeldict[a].items():
            arr.append((a,k,channeldict[a][k]))
    values=np.array(arr,dtype=dtype)
    c=list(np.sort(values,order='amplicon'))
    d=map(lambda (x,y,z):y,c)
    lut=dict()
    for x in set(d):
        lut[x]=d.count(x)
    return lut

###################
def amplicon_worker((channel_id, data,seqlen,readstarttime,seqids,threedarray,startdict,enddict,procampres,procampresrej,procampresseq,currentseq,procampresdone,procampresfile,procampres2d,args,customdepthslist,procamprestime,eventdict,channeldict)):
    squiggle = extractsquig(data.events)
    currentlysequence=0
    channel_id_orig = channel_id

    if isinstance(channel_id, basestring):
        channel_id = int(channel_id)
    try:
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
    if channel_id in currentseq:
        procampresdone[currentseq[channel_id]]+=1
        del currentseq[channel_id]
    if ((time.time()-readstarttime) > args.time):
        if channel_id in procamprestime:
            procamprestime[channel_id]+=1
        else:
            procamprestime[channel_id]=1
        return 'timeout',channel_id_orig,readnumber,data.events[0].start
    else:
        try:
            readprediction=dict()
            if (args.verbose is True):
                print "Read start time",readstarttime
            try:
                (seqmatchnameR,distanceR,frR,rsR,reR,qsR,qeR) = squiggle_search2(squiggle[0:len(squiggle)],channel_id,readnumber,args,seqids,threedarray,len(squiggle))
            except:
                print "MLPY FAILURE"
                sys.exit()
            if (args.verbose is True):
                print seqmatchnameR,distanceR,frR,rsR,reR,qsR,qeR
            try:
                if (frR == "F"):
                    amplicon, value = min(startdict.items(), key=lambda (_, v): abs(v - rsR))
                else:
                    amplicon, value = min(enddict.items(), key=lambda (_, v): abs(v - rsR))
            except:
                print "Cryptic Code Failure"
                sys.exit()
            if (args.verbose is True):
                print amplicon, value
            try:
                try:
                    combiname = amplicon
                except:
                    print "combiname fail"
                try:
                    dict2=dict(channeldict)
                    curseq=processchanneldict(dict2)
                except:
                    print "Curseq assign failed"
                try:
                    if combiname in curseq:
                        currentlysequence=curseq[combiname]
                except:
                    print "setting currentlysequence failed"
            except:
                print "channel dict fail"
            if (args.verbose is True):
                print "Args Depth ",args.depth
                print "procampres ", procampres[combiname]
            try:
                if (args.goal == "file"):
                    testset = procampresfile
                elif (args.goal == "read"):
                    testset = procampresdone
                elif (args.goal == "2d"):
                    testset = procampres2d
                else:
                    print "Test set not recognised - quitting for sanity. -g should be set to either 'file','2d' or 'read' i.e -g read -g 2d -g file"
                    exit()
            except:
                print "Checking goal failed"
            try:
                if (args.depth >= 0):
                    correctionfactor=0
                    if args.precision:
                        correctionfactor = int(currentlysequence/2)
                    if (testset[combiname] + correctionfactor -args.deptherror >= args.depth):
                        if not args.inhibit:
                            result = "Skip"
                            if (args.verbose is True):
                                print "Skip", procampres[combiname]
                            procampresrej[combiname] +=1
                        else:
                            result = "Sequence"
                            procampresseq[combiname] +=1
                            currentseq[channel_id]=combiname
                    else:
                        result = "Sequence"
                        #result = "Skip"
                        if (args.verbose is True):
                            print "Sequence", procampres[combiname]
                        procampresseq[combiname] +=1
                        currentseq[channel_id]=combiname
                elif (len(customdepthslist)> 0):
                    correctionfactor=0
                    if args.precision:
                        correctionfactor = int(currentlysequence)
                    if ((testset[combiname]+correctionfactor-args.deptherror) >= int(customdepthslist[combiname-1])):
                        if not args.inhibit:
                            result = "Skip"
                            if (args.verbose is True):
                                print "Skip", procampres[combiname]
                            procampresrej[combiname] +=1
                        else:
                            result = "Sequence"
                            procampresseq[combiname] +=1
                            currentseq[channel_id]=combiname
                    else:
                        result = "Sequence"
                        if (args.verbose is True):
                            print "Sequence", procampres[combiname]
                        procampresseq[combiname] +=1
                        currentseq[channel_id]=combiname
            except:
                print "Sequence Decision Failed"
            try:
                procampres[combiname] += 1
                if (combiname not in readprediction):
                    readprediction[combiname]=dict()
                if (0 not in readprediction[combiname]):
                    readprediction[combiname][0]=dict()
                if (channel_id not in readprediction[combiname][0]):
                    readprediction[combiname][0][channel_id]=dict()
                readprediction[combiname][0][channel_id]["name"]=channel_id
                readprediction[combiname][0][channel_id]["matchdistance"]=distanceR
            except:
                print "Dictionary manipulation up the khyber"
            try:
                if args.simcheck:
                    corrected_channel = channel_id
                    if (result == "Sequence"):
                        f = h5py.File(args.watchdir+'/test_ch'+str(corrected_channel)+'_file'+str(readnumber)+'_strand_2d_'+str(combiname)+'.fast5','w')
                        f.create_group("/UniqueGlobalKey/channel_id").attrs['channel_number'] = str(corrected_channel)
                        f.create_group("/Analyses/EventDetection_000/Reads/Read_"+str(readnumber)).attrs['hairpin_found'] = 1
                    else:
                        f = h5py.File(args.watchdir+'/test_ch'+str(corrected_channel)+'_file'+str(readnumber)+'_strand_1d_'+str(combiname)+'.fast5','w')
                        f.create_group("/UniqueGlobalKey/channel_id").attrs['channel_number'] = str(corrected_channel)
                        f.create_group("/Analyses/EventDetection_000/Reads/Read_"+str(readnumber)).attrs['hairpin_found'] = 0
                    dataarray=[(event.start,event.length,event.mean,event.sd) for event in data.events]
                    numpydata = np.array(dataarray)
                    comp_type = np.dtype([('start', np.uint64), ('length', np.uint32), ('mean', np.float64), ('variance', np.float64)])
                    dataout = np.array(dataarray, dtype = comp_type)
                    dataset = f.create_dataset("/Analyses/EventDetection_000/Reads/Read_"+str(readnumber)+"/Events",(len(dataout),), comp_type)
                    if (args.verbose is True):
                        print len(dataout)
                    dataset[...] = dataout
                    #f.create_dataset('/Analyses/EventDetection_000/Reads/Read_'+str(data.read_id)+'/Events', data=numpydata)
                    f.close()
            except:
                print "File write failed"
            if result == "Sequence":
                try:
                    #channeldict[channel_id]={combiname:data.events[0].start}
                    channeldict[channel_id]={combiname:time.time()}
                except:
                    print "Channeldict function failed"
            return result,channel_id_orig,readnumber,data.events[0].start,data.events[0].mean,eventdict,combiname
            #return result,channel_id,data.read_id,data.events[0].start
        except Exception, err:
            err_string="Time Warping Stuff : %s" % ( err)
            print >>sys.stderr, err_string
            result = "TWF"
            print "Result is ", result
            return result,channel_id_orig,readnumber,data.events[0].start

#############################

def checkevents(meansquiggle,meantime,channel_id,eventdict):
    counter=0
    #correct for position 1 frame shift
    #channel = channel_id - 1
    channel = channel_id
    #print "Channel is",channel
    subevents = eventdict[channel]
    if (args.verbose is True):
        print "TYPE: ", type(subevents)
    for event in subevents:
        if (args.verbose is True):
            print "EVENT FROM DICT" ,event
            print "mean" , subevents[event]["mean"]
            print "amplicon", subevents[event]["prediction"]
            print "meansquiggle", meansquiggle
        if subevents[event]["mean"] in meansquiggle:
            if (args.verbose is True):
                print "FOUND match IT!",channel_id
            indexpos = meansquiggle.index(subevents[event]["mean"])
            if meantime[indexpos] == subevents[event]["start"]:
                if (args.verbose is True):
                    print "And its the right time too!"
                    print "meantime",meantime[indexpos],"start",subevents[event]["start"],"prediction",subevents[event]["prediction"],"channel",channel,"indexpos",indexpos,"event",event
                #del subevents[event]
                #eventdict[channel_id]=subevents
                #print "Trying to return 1", event
                return ("1",event)
        else:
            continue
	print "read not matched"
    return ("0",event)
            #print "Trying to return 0", event


def process_fast5((filepath,eventdict,procampresfile,procampres2d)):
    hdf = h5py.File(filepath, 'r')
    #print filepath
    channel_id=hdf["UniqueGlobalKey/channel_id"].attrs['channel_number']
    if (args.verbose is True):
        print "Processing:", channel_id
    corrected_channel = int(channel_id)
    corrected_channel = str(corrected_channel)
    #print "Corrected Channel Type: ",type(corrected_channel),"Value: ",corrected_channel
    #print "***badger***"
    for element in hdf['Analyses/EventDetection_000/Reads']:
        #print element
        if (args.verbose is True):
            print element
        for thing in hdf['Analyses/EventDetection_000/Reads/'+element]:
            #Here we want to recover a list of means for the read.
            hairpin_found = hdf['Analyses/EventDetection_000/Reads/'+element].attrs['hairpin_found']
            #print "***HAIRPIN***",hairpin_found
            meansquiggle=list()
            meantime=list()
            for event in hdf['Analyses/EventDetection_000/Reads/'+element+'/'+thing][:350]:
            #for event in hdf['Analyses/EventDetection_000/Reads/'+element+'/'+thing][45:55]:
                # print but[2]
                meansquiggle.append(float(event[2]))
                meantime.append(float(event[0]))
            #print meansquiggle
            #print meantime
            #print "Event Dict" , eventdict
            #for thing in eventdict:
            #    print "Eventdict key type: ", type(thing),"value: ",thing
            #print "meansqig", meansquiggle
            #print "meantime", meantime
            #print "channel_id",channel_id
            if corrected_channel in eventdict:
                test,event = checkevents(meansquiggle,meantime,corrected_channel,eventdict)
                if test:
                    if hairpin_found == 1:
                        procampres2d[eventdict[corrected_channel][event]["prediction"]]+=1
                    procampresfile[eventdict[corrected_channel][event]["prediction"]]+=1
                    subevents = eventdict[corrected_channel]
                    if (event in subevents.keys()):
                        del subevents[event]
                    eventdict[corrected_channel]=subevents
                else:
                    print "READ DOES NOT MATCH!!"
            else:
                print "not found"

    hdf.close()
    #Now we want to move the file to a new folder for processing by metrichor
    #print "Trying to move file."
    target = os.path.join(args.watchdir,'done')
    filetocheck = os.path.split(filepath)
    #print target
    #print filetocheck
    destfile = os.path.join(target,filetocheck[1])
    #print destfile
    if not os.path.exists(target):
        os.makedirs(target)
    try:
        #os.symlink(sourcefile, destfile)
        shutil.move(filepath,destfile)
    except Exception, err:
        print "File Move Failed",err

####################
class MyHandler(FileSystemEventHandler):


    def __init__(self,args):
        #We dont care about already existing files when we start this running
        #self.creates=file_dict_of_folder(args.watchdir)
        self.args=args
        self.creates=dict()
        self.processed=dict()
        self.running = True
        #threading.set_profile(stats)
        t = threading.Thread(target=self.processfiles)

        t.daemon = True
        try:
            t.start()
        except (KeyboardInterrupt, SystemExit):
            t.stop()
            die_nicely(oper)

    def lencreates(self):
        return len(self.creates)

    def lenprocessed(self):
        return len(self.processed)


    def processfiles(self):
        everyten=0
        while self.running:
            #ts = time.time()
            #blinkgo(bstick)
            #blinker.setjob("single","dimgrey")
            #print datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'), "CACHED:", self.lencreates(), "PROCESSED:",  self.lenprocessed()
            #if (args.verbosegen is True):
            #    print "obs:",procampres
            #    print "rej:",procampresrej
            #    print "seq:",procampresseq
            #print time.ctime(),": Obs:",sum_dict(procampres),"Rej:",sum_dict(procampresrej),"Seq:",sum_dict(procampresseq),"Done:",sum_dict(procampresdone),"File:",sum_dict(procampresfile),"2D:",sum_dict(procampres2d)
            #print "Obs Details:",procampres
            #print "Rej Details:",procampresrej
            #print "Seq Details:",procampresseq
            #print "DoneDetails:",procampresdone
            #print "FileDetails:",procampresfile
            #print "2D  Details:",procampres2d
            #print channeldict
            #dict2=dict(channeldict)
            #print "CurSDetails:",processchanneldict(dict2)
            ## Now we try and clean the channeldict of 'old reads'
            ##channeldictexpire(channeldict)
            for fast5file, createtime in sorted(self.creates.items(), key=lambda x: x[1]):
                #tn=time.time()
                if (args.verbose is True):
                    print fast5file,time.ctime()
                #lets try letting the os worry about if the file is finished writing
                delaytime = 5
                #if args.simcheck:
                #    delaytime = int(ampliconmeanlength/args.speed)
                if ( int(createtime)+delaytime < time.time() ): # file created 5 sec ago, so should be complete. For simulations we make the time longer.
                    if (fast5file not in self.processed.keys() ):
                        if (args.verbose is True):
                            print "File seen :" + fast5file
                        del self.creates[fast5file]
                        self.processed[fast5file]=time.time()
                        try:
                            process_fast5((fast5file,eventdict,procampresfile,procampres2d))
                        except Exception, err:
                            err_string="Error processing fast5 file for checking read: %s : %s" %(fast5file, err)
                    if (fast5file in self.processed.keys() and fast5file in self.creates.keys() ):
                        del self.creates[fast5file]
            time.sleep(5)

    def on_created(self, event):
        #print "On Created Called"
        if ("muxscan" not in event.src_path and event.src_path.endswith(".fast5")):
            self.creates[event.src_path] = time.time()
        #self.total[event.src_path] = time.time()

    def on_moved(self,event):
        if ("muxscan" not in event.dest_path and event.dest_path.endswith(".fast5")):
            self.creates[event.dest_path] = time.time()




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


    def __init__(self,args,p,seqids,threedarray,endcollection,startcollection):
        self.current_unblock_map = LockedDict()
        self.startcollection=startcollection
        self.endcollection=endcollection
        #print threedarray
        #self.eventdict=dict()

    def mycallback(self, actions):
        try:
            actions[0]
        except NameError:
            print "not defined"

        if actions[0] == "Skip":
            self.current_unblock_map[actions[1]]=actions[2]
            #blinkreject(bstick)
            #blinker.setjob("solid","blue")
        elif actions[0] == "timeout":
            if not args.inhibit:
                self.current_unblock_map[actions[1]]=actions[2]
            #blinkTO(bstick)
            #blinker.setjob("solid","red")
            print "Read timeout"
            #logging.info('%s,%s,%s,%s', actions[1], actions[2], 'TOT',actions[3])
        #elif actions[0] == "evenskip":
        #    logging.info('%s,%s,%s,%s', actions[1], actions[2], 'EVE',actions[3])
        else:
            #blinker.setjob("solid","green")
            #blinkaccept(bstick)
            next
        if actions[0] is not "TWF":
            #print "ACTION SEEN ",actions[0]
            try:
                logging.info('Channel:%s,Read:%s,Time:%s,Obs:%s,Rej:%s,Seq:%s',actions[1],actions[2],actions[3],procampres,procampresrej,procampresseq)
                if (args.verbose is True):
                    print "Channel id : " , actions[1]
                    print "Read id : " , actions[2]
                    print "Start Time : " , actions[3]
                    print "Mean Current : " , actions[4]
                if (actions[1] not in eventdict):
                    eventdict[actions[1]]=dict()
                chanvals = eventdict[actions[1]]
                if (actions[2] not in chanvals):
                    chanvals[actions[2]]=dict()
                chanvals[actions[2]]["start"]=actions[3]
                chanvals[actions[2]]["mean"]=actions[4]
                chanvals[actions[2]]["prediction"]=actions[6]
                eventdict[actions[1]]=chanvals
            except:
                #print "Logging failed"
                next
        else:
            print "TWF FOUND HELP"


    def data_received(self, channels):
        d=list()
        for channel_id, data in channels.iteritems():
            d.append([channel_id,data,seqlen,time.time()])
        procdata=tuple(d)
        results = p.map_async(mp_worker, (procdata),chunksize=1,callback=self.mycallback)


    def apply_async_with_callback(self, channels):
        d=list()
        for channel_id, data in channels.iteritems():
            p.apply_async(amplicon_worker, args = ((channel_id,data,seqlengths,time.time(),seqids,threedarray,self.startcollection,self.endcollection,procampres,procampresrej,procampresseq,currentseq,procampresdone,procampresfile,procampres2d,args,customdepthslist,procamprestime,eventdict,channeldict), ), callback = self.mycallback)

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

def check_completion(procampresseq,args,customdepthslist):
    testdict=dict(procampresseq)
    result = True
    #print customdepthslist
    if (args.depth > 0):
        for i in testdict:
            if testdict[i]>=(args.depth-args.deptherror):
                result = True
            else:
                result = False
                break
    elif (len(customdepthslist) > 0):
        for i in testdict:
            #print i, testdict[i],(int(customdepthslist[i-1])-args.deptherror)
            if testdict[i]>=(int(customdepthslist[i-1])-args.deptherror):
                result = True
            else:
                result = False
                break
    return(result)

def sum_dict(dictionarytosum):
    testdict=dict(dictionarytosum)
    cumulative = 0
    for i in testdict:
        cumulative = cumulative + testdict[i]
    return(cumulative)


def run_analysis(args,procampres2d,customdepthslist):

    host = "ws://"+str(args.ip)+":"+str(args.port)+"/"
    #host = "ws://127.0.0.1:5999/"
    setup_conditions = {"ignore_first_events": 100, "padding_length_events": 0,
                    "events_length": 250, "repetitions": 1}

    state=RunningState()
    event_handler = MyHandler(args)
    observer = Observer()
    observer.schedule(event_handler, path=args.watchdir, recursive=False)
    observer.daemon=True
    try:
        observer.start()

        with ReadUntil(host=host,
                   setup_conditions=setup_conditions,
                   data_received=analyser.apply_async_with_callback,
                   connection_closed=state.closed) as my_client:
            my_client.start()
            messagesent=0
            print "Client connection started. Beginning unblock loop..."
            while state.keep_running:
                unblock_now = analyser.next_unblock_map()
                if (args.verbose is True):
                    print eventdict
                if (args.goal == "file"):
                    testset = procampresfile
                elif (args.goal == "read"):
                    testset = procampresdone
                elif (args.goal == "2d"):
                    testset = procampres2d
                else:
                    print "Test set not recognised - quitting for sanity. -g should be set to either 'file','2d' or 'read' i.e -g read -g 2d -g file"
                    die_nicely(oper)
                if check_completion(testset,args,customdepthslist):
                    #print "The minoTours work is done. You should stop sequencing now!"
                    if not args.inhibit and args.stoprun:
                        print "Sequencing complete. Now waiting for "+str(args.writetime)+" seconds to ensure reads are written correctly to disk before stopping minKNOW."
                        if (len(args.customdepth) > 0):
                            endmessage = '{"id":"1", "method":"user_message","params":{"content":"minoTour read until has detected coverage of at least '+str(args.customdepth)+'x on each amplicon and will stop minKNOW in '+str(args.writetime)+' seconds."}}'
                        else:
                            endmessage = '{"id":"1", "method":"user_message","params":{"content":"minoTour read until has detected coverage of at least '+str(args.depth)+'x on each amplicon and will stop minKNOW in '+str(args.writetime)+' seconds."}}'
                        results = execute_command_as_string(endmessage,args.ip,8000)
                        stoprun = '{"id":"1", "method":"stop_experiment","params":"null"}'
                        stopprotocol = '{"id":"1", "method":"stop_script","params":{"name":"MAP_48Hr_Sequencing_Run.py"}}'
                        time.sleep(args.writetime)
                        if args.stoprun:
                            results = execute_command_as_string(stoprun,args.ip,8000)
                            results = execute_command_as_string(stopprotocol,args.ip,8000)
                            endmessage2 = '{"id":"1", "method":"user_message","params":{"content":"The minoTours work is done - run stopped."}}'
                            print "The minoTours work is done. Sequencing will be stopped now!"
                            results = execute_command_as_string(endmessage2,args.ip,8000)
                        else:
                            endmessage2 = '{"id":"1", "method":"user_message","params":{"content":"The minoTours work is done - however the run has not been stopped automatically."}}'
                            print "The minoTours work is done but sequencing has not been stopped. To enable run stopping use the -s option."
                            results = execute_command_as_string(endmessage2,args.ip,8000)
                        die_nicely(oper)
                    else:
                        if (messagesent < 3):
                            endmessage = '{"id":"1", "method":"user_message","params":{"content":"minoTour software has detected coverage of at least '+str(args.depth)+'x on each amplicon and suggests that you stop minKNOW now."}}'
                            results = execute_command_as_string(endmessage,args.ip,8000)
                            print "You have reached your goal now. Sequencing should be stopped now!"
                            messagesent += 1
                        if (messagesent == 3):
                            die_nicely(oper)
                #print "test ublock len"
                if len(unblock_now)>0:
                    #if (args.verbosegen is True):
                    #blinker.setjob("light","green")
                    #blinkgo(bstick)
                    print "Unblocking ",len(unblock_now.keys())
                    #print "Unblocked"
                    my_client.unblock(unblock_now)
                else:
                    next

                ts = time.time()
                print datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'), "CACHED:", event_handler.lencreates(), "PROCESSED:",  event_handler.lenprocessed()
                if (args.verbosegen is True):
                    print "obs:",procampres
                    print "rej:",procampresrej
                    print "seq:",procampresseq
                print time.ctime(),": Obs:",sum_dict(procampres),"Rej:",sum_dict(procampresrej),"Seq:",sum_dict(procampresseq),"Done:",sum_dict(procampresdone),"File:",sum_dict(procampresfile),"2D:",sum_dict(procampres2d)
                print "Obs Details:",procampres
                print "Rej Details:",procampresrej
                print "Seq Details:",procampresseq
                print "DoneDetails:",procampresdone
                print "FileDetails:",procampresfile
                print "2D  Details:",procampres2d
                #print channeldict
                dict2=dict(channeldict)
                print "CurSDetails:",processchanneldict(dict2)
                ## Now we try and clean the channeldict of 'old reads'
                channeldictexpire(channeldict)


            # Throttle rate at which we make unblock controls. Although
            # unblocks should be timely, it is more efficient on the network
            # and on the hardware to unblock a bigger list of channels at once.
                time.sleep(1)
            print "...unblock loop ended. Connection closed."
    except (KeyboardInterrupt, SystemExit):
        print "catching a ctrl-c event"
        my_client.stop()
        observer.stop()
        observer.join()
        die_nicely(oper)


if __name__ == "__main__":
    #Supposedly handling multiprocessing problems
    multiprocessing.freeze_support()

    #Initialising a manager to support multiprocessing shared dictionaries
    manager = multiprocessing.Manager()
    ##run('main_function()', 'stats')
    ##stats=Stats('stats')
    ##threading.setprofile('stats')
    global oper

    oper = platform.system()
    if oper is 'Windows':  # MS
        oper = 'windows'
    else:
        oper = 'linux'  # MS


    ## linux version
    if (oper is "linux"):
        config_file = os.path.join(os.path.sep, os.path.dirname(os.path.realpath('__file__')), 'aReadUntil.config')

    ## linux version
    if (oper is "windows"):
        config_file = os.path.join(os.path.sep, os.path.dirname(os.path.realpath('__file__')), 'aReadUntil.config')

    __version__ = "1.2"
    __date__ = "8th April 2016"

    #global #blinker
    #blinker = ThreadingExample(r_led_count=8, max_rgb_value=128)
    #blinker.setjob("light","white")

    #Parsing arguments
    parser = configargparse.ArgParser(description='aReadUntil: A program providing read until with the Oxford Nanopore minION device. This program will ultimately be driven by minoTour to enable selective remote sequencing. This program is based on original code generously provided by Oxford Nanopore Technologies. Note that whilst some parameters can be set via a config file, the explicit parameters to stop a run (-s), prevent read until working (-i) and switch to presicion mode (-precison) can only be set via the command line.',default_config_files=[config_file])
    parser.add('-fasta', '--reference-fasta-file', type=str, dest='fasta', required=True, default=None, help="The fasta format file describing the reference sequence for your organism.")
    parser.add('-c', '--custom-genome',action='store_true', help="This will use a reduced search space genome to match against.",default=False,dest='custom')
    parser.add('-ids', '--ids', type=str, required=True, default=None, help="A file containing a list of amplicon positions defined for the reference sequence. 1 amplicon per line in the format fasta_sequence_name:start-stop e.g EM_079517:27-1938", dest='ids')
    parser.add('-d', '--depth',type=int, required=False, default=None, help = 'The desired coverage depth for each amplicon. Note this is unlikely to be achieved for each amplicon and should probably be an overestimate of the minimum coverage required.', dest='depth')
    parser.add('-cd', '--custom-depth',type=str, required=False,default="", help= 'A comma seperated list of custom depths for each amplicon. You must provide a coverage depth for each amplicon in the order they are presented in the ids file.',dest='customdepth')
    parser.add('-e', '--error',type=int, required=True, default=0, help = 'Set an error range for coverage depth.', dest='deptherror')
    parser.add('-procs', '--processor-number', type=int, dest='procs',required=True, help = 'The number of processors to run this on.')
    parser.add('-t', '--time', type=int, dest='time', required=True, default=300, help="This is an error catch for when we cannot keep up with the rate of sequencing on the device. It takes a finite amount of time to process through the all the channels from the sequencer. If we cannot process through the array quickly enough then we will \'fall behind\' and lose the ability to filter sequences. Rather than do that we set a threshold after which we allow the sequencing to complete naturally. The default is 300 seconds which equates to 9kb of sequencing at the standard rate.")
    parser.add('-m', '--template-model',type=str, required=True, help = 'The appropriate template model file to use', dest='temp_model')
    parser.add('-g', '--goal',type=str, required=True, help = 'The measure by which reads will be counted - either based on the presence of files ( -g file) or potential 2D files generated (-g 2d) or new reads generated ( -g read )', dest='goal')
    parser.add('-precision', action='store_true', help="This option will attempt to obtain exactly the number of reads required per amplicon. It is provided as a novelty to illustrate the theoretical level of control of the device. In reality it will slow down the time taken to reach a specific goal due to the possibility of reads failing and the delay in writing true reads to disk.", default=False,dest='precision')
    parser.add('-seq', '--seq-speed', type=int, default=30, required=False, help="This is the assumed sequencing speed. The default is set at 30b/s (the speed of the simulator). This should be configured to the appropriate value for your chemistry.", dest='speed')
    parser.add('-i', action='store_true', help="This will prevent read until from working but will otherwise report what is happening in the sequencer.", default=False,dest='inhibit')
    parser.add('-s', action='store_true', help="This will enable read until to stop your sequencing when it is complete.", default=False,dest='stoprun')
    parser.add('-ip', '--ip-address', type=str ,dest='ip',required=False,default="127.0.0.1", help="The IP address of the minKNOW machine.")
    parser.add('-p', '--port', type=int, dest='port', default=None,required=True, help='The port that ws_event_sampler is running on.' )
    parser.add('-wt', '--write-time', type=int, dest='writetime', default=15, help='If you are automatically stopping the minKNOW run, the stop command will wait n seconds after the last read has completed to ensure all files are written. Default value is 15.' )
    parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows). This folder must already exist. The script will not create it for you. This is to prevent the wrong folder being monitored for files which would disrupt read until.", dest='watchdir')
    parser.add('-log', '--log-file', type=str, dest='logfile', default='readuntil.log', help="The name of the log file that data will be written to regarding the decision made by this program to process read until.")
    parser.add('-v', '--verbose-detail', action='store_true', help="Print more detailed coverage info.", default=False, dest='verbosegen')
    parser.add('-v2', '--verbose-true', action='store_true', help="Print detailed messages while processing files.", default=False, dest='verbose')
    parser.add('-sim', '--sim',action='store_true', help="This action will write artificial fast5 files to a folder for testing purposes.",default=False,dest='simcheck')
    parser.add_argument('-ver', '--version', action='version',version=('%(prog)s version={version} date={date}').format(version=__version__,date=__date__))
    #Weak attempt to solve sharing arguments around - really should be passing these.
    global args
    args = parser.parse_args()

    check_files((args.fasta,args.temp_model))
    checkfasta(args.fasta)




    #Checking to see if depths have been set for the amplicons.
    if (len(args.customdepth) > 0 and args.depth > 0):
        print "You can only set depth (-d) or custom depth (-cd), not both! Please try again!"
        #blinker.setjob("warning","red")
        #time.sleep(2)
        die_nicely(oper)
    if (len(args.customdepth) < 1 and args.depth < 1):
        print "You must set either depth (-d) or custom depth (-cd). Please try again!"
        #blinker.setjob("warning","red")
        #time.sleep(2)
        die_nicely(oper)

    #Check to see if folder exists.

    if not os.path.isdir(args.watchdir):
        print "**! Sorry, but the folder "+args.watchdir+" cannot be found.\n\n**!  Please check you have entered the path correctly and try again.\n\n**!  This script will now terminate.\n"
        die_nicely(oper)


    #global p

    p = multiprocessing.Pool(args.procs)


    #Communicating with MinKNOW
    if not args.inhibit:
        if (len(args.customdepth) > 0):
            if args.stoprun:
                message = "minoTour software is implementing read until and will send a stop sequencing message when the run is complete (defined by "+str(args.customdepth)+"x coverage on each amplicon) AND SHUT DOWN minKNOW!!!!"
            else:
                message = "minoTour software is implementing read until and will alert you when the run is complete (defined by "+str(args.customdepth)+"x coverage on each amplicon)."
            startmessagenew = '{"id":"1", "method":"user_message","params":{"content":"minoTour software is implementing read until on this version of minKNOW and will send a stop sequencing message when the run is complete (defined by '+str(args.customdepth)+'x coverage on each amplicon) AND SHUT DOWN minKNOW!!!! If you don\'t want this to happen, cancel your sequencing run or your internet connection now!"}}'
        else:
            if args.stoprun:
                message = "minoTour software is implementing read until and will send a stop sequencing message when the run is complete (defined by "+str(args.depth)+"x coverage on each amplicon) AND SHUT DOWN minKNOW!!!!"
            else:
                message = "minoTour software is implementing read until and will alert you when the run is complete (defined by "+str(args.depth)+"x coverage on each amplicon)."
            startmessagenew = '{"id":"1", "method":"user_message","params":{"content":"minoTour software is implementing read until on this version of minKNOW and will send a stop sequencing message when the run is complete (defined by '+str(args.depth)+'x coverage on each amplicon) AND SHUT DOWN minKNOW!!!! If you don\'t want this to happen, cancel your sequencing run or your internet connection now!"}}'
    else:
        if (len(args.customdepth) > 0):
            message = "minoTour software is monitoring read until and will alert you when you have reached your target of "+str(args.customdepth)+"x coverage on each amplicon."
            startmessagenew = '{"id":"1", "method":"user_message","params":{"content":"minoTour software is monitoring read until on this version of minKNOW and will alert you when you have reached your target of '+str(args.customdepth)+'x coverage on each amplicon."}}'
        else:
            message = "minoTour software is monitoring read until and will alert you when you have reached your target of "+str(args.depth)+"x coverage on each amplicon."
            startmessagenew = '{"id":"1", "method":"user_message","params":{"content":"minoTour software is monitoring read until on this version of minKNOW and will alert you when you have reached your target of '+str(args.depth)+'x coverage on each amplicon."}}'
    print message
    results = execute_command_as_string(startmessagenew,args.ip,8000)


    readuntilrange = 800

    amplicon_file = open(args.ids, "r")
    amplicons = []
    for line in amplicon_file.readlines():
        amplicons.append(line.rstrip())
    if (args.verbose is True):
        print amplicons
    amplicon_file.close()

    ##Checking that amplicon ids are present within the reference sequence.
    idlist=list()
    for record in SeqIO.parse(args.fasta, 'fasta'):
        print record.id
        idlist.append(record.id)

    idsconcat = " ".join(idlist)
    for amplicon in amplicons:
        if amplicon.split(':')[0] not in idsconcat:
            print "!** At least one amplicon is not in your reference sequence.\n\r Please check:"
            print amplicon
            print "!** This program will now exit.\n"
            sys.exit()

    #Creating a dictionary of depths to test against for each amplicon
    global customdepthslist
    customdepthslist=list()
    if (len(args.customdepth) > 0):
        customdepthslist = args.customdepth.split(',')
    else:
        for amplicon in amplicons:
            customdepthslist.append(0)
    fasta_file = args.fasta #fasta file
    model_file_template = args.temp_model #model file
    model_kmer_means_template,kmer_len=process_model_file(model_file_template) #getting lookup dict of kmers and means
    #global kmerhashT #declaring globally a variable to store the kmerhash representing the refernce
    #kmerhashT = process_ref_fasta(fasta_file,model_kmer_means_template,kmer_len) #assign the reference to the hash
    #global seqids,threedarray
    global seqlengths #global declaration of sequence lengths
    seqlengths = get_seq_len(fasta_file) # setting seqeunce lengths
    print seqlengths
    #get_amplicons(amplicons,seqlengths,args) #this seems pointless
    ampdict=[]
    #global ampstartdict #holds the start positions of each amplicon
    ampstartdict=dict()
    #global ampenddict # holds the end position of eeach amplicon
    ampenddict=dict()
    correctedampdict=[] #a dictionary holding corrected amplicon positions for a truncated reference
    #global correctedampstartdict
    correctedampstartdict=dict() #a dictionary holding corrected amplicon start positions
    #global correctedampenddict
    correctedampenddict=dict() #a dictionary holding corrected amplicon end positions
    counter = 0
    #A dict of seen sequences
    global procampres
    procampres=manager.dict()
    #A dict of rejected sequences
    global procampresrej
    procampresrej=manager.dict()
    #A dict of sequences where we have seen the read before
    global procampresseq
    procampresseq=manager.dict()
    #A dict of sequences where reads have been seen before
    global procampresdone
    procampresdone=manager.dict()
    #A dict of files that have been seen
    global procampresfile
    procampresfile=manager.dict()
    #A dict of when files were seen
    global procamprestime
    procamprestime=manager.dict()
    #A dict of 2D files seen
    global procampres2d
    procampres2d=manager.dict()
    #A dict of currently sequencing things
    global currentseq
    currentseq=manager.dict()
    #This dict holds a lookup table for events from read until that have been seen.
    global eventdict
    eventdict=manager.dict()
    #This dict will hold a list of channels, the amplicon type and the start time
    global channeldict
    channeldict=manager.dict()
    meanlengthcount=0
    amplength=0
    for amplicon in amplicons:
        counter+=1
        meanlengthcount+=1
        ampstart = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
        ampstop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
        amplength += (ampstop-ampstart+1)
        ampstartdict[counter]=ampstart
        ampenddict[counter]=ampstop
        ampdict.append((counter,ampstart,ampstop))
        procampres[counter]=0
        procampresrej[counter]=0
        procampresseq[counter]=0
        procampresdone[counter]=0
        procampresfile[counter]=0
        procampres2d[counter]=0
    global readtimeout
    readtimeout = (amplength/meanlengthcount)*5/args.speed
    print "Mean amplicon length:",amplength/meanlengthcount
    global ampliconmeanlength
    ampliconmeanlength = amplength/meanlengthcount
    print "Autocalculated time threshold(s):",readtimeout
    if (args.verbose is True):
        print "******AMP DICTIONARY*******"
        print type(ampstartdict)
        print ampstartdict
    readprediction=dict()
    if (args.verbose is True):
        print procampres
    if (args.verbose is True):
        print "We want to build a custom reference that is smaller than the original reference."
        print "First get a list of all the positions we will need to search"
    keypositions=dict()
    for amplicon in amplicons:
        if (amplicon.split(':',1)[0] not in keypositions):
            keypositions[amplicon.split(':',1)[0]]=[]
        ampstart = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
        ampstop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
        keypositions[amplicon.split(':',1)[0]].append(ampstart)
        keypositions[amplicon.split(':',1)[0]].append(ampstop)
    for sequence in keypositions:
        keypositions[sequence].sort()
    ranges=dict()
    for sequence in keypositions:
        if (sequence not in ranges):
            ranges[sequence]=[]
        for position in keypositions[sequence]:
            start  = position - readuntilrange
            if start < 1:
                start = 1
            stop = position + readuntilrange
            ranges[sequence].append((start,stop))
    for sequence in ranges:
        ranges[sequence] = list(merge_ranges((ranges[sequence])))
    correctedamplicons=list()
    counter=0
    for amplicon in amplicons:
        counter+=1
        #print counter
        ampstart = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
        ampstop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
        correctedamplicons.append(amplicon.split(':',1)[0]+':'+str(correctposition(ampstart,ranges,amplicon.split(':',1)[0]))+'-'+str(correctposition(ampstop,ranges,amplicon.split(':',1)[0])))
        correctedampstartdict[counter]=correctposition(ampstart,ranges,amplicon.split(':',1)[0])
        correctedampenddict[counter]=correctposition(ampstop,ranges,amplicon.split(':',1)[0])
    #global kmerhashTRU
    #customfasta = get_custom_fasta(fasta_file,ranges,args)
    if args.custom is True:
        try:
            seqids,threedarray = get_custom_fasta(fasta_file,ranges,args,model_kmer_means_template,kmer_len) #assign the reference to the hash
        except:
            print "\n\r!** There is a problem in generating your custom reference. This is typically because your seqeunce IDs file is using a different sequence ID to that in your reference.\n\r If your reference sequence fasta is\n\r\n\r >J02459\n\rATGTGA...\n\r\n\r then your amplicon IDS should be in the format:\n\rJ02459:52-1980\n\r"
            print "This script will now exit.\n\r"
            sys.exit()
        ends,starts = correctedampenddict,correctedampstartdict
    else:
        seqids,threedarray = process_ref_fasta(fasta_file,model_kmer_means_template,kmer_len) #assign the reference to the hash
        ends,starts = ampenddict,ampstartdict

    #global seqidscust,threedarraycust
    global current_time
    current_time = time.time()
    if (args.verbose is True):
        print current_time

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
    #blinker.setjob("blink","green")
    if args.inhibit:
        print "This script will not implement read until. It will just report whatever is happening via the read unitl API."
        if not query_yes_no("Are you happy to proceed?"):
            #blinker.setjob("warning","red")
            time.sleep(2)
            #blinker.setjob("lightoff","black")
            time.sleep(1)
            sys.exit()
    else:
        print "This script WILL implement read until.\nIf you proceed it is at your own risk.\n\n"
        if not query_yes_no("Seriously - are you happy to proceed? Entering yes will make it your fault..."):
            #blinker.setjob("warning","red")
            time.sleep(2)
            #blinker.setjob("lightoff","black")
            time.sleep(1)
            die_nicely(oper)
    sample_id = '{"id":"1","method":"get_engine_state","params":{"state_id":"sample_id"}}'
    sampleid = execute_command_as_string(sample_id,args.ip,8000)
    print "Sample Name is: ",sampleid["result"]
    logging.basicConfig(format='%(levelname)s:%(message)s',filename=args.logfile, filemode='w', level=logging.INFO    )
    current_time = time.time()
    print current_time

    print starts
    print ends
    analyser = MyAnalyser(args,p,seqids,threedarray,ends,starts)

    while 1:
        try:
            try:
                print "Running Analysis"
                #blinker.setjob("single","dimgrey")
                run_analysis(args,procampres2d,customdepthslist)
            except socket_error as serr:
                if serr.errno != errno.ECONNREFUSED:
                    raise serr
                    #blinkwarn(bstick)
                    #blinker.setjob("warning","red")
            #blinker.setjob("solid","brown")
            print "Hanging around, waiting for the server..."
            time.sleep(5) # Wait a bit and try again
        except (KeyboardInterrupt, SystemExit):
            print "caught ctrl-c in run loop"
            p.close()
            p.terminate()
            die_nicely(oper)
