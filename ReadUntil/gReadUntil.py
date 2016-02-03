from read_until import ReadUntil
import time
import errno
from socket import error as socket_error
import threading
import MySQLdb
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
#import _ucrdtw
#from fastdtw import fastdtw, dtw

parser = configargparse.ArgParser(description='gReadUntil.py: A program providing read until for genome sequences with the Oxford Nanopore minION device. This program will ultimately be driven by minoTour to enable selective remote sequencing. This program is partly based on original code generously provided by Oxford Nanopore Technologies.')
parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file describing the reference sequence for your organism.")
parser.add('-ids', nargs = '*', dest='ids',required=True, help = 'Positional IDs to enrich for in the form seqid:start-stop . Can be space seperated eg: J02459:10000-15000 J02459:35000-40000')
parser.add('-procs', '--proc_num', type=int, dest='procs',required=True, help = 'The number of processors to run this on.')
parser.add('-t', '--time', type=int, dest='time', required=True, default=300, help="This is an error catch for when we cannot keep up with the rate of sequencing on the device. It takes a finite amount of time to process through the all the channels from the sequencer. If we cannot process through the array quickly enough then we will \'fall behind\' and lose the ability to filter sequences. Rather than do that we set a threshold after which we allow the sequencing to complete naturally. The default is 300 seconds which equates to 9kb of sequencing at the standard rate.")
parser.add('-m', '--model',type=str, required=True, help = 'The appropriate template model file to use', dest='temp_model')
parser.add('-l', '--model_length',type=int, required=True, help = 'The word size of the mode file - e.g 5,6 or 7', dest='model_length')
parser.add('-i', '--index', type=int, dest ='indexpos', default=1, required=False, help = 'The index position of mean events in the reference model file.')
parser.add('-ip', '--ip-address', type=str ,dest='ip',required=False,default="127.0.0.1", help="The IP address of the machine running minKNOW.")
parser.add('-p', '--port', type=int, dest='port', default=None,required=True, help='The port that ws_event_sampler is running on.' )
parser.add('-log', '--log-file', type=str, dest='logfile', default='readuntil.log', help="The name of the log file that data will be written to regarding the decision made by this program to process read until.")
args = parser.parse_args()

######################################################
def get_seq_len(ref_fasta):
	seqlens=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		seq=record.seq
		seqlens[record.id]=len(seq)
	return seqlens



######################################################
def process_model_file(model_file):
	model_kmers = dict()
	with open(model_file, 'rb') as csv_file:
		reader = csv.reader(csv_file, delimiter="\t")
    		d = list(reader)
		for r in range(1, len(d)):
			kmer = d[r][0]
			mean = d[r][args.indexpos]
			print mean
			if (float(mean) <= 25):
				print "I'm almost certain you are not looking at means here - you need to fix this!"
				exit()
			model_kmers[kmer]=mean
	return 	model_kmers

######################################################
def process_ref_fasta2(ref_fasta,model_kmer_means):
	print "processing the reference fasta."
	kmer_len=args.model_length
	kmer_means=dict()

	for record in SeqIO.parse(ref_fasta, 'fasta'):
		kmer_means[record.id]=dict()
		kmer_means[record.id]["F"]=list()
#		kmer_means[record.id]["R"]=list()

		seq = record.seq
		for x in range(len(seq)+1-kmer_len):
			kmer = str(seq[x:x+kmer_len])
			kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))

#		seq = revcomp = record.seq.reverse_complement()
#		for x in range(len(seq)+1-kmer_len):
#			kmer = str(seq[x:x+kmer_len])
#			kmer_means[record.id]["R"].append(float(model_kmer_means[kmer]))

	return kmer_means
#######################################################################

def process_ref_fasta(ref_fasta,model_kmer_means):
	print "processing the reference fasta."
	kmer_len=args.model_length
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


def runProcess(exe):
	p=subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	while(True):
		retcode= p.poll()
		line=p.stdout.readline()
		yield line
		if(retcode is not None):
			break

#######################################################################
def squiggle_search(squiggle,kmerhash,channel_id,read_id,seqlen):
	result=[]
	for id in kmerhash:
#		query = sklearn.preprocessing.scale(squiggle,axis=0,with_mean=True,with_std=True,copy=True)
#		compa = sklearn.preprocessing.scale(kmerhash[id]["F"],axis=0,with_mean=True,with_std=True,copy=True)
#		compb = sklearn.preprocessing.scale(kmerhash[id]["R"],axis=0,with_mean=True,with_std=True,copy=True)
		#starttime = time.time()
#		dist, cost, path = mlpy.dtw_subsequence(query,compa)
#		result.append((dist,id,"F",path))
#		dist, cost, path = mlpy.dtw_subsequence(query,compb)
#		result.append((dist,id,"R",path))
		#Here we are going to try to call a gpu based time warp for this data. To do this we need to write out a file to query with
		#Ideally this should have a unique name. We shall call it channel_id_read_id.
		#To do this we need the read_id and channel_id
		queryfile=str(channel_id)+"_"+str(read_id)+"_query.bin"
		#We are going to normalise this sequence with the sklearn preprocessing algorithm to see what happens.
		queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
		with open(queryfile, "wb") as f:
			f.write(ar.array("f", queryarray))
		subjectfile = id+"_"+"F"+"_subject.bin"
		subjectfile = re.sub('\|','_',subjectfile)
		seqlen2 = str(seqlen[id])
		commands = queryfile+' '+subjectfile+' 200 '+seqlen2+' 0.05'
		current = str(multiprocessing.current_process())
		currentnum=int(re.search(r'\d+', current).group())
		gpucode=str()
		if (currentnum % 2 == 0):
			#print "Even"
			gpucode='./GPU-DTW '
		else:
			#print "Odd"
			gpucode='./GPU-DTW '
		#print "Running forward";
		runcommand = gpucode+commands
		location = ()
		distance = ()
		for line in runProcess(runcommand.split()):
			#print line.rstrip('\n')
			if "Location" in line:
				location = int(line.split(': ',1)[1].rstrip('\n'))
		#		print "Location",location
			if "Distance" in line:
				distance = float(line.split(': ',1)[1].rstrip('\n'))
		#		print "Distance",distance
		result.append((distance,id,"F",location))
#		subjectfile2 = id+"_"+"R"+"_subject.bin"
#		subjectfile2 = re.sub('\|','_',subjectfile2)
#		seqlen2 = str(seqlen[id])
#		commands = queryfile+' '+subjectfile2+' 200 '+seqlen2+' 0.05'
#		#print "Running Reverse"
#		runcommand = gpucode+commands
#		location = ()
#		distance = ()
#		for line in runProcess(runcommand.split()):
#			#print line.rstrip('\n')
#			if "Location" in line:
#				location = int(line.split(': ',1)[1].rstrip('\n'))
#		#		print "Location",location
#			if "Distance" in line:
#				distance = float(line.split(': ',1)[1].rstrip('\n'))
#		#		print "Distance",distance
#		result.append((distance,id,"R",location))
		os.remove(queryfile)



	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3]


#######################################################################
def squiggle_search2(squiggle,channel_id,read_id,kmerhash,seqlen):
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
def go_or_no(seqid,direction,position,seqlen):
	for sequence in args.ids:
		start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
		stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
		length = seqlen[seqid]
		#We note that the average template read length is 6kb for the test lambda dataset. Therefore we are interested in reads which start at least 3kb in advance of our position of interest
		balance = 1000
		if seqid.find(sequence.split(':', 1 )[0]) > 0:
			#print "Found it"
			if direction == "F":
				#print "Forward Strand"
				if position >= ( start - balance ) and position <= stop:
					return "Sequence"
			elif direction == "R":
				if position >= ( length - stop - balance) and position <= ( length - start ):
					#print "Reverse Strand"
					return "Sequence"
	return "Skip"

###################
def mp_worker((channel_id, data,kmerhash,seqlen,readstarttime)):
    #for ref in kmerhash:
    	#print ref
    #print "worker running"
    #print "channel_id:",channel_id
    #print "readstarttime:",readstarttime
    #for d in data:
     #   print d
    #print type(data.read_id)
    if ((time.time()-readstarttime) > args.time):
    	#print "We have a timeout"
        #logging.info('%s,%s,%s,%s', channel_id, data.read_id, 'TOT',data.events[0].start)
        return 'timeout',channel_id,data.read_id,data.events[0].start
    elif channel_id %2 == 0:
    #elif channel_id >= 252:
        #print "Even numbered channel so skip"
        return 'evenskip',channel_id,data.read_id,data.events[0].start
    else:
        #print "Odd numbered channel - yay!"
        try:
            #print "Read start time",readstarttime
            #print "Elapsed time since read=",(time.time()-readstarttime)
            squiggle = extractsquig(data.events)
            #print data.events[0].start
            #result = 'bernard'
            squiggleres = squiggle_search2(squiggle,channel_id,data.read_id,kmerhash,seqlen)
            #print squiggleres
            result = go_or_no(squiggleres[0],squiggleres[2],squiggleres[3],seqlen)
            #print "result:",result
            return result,channel_id,data.read_id,data.events[0].start,squiggleres
        except Exception, err:
            err_string="Time Warping Stuff : %s" % ( err)
            print >>sys.stderr, err_string
####################




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
    model_kmer_means=process_model_file(model_file)
    #model_kmer_means = retrieve_model()
    #global kmerhash
    kmerhash = process_ref_fasta(fasta_file,model_kmer_means)

    print type (kmerhash)


    def __init__(self):
        self.current_unblock_map = LockedDict()
        #self.p = multiprocessing.Pool(48)
        #for id in kmerhash:
		#for ref in kmerhash[id]:
		#	print id,ref

    def mycallback(self, actions):
    	#print "Callback CALLED"
    	#print type(actions)
    	#print actions[0]
        #for action in actions2:
        #print action[3][0]
        if actions[0] == "Skip":
            #print "Should skip"
            self.current_unblock_map[actions[1]]=actions[2]
            logging.info('%s,%s,%s,%s,%s,%s,%s,%s', actions[1], actions[2], 'REJ',actions[3],actions[4][0],actions[4][1],actions[4][2],actions[4][3])
        elif actions[0] == "timeout":
            print "Read timeout"
            logging.info('%s,%s,%s,%s', actions[1], actions[2], 'TOT',actions[3])
            #logging.info('%s,%s,%s', actions[1],actions[2],"TOT")
        elif actions[0] == "evenskip":
            #print "Even Skip"
            logging.info('%s,%s,%s,%s', actions[1], actions[2], 'EVE',actions[3])
        else:
            #print "Sequencing from ",action[1],action[2]
            #print actions[1], actions[2], 'SEQ',actions[3],actions[4][0],actions[4][1],actions[4][2],actions[4][3]
            logging.info('%s,%s,%s,%s,%s,%s,%s,%s', actions[1], actions[2], 'SEQ',actions[3],actions[4][0],actions[4][1],actions[4][2],actions[4][3])
        #If we are building up channels to unblock we might as well unblock them now. But this doesn't work as we cannot access the myclient code.
        #unblock_now = self.next_unblock_map()
        #if len(unblock_now)>40:
        #    print "Unblocking channels: ", unblock_now.keys()
        #    print time.strftime('%Y-%m-%d %H:%M:%S'),
        #    print "Unblocking ",len(unblock_now.keys())
        #    ReadUntil.unblock(unblock_now)

    def data_received(self, channels):

        print "Channels length",len(channels)
        d=list()
        #print "D length",len(d)
        print "Checking Channels"
        for channel_id, data in channels.iteritems():
            #print data
            d.append([channel_id,data,seqlen,time.time()])
        	#print channel_id,data.read_id,data.events[0].start,time.ctime(current_time + (data.events[0].start/10000)	),
        	#print type(data)
        	#print type(data.events)
        procdata=tuple(d)
        #print "D length",len(d)
        #print "procdata length",len(procdata)
        #p = multiprocessing.Pool(48)
	#for element in d:
		#print type (element)
        results = p.map_async(mp_worker, (procdata),chunksize=1,callback=self.mycallback)
	#p.apply_async(mp_worker, (procdatat), callback=procresult)
        #p.close()
        #p.join()
#        seqcount = 0
#        for result in results:
#            #print result
#            if result[0] == "Skip":
#                #print "Should skip"
#                self.current_unblock_map[result[1]]=result[2]
#            else:
#                seqcount += 1
#        print time.strftime('%Y-%m-%d %H:%M:%S'),
#        print "Sequencing from ",seqcount

    #def procresult(result,channel_id,data.read_id):

    def apply_async_with_callback(self, channels):
        #global kmerhash
        #for ref in kmerhash:
        #   print ref
        print "Channels length",len(channels)
        d=list()
        print "Checking Channels"
        for channel_id, data in channels.iteritems():
            #print data
            #d.append([channel_id,data,kmerhash,seqlen,time.time()])
            MyAnalyser.p.apply_async(mp_worker, args = ((channel_id,data,MyAnalyser.kmerhash,MyAnalyser.seqlen,time.time()), ), callback = self.mycallback)
        #procdata=tuple(d)
        #pool = multiprocessing.Pool(8)
        #for d in procdata:

        #pool.close()
        #pool.join()
        #print(result_list)

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

def run_analysis():
    analyser = MyAnalyser()
    host = "ws://"+str(args.ip)+":"+str(args.port)+"/"
    #host = "ws://127.0.0.1:5999/"
    setup_conditions = {"ignore_first_events": 100, "padding_length_events": 0,
                        "events_length": 250, "repetitions": 1}

    state=RunningState()
    with ReadUntil(host=host,
                   setup_conditions=setup_conditions,
                   data_received=analyser.apply_async_with_callback,
                   connection_closed=state.closed) as my_client:
        # Start sending stuff to our analyser
        my_client.start()
        print "Client connection started. Beginning unblock loop..."
        while state.keep_running:
            unblock_now = analyser.next_unblock_map()
            #print "Checking unblock status"
            if len(unblock_now)>0:
                print "Unblocking channels: ", unblock_now.keys()
                print time.strftime('%Y-%m-%d %H:%M:%S'),
                print "Unblocking ",len(unblock_now.keys())
                my_client.unblock(unblock_now)

            # Throttle rate at which we make unblock controls. Although
            # unblocks should be timely, it is more efficient on the network
            # and on the hardware to unblock a bigger list of channels at once.
            time.sleep(15)
        print "...unblock loop ended. Connection closed."

if __name__ == "__main__":
	#global p
	print "***********************************************************************************************"
	print "**** This version of the code will not process reads derived from even numbered channels.  ****"
	print "**** This is in the vain attempt that we might generate a really cool visual!              ****"
	print "***********************************************************************************************"
	logging.basicConfig(format='%(levelname)s:%(message)s',filename=args.logfile, filemode='w', level=logging.INFO	)
	#logging.debug('This message should go to the log file')
	#logging.info('So should this')
	#logging.warning('And this, too')
	#p = multiprocessing.Pool(4)
    # A few extra bits here to automatically reconnect if the server goes down
    # and is brought back up again.
	current_time = time.time()
	print current_time
	#for id in kmerhash:
	#	for ref in kmerhash[id]:
	#		print id,ref
	#		print type(kmerhash[id][ref])
	#		testarray = sklearn.preprocessing.scale(np.array(kmerhash[id][ref][0:10000]),axis=0,with_mean=True,with_std=True,copy=True)
	#		filename = id+"_"+ref+"_subject.bin"
	#		filename = re.sub('\|','_',filename)
	#		with open(filename, "wb") as f:
	#			f.write(ar.array("f", testarray))
	#		filename = id+"_"+ref+"_subject.txt"
	#		filename = re.sub('\|','_',filename)
	#		np.savetxt(filename, testarray, delimiter=',')
	#		print len(testarray)
	#		filename2 = id+"_"+ref+"_testquery.bin"#
	#		filename2 = re.sub('\|','_',filename2)
	#		with open(filename2, "wb") as f:
	#			f.write(ar.array("f", np.array(kmerhash[id][ref])[700:1000]))
	#		filename2 = id+"_"+ref+"_testquery.txt"
	#		filename2 = re.sub('\|','_',filename2)
	#		np.savetxt(filename2, np.array(kmerhash[id][ref])[700:1000], delimiter=',')#
	#		print len(testarray[700:1000])


	while 1:
		try:
			print "Running Analysis"
			run_analysis()
		except socket_error as serr:
			if serr.errno != errno.ECONNREFUSED:
				raise serr
		print "Hanging around, waiting for the server..."
		time.sleep(5) # Wait a bit and try again
