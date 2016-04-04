#from read_until import ReadUntil
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
import subprocess
import re
import logging
import glob
import h5py
import platform

from ruutils import process_model_file

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

parser = configargparse.ArgParser(description='real_read_until: A program providing read until with the Oxford Nanopore minION device. This program will ultimately be driven by minoTour to enable selective remote sequencing. This program is heavily based on original code generously provided by Oxford Nanopore Technologies.')
parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file describing the reference sequence for your organism.")
parser.add('-targets', nargs = '*', dest='targets',required=True, help = 'Positional IDs to enrich for in the form seqid:start-stop . Can be space seperated eg: J02459:10000-15000 J02459:35000-40000')
parser.add('-procs', '--proc_num', type=int, dest='procs',required=True, help = 'The number of processors to run this on.')
parser.add('-m', '--model',type=str, required=True, help = 'The appropriate template model file to use', dest='temp_model')
parser.add('-log', '--log-file', type=str, dest='logfile', default='readuntil.log', help="The name of the log file that data will be written to regarding the decision made by this program to process read until.")
parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
parser.add('-o', '--output', type=str, required=False, default='test_read_until_out', help="Path to a folder to symbolically place reads representing match and not match.", dest='output_folder')
parser.add('-v', '--verbose-true', action='store_true', help="Print detailed messages while processing files.", default=False, dest='verbose')
parser.add_argument('-ver', '--version', action='version',version=('%(prog)s version={version} date={date}').format(version=__version__,date=__date__))
args = parser.parse_args()

######################################################
def get_seq_len(ref_fasta):
	seqlens=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		seq=record.seq
		seqlens[record.id]=len(seq)
	return seqlens

"""######################################################
def process_model_file(model_file):
    model_kmers = dict()
    with open(model_file, 'rb') as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        d = list(reader)
        header = d[0]
        meancolumnid = header.index("level_mean")
        kmercolumnid = header.index("kmer")
        d2 = np.array(d[1:])
        model_kmers=dict(d2[:,[0,2]])
        kmer_length=len(d2[0][0])
    return model_kmers,kmer_length
"""
######################################################
def process_ref_fasta2(ref_fasta,model_kmer_means):
	print "processing the reference fasta."
	kmer_len=5
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

def process_ref_fasta_subset(ref_fasta,model_kmer_means,seqlen,kmerlen):
	print "processing the reference fasta."
	kmer_means=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		chunkcounter=0
		for sequence in args.targets:
			print sequence
			chunkcounter += 1
			print "ID", record.id
			print "length", len(record.seq)
			start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
			stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
			seqname = sequence.split(':', 1)[0]
			#length = seqlen[seqid]
			print chunkcounter,sequence,start,stop,seqname
			if seqname in record.id:
				seqchunkname=record.id + "_" + str(chunkcounter)
				print "We want to extract this chunk " + seqchunkname
				kmer_means[seqchunkname]=dict()
				kmer_means[seqchunkname]["F"]=list()
				kmer_means[seqchunkname]["R"]=list()
				kmer_means[seqchunkname]["Fprime"]=list()
				kmer_means[seqchunkname]["Rprime"]=list()
				print "ID", seqchunkname
				print "length", len(record.seq[start:stop])
				print "FORWARD STRAND"
				seq = record.seq[start:stop]
				#print seq
				for x in range(len(seq)+1-kmer_len):
					#print x
					kmer = str(seq[x:x+kmer_len])
					kmer_means[seqchunkname]["F"].append(float(model_kmer_means[kmer]))
				print "REVERSE STRAND"
				seq = revcomp = record.seq[start:stop]
				for x in range(len(seq)+1-kmer_len):
					kmer = str(seq[x:x+kmer_len])
					kmer_means[seqchunkname]["R"].append(float(model_kmer_means[kmer]))

				kmer_means[seqchunkname]["Fprime"]=sklearn.preprocessing.scale(kmer_means[seqchunkname]["F"], axis=0, with_mean=True, with_std=True, copy=True)
				kmer_means[seqchunkname]["Rprime"]=sklearn.preprocessing.scale(kmer_means[seqchunkname]["R"], axis=0, with_mean=True, with_std=True, copy=True)
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
				print "Location",location
			if "Distance" in line:
				distance = float(line.split(': ',1)[1].rstrip('\n'))
				print "Distance",distance
		result.append((distance,id,"F",location))
		os.remove(queryfile)



	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3]


#######################################################################
def squiggle_search2(squiggle,kmerhash,seqlen):
	result=[]

	for ref in kmerhash:
		#print "ss2",ref
		queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)

		dist, cost, path = mlpy.dtw_subsequence(queryarray,kmerhash[ref]['Fprime'])
		result.append((dist,ref,"F",path[1][0],ref,path[1][-1]))
		dist, cost, path = mlpy.dtw_subsequence(queryarray,kmerhash[ref]['Rprime'])
		result.append((dist,ref,"R",path[1][0],ref,path[1][-1]))


	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4],sorted(result,key=lambda result: result[0])[0][5]

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
def go_or_no(seqid,direction,position,seqlen,args):
	for sequence in args.targets:
		#print sequence
		start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
		stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
		length = seqlen[seqid]
		#We note that the average template read length is 6kb for the test lambda dataset. Therefore we are interested in reads which start at least 3kb in advance of our position of interest
		balance = 3000
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

def mp_worker((channel_id, data,kmerhash,seqlen,readstarttime,kmerhash_subset)):
	#for ref in kmerhash:
		#print ref
	#print "worker running"
	#print "channel_id:",channel_id
	#print "readstarttime:",readstarttime
	#for d in data:
	 #   print d
	#print type(data.read_id)
	if ((time.time()-readstarttime) > args.time):
		print "We have a timeout"
		#logging.info('%s,%s,%s,%s', channel_id, data.read_id, 'TOT',data.events[0].start)
		return 'timeout',channel_id,data.read_id,data.events[0].start
	#elif channel_id %2 == 0:
	#	print "Even numbered channel so skip"
	#	return 'evenskip',channel_id,data.read_id,data.events[0].start
	else:
		try:
			print "Read start time",readstarttime
			#print "Elapsed time since read=",(time.time()-readstarttime)
			squiggle = extractsquig(data.events)
			#print data.events[0].start
			#result = 'bernard'
			squiggleres = squiggle_search2(squiggle,channel_id,data.read_id,kmerhash,seqlen)
			print "Full Length:",squiggleres
			print "Full Length Match Length:", squiggleres[5]-squiggleres[3]
			squiggleres2 = squiggle_search2(squiggle,channel_id,data.read_id,kmerhash_subset,seqlen)
			squiggleres3 = squiggle_search2(squiggle[100:200],channel_id,data.read_id,kmerhash_subset,seqlen)
			squiggleres4 = squiggle_search2(squiggle[50:100],channel_id,data.read_id,kmerhash_subset,seqlen)

			if squiggleres2[5] > squiggleres3[3] > squiggleres2[3] and squiggleres3[3] > squiggleres4[3] > squiggleres2[3]:
				print "!!!!!!!!!!!!!!!! We got a good one! !!!!!!!!!!!!!!!!"
				print "Subset:",squiggleres2
				print "Subset Match Length:", squiggleres2[5]-squiggleres2[3]
				print "SecondHalf:",squiggleres3
				print "SecondHalf:", squiggleres3[5]-squiggleres3[3]
				print "FirstHalf:", squiggleres4, squiggleres4[5]-squiggleres4[3]
				result = "Sequence"
			else:
				result = "Skip"
				#print "This read don't match"
				#print "Subset:",squiggleres2
				#print "Subset Match Length:", squiggleres2[5]-squiggleres2[3]

			#result = go_or_no(squiggleres[0],squiggleres[2],squiggleres[3],seqlen)
			#print "result:",result
			return result,channel_id,data.read_id,data.events[0].start,squiggleres
		except Exception, err:
			err_string="Time Warping Stuff : %s" % ( err)
			print >>sys.stderr, err_string
####################






def process_hdf5((filename,kmerhash_subset,procampres,seqlen,args)):
    #returnlist=list()
	#print filename
	returndict=dict()
	#print "Can we get in to a read here?"
	hdf = h5py.File(filename, 'r')
	for read in hdf['Analyses']['EventDetection_000']['Reads']:
		#print "We're in a read"
		events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
		event_collection=list()
		for event in events:
			event_collection.append(float(event[0]))
		#print event_collection[50:350]
		#print event_collection

		squiggle = event_collection[50:350]

		#squiggleres = squiggle_search2(squiggle,channel_id,data.read_id,kmerhash,seqlen)
        #print squiggleres
		#result = go_or_no(squiggleres[0],squiggleres[2],squiggleres[3],seqlen)
		#print result

		#print data.events[0].start
		#result = 'bernard'
		#print squiggle
		#exit()
		squiggleres = squiggle_search2(squiggle,kmerhash_subset,len(squiggle))
	##	print "Full Length:",squiggleres
	#	print "Full Length Match Length:", squiggleres[5]-squiggleres[3]
		squiggleres2 = squiggle_search2(squiggle,kmerhash_subset,len(squiggle))
		#print "Do we get here?"
		squiggleres3 = squiggle_search2(squiggle[150:300],kmerhash_subset,len(squiggle[150:300]))
		#squiggleres4 = squiggle_search2(squiggle[25:125],kmerhash_subset,len(squiggle[25:125]))
		#if squiggleres2[5] > squiggleres3[3] > squiggleres2[3] and squiggleres3[3] > squiggleres4[3] > squiggleres2[3]:
		if squiggleres2[5] > squiggleres3[3] > squiggleres2[3]:
	#		print "!!!!!!!!!!!!!!!! We got a good one! !!!!!!!!!!!!!!!!"
	#		print "Subset:",squiggleres2
	#		print "Subset Match Length:", squiggleres2[5]-squiggleres2[3]
	#		print "SecondHalf:",squiggleres3
	#		print "SecondHalf:", squiggleres3[5]-squiggleres3[3]
	#		print "FirstHalf:", squiggleres4, squiggleres4[5]-squiggleres4[3]
			print squiggleres[0],squiggleres[2],squiggleres[3],seqlen
			try:
				result = go_or_no(squiggleres[0],squiggleres[2],squiggleres[3],seqlen,args)
				print "result", result
			except Exception,err:
				print err

			#print result
			result = "Sequence"
		else:
			#print "Subset:",squiggleres2
			#print "Subset Match Length:", squiggleres2[5]-squiggleres2[3]
			result = "Skip"

	hdf.close()
	return (result,filename)


def mycallback((result,filename)):
	#print "done"
	print filename,
	filetocheck = os.path.split(filename)
	sourcefile = filename
	if result == "Sequence":
		path1 = os.path.join(args.output_folder,'sequence')
		path2 = os.path.join(path1,'downloads')
		path3 = os.path.join(path2,'pass')
		path4 = os.path.join(path2,'fail')
		if not os.path.exists(path1):
			os.makedirs(path1)
		if not os.path.exists(path2):
			os.makedirs(path2)
		if not os.path.exists(path3):
			os.makedirs(path3)
		if not os.path.exists(path4):
			os.makedirs(path4)

		print "Sequence Found"
		if "pass" in filename:
			destfile = os.path.join(path3,filetocheck[1])
		else:
			destfile = os.path.join(path4,filetocheck[1])
		try:
			#os.symlink(sourcefile, destfile)
			shutil.copy(sourcefile,destfile)
		except Exception, err:
			print "File Copy Failed",err
	else:
		path1 = os.path.join(args.output_folder,'reject')
		path2 = os.path.join(path1,'downloads')
		path3 = os.path.join(path2,'pass')
		path4 = os.path.join(path2,'fail')
		if not os.path.exists(path1):
			os.makedirs(path1)
		if not os.path.exists(path2):
			os.makedirs(path2)
		if not os.path.exists(path3):
			os.makedirs(path3)
		if not os.path.exists(path4):
			os.makedirs(path4)
		print "No Match"
		if "pass" in filename:
			destfile = os.path.join(path3,filetocheck[1])
		else:
			destfile = os.path.join(path4,filetocheck[1])
		try:
			#os.symlink(sourcefile, destfile)
			shutil.copy(sourcefile,destfile)
		except Exception, err:
			print "File Copy Failed",err

if __name__ == "__main__":
	print "***********************************************************************************************"
	print "**** This code will open a collection of reads and simulate read until on them. It will    ****"
	print "**** copy reads into a secondary folder for subsequent processing by another analysis      ****"
	print "**** package.                                                                              ****"
	print "***********************************************************************************************"

	p = multiprocessing.Pool(args.procs)
	manager = multiprocessing.Manager()
	procampres=manager.dict()
	fasta_file = args.fasta
	seqlen = get_seq_len(fasta_file)

	model_file = args.temp_model
	model_kmer_means,kmer_len=process_model_file(model_file)

	kmerhash_subset = process_ref_fasta_subset(fasta_file,model_kmer_means,seqlen,kmer_len)

	d=list()
	filenamecounter=0
	for filename in glob.glob(os.path.join(args.watchdir, '*.fast5')):
		filenamecounter+=1
		#print filename
		d.append([filename,kmerhash_subset,procampres,seqlen,args])
	for filename in glob.glob(os.path.join(args.watchdir, "pass",'*.fast5')):
		filenamecounter+=1
		#print filename
		d.append([filename,kmerhash_subset,procampres,seqlen,args])
	for filename in glob.glob(os.path.join(args.watchdir, "fail",'*.fast5')):
		filenamecounter+=1
		#print filename
		d.append([filename,kmerhash_subset,procampres,seqlen,args])
	procdata=tuple(d)

	results=[]
	for x in (procdata):
		r = p.apply_async(process_hdf5, (x,), callback=mycallback)
		results.append(r)
	for r in results:
		r.wait()
