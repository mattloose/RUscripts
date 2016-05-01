#!C:\anaconda python
import sys, os, re
import time
import threading, thread
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from StringIO import StringIO
import string
import mlpy
#import sklearn.preprocessing
import random
import math
import csv
import numpy as np
import array as ar
import configargparse
import subprocess
import shutil
import glob
import h5py
from itertools import islice
from collections import OrderedDict
import psutil
import multiprocessing
import platform
import bisect
sys.path.append("ReadUntil")
from ruutils import process_model_file,check_files,checkfasta









###########################################################
def make_hdf5_object_attr_hash(hdf5object, fields):
	att_hash=dict()
	for field in fields:
		if (field in hdf5object.attrs.keys() ):
			#print "filed: ",field (args.ref_fasta is not None), hdf5object.attrs[field]
			att_hash[field]=hdf5object.attrs[field]
	return att_hash

######################################################
def scale(a):  # MS
   mu = np.mean(a, None)
   sigma = np.std(a)
   if sigma == 0: return 0
   else: return (a - mu) / sigma


######################################################
def get_amplicons():
    #print "Groking amplicons"
    if (args.verbose is True):
        print "ids is of type", type(amplicons)
    for sequence in amplicons:
        if (args.verbose is True):
            print sequence
        start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
        stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
        seqname = sequence.split(':',1)[0]
        if (args.verbose is True):
            print start
            print stop
        REVERSE_stop = seqlengths[seqname]-start
        REVERSE_start = seqlengths[seqname]-stop
        if (args.verbose is True):
            print REVERSE_stop
            print REVERSE_start

######################################################
def get_seq_len(ref_fasta):
	seqlens=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		seq=record.seq
		seqlens[record.id]=len(seq)
	return seqlens
#######################################################################
def raw_squiggle_search2(squiggle,hashthang):
    result=[]
    #print args.speedmode
    for ref in hashthang:
        try:
            queryarray = scale(squiggle)
            mx = np.max(queryarray)
            scalingFactor = 1 # iqr # 3 # 1.2 # MS
            queryarray *= scalingFactor
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Fprime'])
            result.append((dist,ref,"F",path[1][0],path[1][-1],path[0][0],path[0][-1]))
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Rprime'])
            result.append((dist,ref,"R",(len(hashthang[ref]['Rprime'])-path[1][-1]),(len(hashthang[ref]['Rprime'])-path[1][0]),path[0][0],path[0][-1]))
        except Exception,err:
            print "Warp Fail"
    return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4],sorted(result,key=lambda result: result[0])[0][5],sorted(result,key=lambda result: result[0])[0][6]

def raw_squiggle_search3(squiggle,hashthang):
    result=[]
    for ref in hashthang:
        try:
            queryarray=scale(squiggle)
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Fprimewin'])
            result.append((dist,ref,"F",path[1][0],path[1][-1],path[0][0],path[0][-1],cost.mean()))
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Rprimewin'])
            result.append((dist,ref,"R",(len(hashthang[ref]['Rprimewin'])-path[1][-1]),(len(hashthang[ref]['Rprimewin'])-path[1][0]),path[0][0],path[0][-1],cost.mean()))
        except Exception,err:
            print "Warp Fail"
    return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4],sorted(result,key=lambda result: result[0])[0][5],sorted(result,key=lambda result: result[0])[0][6]


######################################################
def get_custom_fasta(ref_fasta,model_kmer_means,subsectionlist):
    print "Generating a custom fasta"
    sequencedict=dict()
    for sequence in subsectionlist:
        print sequence
        for record in SeqIO.parse(ref_fasta, 'fasta'):
            if (record.id == sequence):
                if (sequence not in sequencedict):
                    sequencedict[sequence]=list()
                for sections in subsectionlist[sequence]:
                    start = sections[0]
                    end = sections[1]
                    #print start,end
                    #print record.seq[sections[0]-1:sections[1]-1]
                    if (len(sequencedict[sequence])>0):
                        sequencedict[sequence]=str(sequencedict[sequence])+str(record.seq[sections[0]-1:sections[1]-1])
                    else:
                        sequencedict[sequence]=str(record.seq[sections[0]-1:sections[1]-1])
                    #print record.seq[int(sections[0]),int(sections[1])]
                #print sequencedict[sequence]
    print "processing the custom fasta"
    kmer_len=len(model_kmer_means.keys()[0])
    kmer_means=dict()
    for sequence in sequencedict:
        #print sequence
        kmer_means[sequence]=dict()
        kmer_means[sequence]["F"]=list()
        kmer_means[sequence]["R"]=list()
        kmer_means[sequence]["Fprime"]=list()
        kmer_means[sequence]["Rprime"]=list()
        seq = Seq(sequencedict[sequence], generic_dna)
        #seq = sequencedict[sequence]
        if (args.verbose is True):
            print "ID", sequence
            print "length", len(seq)
            print "FORWARD STRAND"
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[sequence]["F"].append(float(model_kmer_means[kmer]))
        if (args.verbose is True):
            print "REVERSE STRAND"
        seq = revcomp = seq.reverse_complement()
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[sequence]["R"].append(float(model_kmer_means[kmer]))
        #kmer_means[sequence]["Fprime"]=sklearn.preprocessing.scale(kmer_means[sequence]["F"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[sequence]["Fprime"]=scale(kmer_means[sequence]["F"])
        #kmer_means[sequence]["Rprime"]=sklearn.preprocessing.scale(kmer_means[sequence]["R"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[sequence]["Rprime"]=scale(kmer_means[sequence]["R"])
    return kmer_means

######################################################
def check_basecalled(hdf):
    '''
    Function to check if an hdf file is basecalled.
    '''
    for element in hdf:
        for element2 in hdf[element]:
            for element3 in hdf[element][element2]:
                for element4 in hdf[element][element2][element3]:
                    if any("Model" in s for s in [element,element2,element3,element4]):
                        return True
    return False

######################################################
def process_ref_fasta_raw(ref_fasta,model_kmer_means):
    print "processing the reference fasta."
    kmer_len=len(model_kmer_means.keys()[0])
    kmer_means=dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id]=dict()
        kmer_means[record.id]["F"]=list()
        kmer_means[record.id]["R"]=list()
        kmer_means[record.id]["Fprime"]=list()
        kmer_means[record.id]["Rprime"]=list()
        #if (args.verbose is True):
        print "ID", record.id
        print "length", len(record.seq)
        print "FORWARD STRAND"
        seq = record.seq
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))
            #print kmer, float(model_kmer_means[kmer])
        #if (args.verbose is True):
        print "REVERSE STRAND"
        seq = revcomp = record.seq.reverse_complement()
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["R"].append(float(model_kmer_means[kmer]))
            #print kmer, float(model_kmer_means[kmer])
        #kmer_means[record.id]["Fprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["F"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]["Fprime"]=scale(kmer_means[record.id]["F"])
        #kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]["Rprime"]=scale(kmer_means[record.id]["R"])
    return kmer_means
#######################################################################
def process_hdf5((filename,procampres,kmerhashTRU,correctedamplicons,correctedampstartdict,correctedampenddict)):
        readprediction=dict()
        #print "yo"
        if (args.verbose is True):
            print filename
        #print startpos,endpos
        hdf = h5py.File(filename, 'r')
        for read in hdf['Analyses']['EventDetection_000']['Reads']:
            events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
            event_collection=list()
            for event in events:
                event_collection.append(float(event['mean']))
            if (len(event_collection) >= 300):
                (seqmatchnameR,distanceR,frR,rsR,reR,qsR,qeR) = raw_squiggle_search2(event_collection[50:300],kmerhashTRU)
                if (args.verbose is True):
                    print seqmatchnameR, distanceR, frR, rsR, reR, qsR, qeR

                if distanceR < 40000: #This parameter is an unused test at this point - it is testing on valid distances
                    if (frR == "F"): #This searches for the closest match on the forward strand.
                        amplicon, value = min(correctedampstartdict.items(), key=lambda (_, v): abs(v - rsR))
                    else: #This searches for the closest match on the reverse strand.
                        amplicon, value = min(correctedampenddict.items(), key=lambda (_, v): abs(v - rsR))
                        if (args.verbose is True):
                            print amplicon, value

                    procampres[amplicon] += 1
                    if (amplicon not in readprediction):
                        readprediction[amplicon]=dict()
                    if (0 not in readprediction[amplicon]):
                        readprediction[amplicon][0]=dict()
                    if (filename not in readprediction[amplicon][0]):
                        readprediction[amplicon][0][filename]=dict()
                    readprediction[amplicon][0][filename]["name"]=filename
                    readprediction[amplicon][0][filename]["matchdistance"]=distanceR
                else:
                    print "Poor match"

        hdf.close()
        #print procampres
        #print readprediction
        return readprediction

######################

def merge_ranges(ranges):
    """
    Merge overlapping and adjacent ranges and yield the merged ranges
    in order. The argument must be an iterable of pairs (start, stop).

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop

######################

def correctposition(value,ranges,sequence):
    correction = 0
    for range in ranges[sequence]:
        #print range
        if value >= range[0] and value <= range[1]:
            #print "Found it"
            return value - range[0] + correction + 1
        correction = correction + (range[1]-range[0] + 1)

######################

if __name__ == "__main__":
    multiprocessing.freeze_support()


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

    __version__ = "1.1"
    __date__ = "1st May 2016"

    parser = configargparse.ArgParser(description='ampliconSPLIT: A program designed to identify and group individual amplicons from minION reads prior to base calling. The depth setting limits the number of reads copied to each sub folder. Developed by Matt Loose @mattloose or matt.loose@nottingham.ac.uk for help!',default_config_files=[config_file])
    parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file for the reference sequence for your organism.")
    parser.add('-ids', '--reference_amplicon_positions', type=str, required=True, default=None, help="A file containing a list of amplicon positions defined for the reference sequence. 1 amplicon per line in the format fasta_sequence_name:start-stop e.g EM_079517:27-1938", dest='ids')
    parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
    parser.add('-o', '--output-dir', type=str, required=True, default="prefiltered", help="The path to the destination folder for the preprocessed reads" , dest="targetpath")
    parser.add('-d', '--depth',type=int, required=True, default=None, help = 'The desired coverage depth for each amplicon. Note this is unlikely to be achieved for each amplicon and should probably be an overestimate of the minimum coverage required.', dest='depth')
    parser.add('-procs', '--proc_num', type=int, dest='procs',required=True, help = 'The number of processors to run this on.')
    parser.add('-t', '--template_model',type=str, required=True, help = 'The appropriate template model file to use. This file can be generated uing the getmodels.py script.', dest='temp_model')
    parser.add('-v', '--verbose-true', action='store_true', help="Print detailed messages while processing files.", default=False, dest='verbose')
    parser.add_argument('-ver', '--version', action='version',version=('%(prog)s version={version} date={date}').format(version=__version__,date=__date__))
    args = parser.parse_args()

    check_files((args.fasta,args.temp_model))
    checkfasta(args.fasta)

    if not os.path.isdir(args.watchdir):
        print "**! Sorry, but the folder "+args.watchdir+" cannot be found.\n\n**!  Please check you have entered the path correctly and try again.\n\n**!  This script will now terminate.\n"
        sys.exit()



    p = multiprocessing.Pool(args.procs)
    readuntilrange = 900
    manager = multiprocessing.Manager()



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


    fasta_file = args.fasta
    model_file_template = args.temp_model

    model_kmer_means_template,kmer_len=process_model_file(model_file_template)

    kmerhashT = process_ref_fasta_raw(fasta_file,model_kmer_means_template)

    seqlengths = get_seq_len(fasta_file)
    get_amplicons()

    ampdict=[]
    ampstartdict=dict()
    ampenddict=dict()
    correctedampdict=[]
    correctedampstartdict=dict()
    correctedampenddict=dict()
    counter = 0
    procampres=manager.dict()

    for amplicon in amplicons:
        #print amplicon
    	counter+=1
    	ampstart = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
    	ampstop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
    	ampstartdict[counter]=ampstart
    	ampenddict[counter]=ampstop
    	ampdict.append((counter,ampstart,ampstop))
        procampres[counter]=0

    procampres["DO"]=0
    procampres["HF"]=0
    procampres["NH"]=0
    procampres["BF"]=0


    print "******AMP DICTIONARY*******"
    print type(ampstartdict)
    print ampstartdict
    readprediction=dict()

    print procampres

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
            #catch negative values
            if start < 1:
                start = 1
            stop = position + readuntilrange
            ranges[sequence].append((start,stop))

    for sequence in ranges:
        ranges[sequence] = list(merge_ranges((ranges[sequence])))






    correctedamplicons=list()
    for amplicon in amplicons:
        ampstart = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
    	ampstop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
        correctedamplicons.append(amplicon.split(':',1)[0]+':'+str(correctposition(ampstart,ranges,amplicon.split(':',1)[0]))+'-'+str(correctposition(ampstop,ranges,amplicon.split(':',1)[0])))

    counter = 0
    for amplicon in correctedamplicons:
    	counter+=1
    	ampstart = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
    	ampstop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
    	correctedampstartdict[counter]=ampstart
    	correctedampenddict[counter]=ampstop
    	correctedampdict.append((counter,ampstart,ampstop))

    kmerhashTRU = get_custom_fasta(fasta_file,model_kmer_means_template,ranges)


    print "Attempting to match reads and split into folders based on 250 events, excluding the first 50."

    d=list()
    filenamecounter=0
    for dirpath, dirnames, files in os.walk(os.path.join(args.watchdir)):
        for f in files:
            if (f.endswith('.fast5')):
                filename = os.path.join(dirpath,f)
                filenamecounter+=1
                d.append([filename,procampres,kmerhashTRU,correctedamplicons,correctedampstartdict,correctedampenddict])
    procdata=tuple(d)

    procampres["TF"]=filenamecounter



    results = p.map(process_hdf5, (procdata),chunksize=1)
    p.close()
    masterreadprediction=dict()
    for element in results:
        for amplicon in element:
            if (amplicon not in masterreadprediction):
                masterreadprediction[amplicon]=dict()
            for quality in element[amplicon]:
                if (quality not in masterreadprediction[amplicon]):
                    masterreadprediction[amplicon][quality]=dict()
                for filename in element[amplicon][quality]:
                    if (filename not in masterreadprediction[amplicon][quality]):
                        masterreadprediction[amplicon][quality][filename]=dict()
                    masterreadprediction[amplicon][quality][filename]["name"]=element[amplicon][quality][filename]["name"]
                    masterreadprediction[amplicon][quality][filename]["matchdistance"]=element[amplicon][quality][filename]["matchdistance"]

    print "Amplicon Read Counts"
    for amplicon in masterreadprediction:
        numberofreads = 0
        for i in range(5):
            try:
                if len(masterreadprediction[amplicon][i].keys()) > 0:
                    numberofreads += len(masterreadprediction[amplicon][i].keys())
            except Exception, err:
                print "",
        print "Amplicon Number:",amplicon,"Reads:",numberofreads

    #print masterreadprediction

    print "Copying Amplicon Data"
    for amplicon in masterreadprediction:
        print "Amplicon Number",amplicon
    	counter = 0
        for i in range(5):
            try:
                if (len(masterreadprediction[amplicon][i].keys())>0):
                    if (args.verbose is True):
                        print len(masterreadprediction[amplicon][i].keys())
                    if (counter < args.depth):
                        ordered0 = OrderedDict(sorted(masterreadprediction[amplicon][i].iteritems(), key=lambda x: x[1]['matchdistance']))
                        for read in ordered0:
                            if (args.verbose is True):
                                print "Checking if read is basecalled"
                                print read
                            hdf = h5py.File(read, 'r')
                            readstatus=False
                            if check_basecalled(hdf) is True:
                                readstatus=True
                            hdf.close()
                            if (args.verbose is True):
                                print read, ordered0[read]["matchdistance"]
                            if not os.path.exists(args.targetpath):
                                os.makedirs(args.targetpath)
                            destdir = os.path.join(args.targetpath)
                            if not os.path.exists(destdir):
                                os.makedirs(destdir)
                            destdir2 = os.path.join(destdir,str(amplicon))
                            if not os.path.exists(destdir2):
                                os.makedirs(destdir2)
                            if readstatus is True:
                                destdir2 = os.path.join(destdir,str(amplicon),"downloads")
                            else:
                                destdir2 = os.path.join(destdir,str(amplicon))
                            if not os.path.exists(destdir2):
                                os.makedirs(destdir2)
                            try:
                                filetocheck = os.path.split(read)
                                sourcefile = read
                                destfile = os.path.join(destdir2,filetocheck[1])
                                if (args.verbose is True):
                                    print "sourcefile is:",sourcefile
                                    print "destfile is:",destfile
                                try:
                                    shutil.copy(sourcefile,destfile)
                                except Exception, err:
                                    print "File Copy Failed",err
                            except Exception, err:
                                print "Unknown File Error"
                                print "This might help:", err
                            counter += 1
                            if counter >= args.depth:
                                break
            except Exception, err:
                if (args.verbose is True):
                    print "No reads of class "+str(i)

    sys.exit()
