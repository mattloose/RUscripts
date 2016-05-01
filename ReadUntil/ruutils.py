import csv
import numpy as np
import sys,os,re
import urllib2
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from StringIO import StringIO
import sklearn.preprocessing
import multiprocessing
import ctypes
import mlpy
import time
import h5py
import psutil




def check_files(listoffiles):
    """
    Supplied a list of files, this routine checks to see if they exist. If the files cannot be found, the program exits.
    """
    print "\n\rChecking Files..."
    for file in listoffiles:
        if not os.path.isfile(file):
            print "\n\r**! One of the files you supplied cannot be found. Please check:\n\r\n\r"+str(file)+"\n\r\n\r**! This program will now exit.\n\r"
            sys.exit()
    print "\n\rAll OK.\n\r"

def checkfasta(fasta):
    try:
        for record in SeqIO.parse(fasta, 'fasta'):
            print record.id
    except:
        print "Fasta does not appear to be valid.\n\r"
        print "This program will now exit. Please try again.\n\r"
        sys.exit()

def correctposition(value,ranges,sequence):
    correction = 0
    for range in ranges[sequence]:
        if value >= range[0] and value <= range[1]:
            return value - range[0] + correction + 1
        correction = correction + (range[1]-range[0] + 1)


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


def scaleLocally(a, sz):
	normWinSz=sz
	n = (normWinSz/2)+1 # eg win 64 -> n == 33
	start = scale(a[ : normWinSz+1 ])[:n]
	end   = scale(a[ -(normWinSz+1) : ])[-n:]
	mid = [ scale(a[i-n:i+n])[n] for i in range(n,len(a)-n) ]
	a = np.hstack([start, mid, end])
	#print len(a)
	#exit()

	#a = np.hstack(map(scale, split(a, sz)))
	return a

def scale(a):
	'''
	return my_scale(a)
	'''
	return sklearn.preprocessing.scale(a, axis=0, with_mean=True, with_std=True, copy=True)

###########################################################
def file_dict_of_folder(path):

	file_list_dict=dict()
	ref_list_dict=dict()

	print "File Dict Of Folder Called"
	if os.path.isdir(path):
		print "caching existing fast5 files in: %s" % (path)
		for path, dirs, files in os.walk(path) :
			for f in files:
				#if (("downloads" in path )):
				if ("muxscan" not in f and f.endswith(".fast5") ):
					file_list_dict[os.path.join(path, f)]=os.stat(os.path.join(path, f)).st_mtime

	print "found %d existing fast5 files to process first." % (len(file_list_dict) )
    	return file_list_dict




def genome_worker((channel_id, data,readstarttime,args,seqlen,seqids,threedarray)):
    #print "mpworker"
    #current = multiprocessing.current_process()
    #print "Process is",current.name,"running on",channel_id
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
        #if args.verbose is True:
        print "We have a timeout ",channel_id
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
            #return "skip",channel_id,readnumber,data.events[0].start,"squiggleres"
            #squiggleres = squiggle_search2(squiggle,channel_id,readnumber,kmerhash,seqlen,args)
            squiggleres = squiggle_search2(squiggle,channel_id,readnumber,args,seqids,threedarray,seqlen)
            if args.verbose is True: print squiggleres
            #result = go_or_no(squiggleres[0],squiggleres[2],squiggleres[3],seqlen,args)
            seqid,_,direction,position,__,___,____=squiggleres
            #seqid,_,direction,position,__=squiggleres
            result = go_or_no(seqid,direction,position,seqlen,args)
            if args.verbose is True: print "Result ",result
            #print result,channel_id,readnumber,data.events[0].start,squiggleres
            return result,channel_id,readnumber,data.events[0].start,squiggleres
        except Exception, err:
            err_string="Time Warping Stuff : %s" % ( err)
            print >>sys.stderr, err_string

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
                if position >= ( start - balance ) and position <= stop:
                #if position >= ( length - stop - balance) and position <= ( length - start ):
                    return "Sequence"
    return "Skip"



######################################################################
def extractsquig(events):
    squiggle=list()
    for event in events:
        squiggle.append(np.float32(event.mean))
    return(squiggle)



def squiggle_search2(squiggle,channel_id,read_id,args,seqids,threedarray,seqlen):
    result=[]
    blocksize=200000
    overlap=blocksize-500
    for ref in seqids:
        refid=seqids.index(ref)
        Rprime,Fprime=threedarray[refid]
        #queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
        queryarray = sklearn.preprocessing.scale(np.array(squiggle,dtype=float),axis=0,with_mean=True,with_std=True,copy=True)
        refsubset = Fprime
        indexes = np.array(xrange(len(refsubset)))
        subrefs = [refsubset[i:i+blocksize]for i in indexes[::overlap]]
        for blockid,ref_ in enumerate(subrefs):
            #current = multiprocessing.current_process()
            tic = time.time()
            dist, cost, path = mlpy.dtw_subsequence(queryarray,ref_)
            #result.append((dist,ref,"F",path[1][0],path[1][-1],path[0][0],path[0][-1]))
            result.append((dist,ref,"F",path[1][0]+(blockid*overlap),path[1][-1]+(blockid*overlap),path[0][0],path[0][-1]))
            #print "Blockid", blockid, time.time()-tic
        refsubset = Rprime
        subrefs = [refsubset[i:i+blocksize]for i in indexes[::overlap]]
        for blockid,ref_ in enumerate(subrefs):
            #print "Blockid", blockid, time.time()
            dist, cost, path = mlpy.dtw_subsequence(queryarray,ref_)
            #result.append((dist,ref,"R",path[1][0]+(blockid*overlap),ref))
            result.append((dist,ref,"R",(len(Rprime)-(path[1][-1]+(blockid*overlap))),(len(Rprime)-(path[1][0]+(blockid*overlap))),path[0][0],path[0][-1]))
    # Note first two elements flipped for return deliberately.
    distanceR,seqmatchnameR,frR,rsR,reR,qsR,qeR=sorted(result,key=lambda result: result[0])[0]
    return seqmatchnameR,distanceR,frR,rsR,reR,qsR,qeR




######################################################
def get_seq_len(ref_fasta):
    seqlens=dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        seq=record.seq
        seqlens[record.id]=len(seq)
    return seqlens

#######################################################################

'''
This version includes some alternative scalings - currently deprecated
######################################################
def process_ref_fasta_raw(ref_fasta,model_kmer_means,args,kmer_len):
    if (args.verbose is True):
        print "processing the reference fasta."
    kmer_means=dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id]=dict()
        kmer_means[record.id]["F"]=list()
        kmer_means[record.id]["R"]=list()
        kmer_means[record.id]["Fprime"]=list()
        kmer_means[record.id]["Rprime"]=list()
        kmer_means[record.id]["Floc"]=list()
        kmer_means[record.id]["Rloc"]=list()
        if (args.verbose is True):
            print "ID", record.id
            print "length", len(record.seq)
            print "FORWARD STRAND"
        seq = record.seq
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))
        if (args.verbose is True):
            print "REVERSE STRAND"
        seq = revcomp = record.seq.reverse_complement()
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            kmer_means[record.id]["R"].append(float(model_kmer_means[kmer]))
        kmer_means[record.id]["Fprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["F"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]["Floc"]=scaleLocally(kmer_means[record.id]["F"],250)
        kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]["Rloc"]=scaleLocally(kmer_means[record.id]["R"],250)
    return kmer_means
'''
def get_custom_fasta(ref_fasta,subsectionlist,args,model_kmer_means,kmer_len):
    if (args.verbose is True):
        print "Generating a custom fasta"
    sequencedict=dict()
    for sequence in subsectionlist:
        if (args.verbose is True):
            print sequence
        for record in SeqIO.parse(ref_fasta, 'fasta'):
            if (record.id == sequence):
                if (sequence not in sequencedict):
                    sequencedict[sequence]=list()
                for sections in subsectionlist[sequence]:
                    start = sections[0]
                    end = sections[1]
                    if (len(sequencedict[sequence])>0):
                        sequencedict[sequence]=str(sequencedict[sequence])+str(record.seq[sections[0]-1:sections[1]-1])
                    else:
                        sequencedict[sequence]=str(record.seq[sections[0]-1:sections[1]-1])
    if (args.verbose is True):
        print "processing the custom fasta"
    kmer_means=dict()
    for sequence in sequencedict:
        kmer_means[record.id]=dict()
        tmp=dict()
        tmp2=dict()
        tmp["F"]=list()
        tmp["R"]=list()
        tmp["Fprime"]=list()
        tmp["Rprime"]=list()
        print "ID", record.id
        print "length", len(record.seq)
        print "FORWARD STRAND"
        seq = Seq(sequencedict[sequence], generic_dna)
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            tmp["F"].append(float(model_kmer_means[kmer]))
        print "REVERSE STRAND"
        seq = revcomp = seq.reverse_complement()
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            tmp["R"].append(float(model_kmer_means[kmer]))
        tmp2["Fprime"]=sklearn.preprocessing.scale(tmp["F"], axis=0, with_mean=True, with_std=True, copy=True)
        tmp2["Rprime"]=sklearn.preprocessing.scale(tmp["R"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]=tmp2
    '''From this dictionary we will return a pair consisting of a list of keys(lookup for sequence name) and a
    3D array each slice of which relates to the seqid,forward and reverse and then the values. This will then
    be used as a numpy shared memory multiprocessing array. We hope.
    Caution - the dictionary returns in the wrong order.
    '''

    items=kmer_means.items()
    '''for k,v in kmer_means.items():
        for x,y in kmer_means[k].items():
            print "idiot check",k,x
            '''
    items_=map(processItems,items)
    seqids,arrays=zip(*items_)
    z=len(seqids)
    print arrays
    r,c=list(arrays)[0].shape
    threedarray=multiprocessing.Array(ctypes.c_double,z*r*c)
    threedarrayshared_array = np.ctypeslib.as_array(threedarray.get_obj())
    a = np.array(arrays,dtype=np.float32)
    threedarrayshared_array = a
    return seqids,threedarrayshared_array


def process_ref_fasta(ref_fasta,model_kmer_means,kmer_len):
    print "processing the reference fasta."
    kmer_means=dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id]=dict()
        tmp=dict()
        tmp2=dict()
        tmp["F"]=list()
        tmp["R"]=list()
        tmp["Fprime"]=list()
        tmp["Rprime"]=list()
        print "ID", record.id
        print "length", len(record.seq)
        print "FORWARD STRAND"
        seq = record.seq
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            tmp["F"].append(float(model_kmer_means[kmer]))
        print "REVERSE STRAND"
        seq = revcomp = record.seq.reverse_complement()
        for x in range(len(seq)+1-kmer_len):
            kmer = str(seq[x:x+kmer_len])
            tmp["R"].append(float(model_kmer_means[kmer]))
        tmp2["Fprime"]=sklearn.preprocessing.scale(tmp["F"], axis=0, with_mean=True, with_std=True, copy=True)
        tmp2["Rprime"]=sklearn.preprocessing.scale(tmp["R"], axis=0, with_mean=True, with_std=True, copy=True)
        kmer_means[record.id]=tmp2
    '''From this dictionary we will return a pair consisting of a list of keys(lookup for sequence name) and a
    3D array each slice of which relates to the seqid,forward and reverse and then the values. This will then
    be used as a numpy shared memory multiprocessing array. We hope.
    Caution - the dictionary returns in the wrong order.
    '''

    items=kmer_means.items()
    '''for k,v in kmer_means.items():
        for x,y in kmer_means[k].items():
            print "idiot check",k,x
            '''
    items_=map(processItems,items)
    seqids,arrays=zip(*items_)
    z=len(seqids)
    r,c=list(arrays)[0].shape
    threedarray=multiprocessing.Array(ctypes.c_double,z*r*c)
    threedarrayshared_array = np.ctypeslib.as_array(threedarray.get_obj())
    a = np.array(arrays,dtype=np.float32)
    threedarrayshared_array = a
    return seqids,threedarrayshared_array


def processItems((seqid,d)):
    result=[]
    for _,l in d.items():
        result.append(l)
    return seqid,np.array(result)


#######################################################################

######################################################
def process_model_file(model_file):
    '''
    This function parses a csv model file into a dict to provide model simulation and return the dict and kmer lenght
    '''
    model_kmers = dict()
    with open(model_file, 'rb') as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        d = list(reader)
        header = d[0]
        #print header
        meancolumnid = header.index("level_mean")
        kmercolumnid = header.index("kmer")

        d2 = np.array(d[1:])
        #for x in d2:
        #    print x
        model_kmers=dict(d2[:,[kmercolumnid,meancolumnid]])
        kmer_length=len(d2[0][kmercolumnid])
    return model_kmers,kmer_length



######################################################
def _urlopen(url, *args):

    """Open a URL, without using a proxy for localhost.

    While the no_proxy environment variable or the Windows "Bypass proxy

    server for local addresses" option should be set in a normal proxy

    configuration, the latter does not affect requests by IP address. This

    is apparently "by design" (http://support.microsoft.com/kb/262981).

    This method wraps urllib2.urlopen and disables any set proxy for

    localhost addresses.

    """
    #print "_urlopen called"
    try:

        host = url.get_host().split(':')[0]
    #    print "Yo mamma",host

    except AttributeError:

        host = urlparse.urlparse(url).netloc.split(':')[0]

    import socket

    #print "socket in"

    # NB: gethostbyname only supports IPv4

    # this works even if host is already an IP address

    addr = socket.gethostbyname(host)

    #if addr.startswith('127.'):
    #    print "returning no proxt"
    #    return _no_proxy_opener.open(url, *args)

    #else:
    #    print "returning urlib"
    return urllib2.urlopen(url, *args)


def execute_command_as_string( data, host = None, port = None):
        '''
        This function is used for sending messages to and from minKNOW.
        '''
        host_name = host

        port_number = port
        #print data

        url = 'http://%s:%s%s' % (host_name, port_number, "/jsonrpc")
        #print url

        req = urllib2.Request(url, data=data, headers={'Content-Length': str(len(data)), 'Content-Type': 'application/json'})
        f = None
        try:
            f = _urlopen(req)
        except Exception, err:
            err_string = "Failed to connect to minKNOW. Likely reasons include minKNOW not running, the wrong IP address for the minKNOW server or firewall issues. Two way control will not be possible with minKNOW."
            print err_string, err
        json_respond = json.loads(f.read())

        f.close()

        return json_respond

def send_message(message,ipadd):
    '''
    This function is used to send messages to and from minKNOW
    '''
    message_to_send = '{"id":"1", "method":"user_message","params":{"content":"%s"}}' %(message)
    results = execute_command_as_string(message_to_send,ipadd,8000)
    return results


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n\r")


def die_nicely(oper):
    print "terminating sub-processes...."

    pid = os.getpid()

    # Tell minup to terminate
    if oper == "windows":
        sys.exit()
        # -- sending minup pid a Ctrl-C signal
        # -- this also cleanly closes subprocesses and threads ....
        #ctypes.windll.kernel32.GenerateConsoleCtrlEvent(0, pid) # 0 => Ctrl-C
    else:
        process = psutil.Process(pid)
        for proc in process.children(recursive=True):
            proc.kill()
        process.kill()
