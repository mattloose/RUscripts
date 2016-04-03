# RUscripts

To get started either download this repository using the ZIP options or (preferably)

git clone https://github.com/mattloose/RUscripts



Important Note:

Running read until will influence the behaviour of your flow cell and change the output of your sequencing experiment. You are strongly advised to run simulations of read until prior to running on a live flow cell. The code as presented here is a demonstration of read until and one implementation. Users run this code entirely at their own risk.

#Getting help with read until.

We are happy to help where we can. Please feel free to contact via twitter ({at}mattloose) or email matt.loose{@}nottingham.ac.uk . We will try to update this read me file whenever significant changes are made to the API or minKNOW versions that might affect read until. For issues with the API itself or the ws_event_sampler program, users will most likely be redirected to ONT.

#Getting started

Scripts for implementing read until and other examples.

To use read until at this time you have to use the command line. To access the command line on windows, click on the start menu and type 'cmd' in the search box. You should see a program called 'cmd'.

Once at the command line, navigate to the folder containing this Readme file.

All scripts are python 2.7.

We assume that the default version of python in your path is the ONT Anaconda installation. To test this, type

    C:\path\to\RUscripts>python

And you should see:

    Python 2.7.5 |Anaconda 1.8.0 (64-bit)| (default, Jul  1 2013, 12:37:52) [MSC v.1500 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

Note that this shows Anaconda 1.8.0

To exit this window type

    >>> exit()

Please note these scripts depend upon APIs which are not provided here. The Read Until API is available from Oxford Nanopore on request to MAP participants.

# minKNOW version compatability statement

These scripts and files are compatible with versions of minKNOW predating the v 0.51.1.62 (released February 22nd 2016). Read Until in this latest version of minKNOW is currently unstable. ONT are fixing this at this time. We provide scripts that compensate for these errors in a simulation, but these will not work for a live read until run at this time.

Currently supported versions of minKNOW are v 0.48.1.3 to 0.51.1.51.

We highlight where the version of minKNOW used might impact results in the documentation below.

*** You must have minKNOW installed on the windows machine you are running this code on. The version of Anaconda it comes with contains some of the required packages. These scripts will also run on Linux or OSX but at this time we provide no explicit documentation to support this.

# Python modules required for running read until using source code.

h5py

configargparse

biopython

mlpy # Must be version 3.5.0 or above - see below

watchdog

sklearn #Available as default with Anaconda


# To install these modules on the ONT Anaconda type:

    C:\path\to\RUscripts>easy_install thrift ws4py h5py configargparse biopython watchdog

# To install mlpy

Installing mlpy with easy_install results in a version less than that required. To overcome this on windows systems visit http://www.lfd.uci.edu/~gohlke/pythonlibs/#mlpy to obtain an mlpy wheel compiled for the appropriate version of windows and python.

On our system this file is called:

mlpy-3.5.0-cp27-none-win_amd64.whl


To install this you need to first install pip.

    C:\path\to\RUscripts>easy_install install pip

Then to install your downloaded wheel:

    C:\path\to\RUscripts>pip install \path\to\mlpy-3.5.0-cp27-non-win_amd64.whl

Now you should have all the prerequisites required to run these scripts on windows using the native python.


# General comments

All scripts here adhere to the following standard arguments:

-h Will print a help statement and exit.
-ver Will print version information and exit.
-v Will print additional debug information and more verbose output to screen while the code is executing.

Additional options are described below for each script.

# Files and Folders in this repository

ampbalance.py - an illustrative script for selecting a fixed number of reads from individual amplicons using squiggle matching alone.

ampbalancetest (folder - contains 110 example fast5 files)

ampliconSPLIT.py - an illustrative script that will split reads into folders based on squiggle matching to individual amplicons.

exampleread\llssbzms2p35x_lambda11ladderup_1208_1_ch11_file27_strand.fast5

getmodels.py - A script for extracting model files to generate a squiggle reference for a given sequence.

J02459.fasta - A reference fasta file for bacteriophage lambda. This is the reference used for all our work.

lambda_amplicons.txt - A file that defines the amplicons used in this study. The file is in the format seq_id:start-stop - so here the seq_id is J02456 (as in the file J02459.fasta) and the start and stop coordinates match our amplicons eg. J02456:52-1980

LICENSE - MIT License

README.md - this readme file.

ReadUntil\aReadUntil.py - A script for amplicon balancing via the read until API.

ReadUntil\gReadUntil.py - A script for selecting a specific region of the genome via the read until API

ReadUntil\lambda_amplicons.fasta - This is a reference file containing copies of individual amplicons for bacteriophage lambda. It is used in testing the behaviour of aReadUntil.py - described below.

ReadUntil\ruutils.py - This file contains a number of useful functions for various scripts. It will never be run itself but is called by many of the other scripts.

test_gReadUntil.py - A script for testing the functionality of gReadUntil.py in the absence of the API itself. 

RUtestset (folder - contains 55 example fast5 files)


# Running Offline Read Until Scripts

# getmodels.py

In order to squiggle match between an incoming read and a reference sequence it is necessary to generate a reference in squiggle space. To do this you require a model file to enable the conversion of kmers in base space to a current value. Currently the appropriate model files are embedded within reads as they are returned by the basecaller.

It is important to have a model file for the chemistry and pore type that you are using.

For read until you only require a template model file. If you wish to run ampbalance (which identifies potential 2D reads from a dataset you also require the complement model file).

It is not necessary to generate a model file each time you carry out a new read until experiment. Instead, new files should be obtained whenever a chemistry or pore change occurs. We suggest that it is good practice to extract a model file from the most recent read data available to you when preparing for a read until run as these models have been silently updated in the past.

Ideally the read from which to extract a model should come from the pass folder of a base called run. To extract both template and complement models, this read must be from a 2D base called run. getmodels.py will output whichever models it finds within the file.

To print the getmodels.py help statement at the prompt type:

    C:\path\to\RUscripts>python getmodels.py -h

which will output:

    usage: getmodels.py [-h] [-r READ] [-v] [-ver]

    getmodels: A program to extract model files from Nanopore reads. Ideally you
    should provide a read that has given a 2D read - i.e comes from the pass
    folder on a 2D basecalling run. For Read Until purposes, extraction of just a
    template model is sufficient. Model files will be appended with a number to
    identify the kmer length of the model.

    optional arguments:
    -h, --help            show this help message and exit
    -r READ, --read READ  Provide a read file to extract the current model from.
    -v, --verbose-true    Print detailed messages while processing files to aid
                        in debugging.
    -ver, --version       Print version information and exit.

To run getmodels.py at the prompt ($) type:

    $python getmodels.py -r /path/to/your/read.fast5

We provide a suitable read in the folder exampleread

So:

    C:\path\to\RUscripts>python getmodels.py -r exampleread/llssbzms2p35x_lambda11ladderup_1208_1_ch11_file27_strand.fast5

will output:

    Processing exampleread/llssbzms2p35x_lambda11ladderup_1208_1_ch11_file27_strand.fast5
    Looking for template model.
    ***: Model Found
    template model file write completed.
    Kmer length is: 6
    File format is:
    Kmer	Mean	Standard Dev
    Looking for complement model.
    ***: Model Found
    complement model file write completed.
    Kmer length is: 6
    File format is:
    Kmer	Mean	Standard Dev

This will have generated two files - one for the template and one for the complement model.

They will be named:

[readtype]_[model_name]_[kmerlength].model

Where [readtype] is either template or complement, [model_name] is the pore/chemistry type (in this case r7.3_e6_70bps_6mer) and [kmerlen] is the locally calculated length of the kmers used in the model.

These model files are used by all other scripts here. Note that the model files will be overwritten by getmodels.py if it is run again on another read with the same pore/chemistry type and kmer length. This is not a problem.

*** We recommend that you update the model regularly.

# ampbalance.py

Ampbalance is a script which will map entire 2D reads to a reference sequence, identify those likely to give a good 2D read by determining if the template and complement reads overlap with one another appropriately and subsequently copy a defined number of those reads up to a specific depth into a folder for further analysis. This script was originally written to facilitate the sequencing and analysis of Ebola virus and to simplify the transfer of read data prior to base calling where network bandwidth was limited (see Quick, J., Loman, N. J., Duraffour, S., Simpson, J. T. & Severi, E. Real-time, portable genome sequencing for Ebola Surveillance. Nature, doi:10.1038/nature16996 (2016).)

Unlike read until applications this script uses the entire read.

Here the script is implemented on reads that have been basecalled to facilitate testing and comparisons.

To print the ampbalance.py help statement at the prompt ($) type:

    C:\path\to\RUscripts>python ampbalance.py -h

which will output:

    usage: ampbalance.py [-h] -fasta FASTA -ids IDS -w WATCHDIR -o TARGETPATH -d
                         DEPTH -procs PROCS [-cautious] [-l LENGTH] -t TEMP_MODEL
                         -c COMP_MODEL [-v] [-ver]

    ampbalance: A program designed to balance amplicons from a specific reference
    sequence post sequencing on ONT minIONs but prebasecalling. Developed by Matt
    Loose @mattloose or matt.loose@nottingham.ac.uk for help! Args that start with
    '--' (eg. --reference_fasta_file) can also be set in a config file
    (/path/to/your/script/amp.config or ) by using .ini or
    .yaml-style syntax (eg. reference_fasta_file=value). If an arg is specified in
    more than one place, then command-line values override config file values
    which override defaults.

    optional arguments:
      -h, --help            show this help message and exit
      -fasta FASTA, --reference_fasta_file FASTA
                            The fasta format file for the reference sequence for
                            your organism.
      -ids IDS, --reference_amplicon_positions IDS
                            A file containing a list of amplicon positions defined
                            for the reference sequence. 1 amplicon per line in the
                            format fasta_sequence_name:start-stop e.g
                            J02459:27-1938
      -w WATCHDIR, --watch-dir WATCHDIR
                            The path to the folder containing the downloads
                            directory with fast5 reads to analyse - e.g.
                            C:\data\minion\downloads (for windows).
      -o TARGETPATH, --output-dir TARGETPATH
                            The path to the destination folder for the
                            preprocessed reads
      -d DEPTH, --depth DEPTH
                            The desired coverage depth for each amplicon. Note
                            this is unlikely to be achieved for each amplicon and
                            should probably be an overestimate of the minimum
                            coverage required.
      -procs PROCS, --proc_num PROCS
                            The number of processors to run this on.
      -cautious, --cautious
                            DTW of long reads on low memory systems can cause
                            unexpected crashes. This option will prevent automatic
                            skipping on any reads over 10,000 events. You can
                            optionally increase this length with the -l parameter.
                            USE WITH CAUTION AS THIS MAY CAUSE A SYSTEM TO CRASH.
      -l LENGTH, --length LENGTH
                            A limit on the length of read that ampbalance will
                            attempt to align using DTW - Long reads can cause
                            problems on low memory systems
      -t TEMP_MODEL, --template_model TEMP_MODEL
                            The appropriate template model file to use
      -c COMP_MODEL, --complement_model COMP_MODEL
                            The appropriate complement model file to use
      -v, --verbose-true    Print detailed messages while processing files.
      -ver, --version       show program's version number and exit

This script is designed to match amplicon sequences of known and approximately uniform length to a reference sequence. The reference sequence should be a single sequence.

We provide an example set of reads using lambda, found in the ampbalancetest folder. A typical usage command would be:

    C:\path\to\RUscripts>python ampbalance.py -fasta J02459.fasta -ids lambda_amplicons.txt -w ampbalancetest -o ampbalanceoutput -d 5 -procs 4 -t template_r7.3_e6_70bps_6mer_6.model -c complement_r7.3_e6_70bps_6mer_6.model -l 3000

This will output the following:

    Reading amplicons
    ******AMP DICTIONARY*******
    <type 'dict'>
    {1: 52, 2: 2065, 3: 4070, 4: 6059, 5: 8012, 6: 10008, 7: 12006, 8: 14011, 9: 16076, 10: 18022, 11: 20053}
    {'DO': 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 'BF': 0, 'HF': 0, 'NH': 0}
    Now we are going to try and open the raw reads and do the same as we have done above...
    {'DO': 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 'BF': 0, 'TF': 109, 'HF': 4, 'NH': 0}
    {'DO': 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0, 6: 0, 7: 1, 8: 0, 9: 0, 10: 0, 11: 0, 'BF': 0, 'TF': 108, 'HF': 5, 'NH': 0}
    .
    . (many more lines here)
    .
    {'DO': 0, 1: 10, 2: 10, 3: 10, 4: 10, 5: 9, 6: 10, 7: 10, 8: 10, 9: 10, 10: 10, 11: 10, 'BF': 0, 'TF': 1, 'HF': 110, 'NH': 0}
    {'DO': 0, 1: 10, 2: 10, 3: 10, 4: 10, 5: 10, 6: 10, 7: 10, 8: 10, 9: 10, 10: 10, 11: 10, 'BF': 0, 'TF': 0, 'HF': 110, 'NH': 0}
    Amplicon Read Counts
    Amplicon Number: 1 Reads: 10
    Amplicon Number: 2 Reads: 10
    Amplicon Number: 3 Reads: 10
    Amplicon Number: 4 Reads: 10
    Amplicon Number: 5 Reads: 10
    Amplicon Number: 6 Reads: 10
    Amplicon Number: 7 Reads: 10
    Amplicon Number: 8 Reads: 10
    Amplicon Number: 9 Reads: 10
    Amplicon Number: 10 Reads: 10
    Amplicon Number: 11 Reads: 10
    Copying Amplicon Data
    Amplicon Number 1
    Amplicon Number 2
    Amplicon Number 3
    Amplicon Number 4
    Amplicon Number 5
    Amplicon Number 6
    Amplicon Number 7
    Amplicon Number 8
    Amplicon Number 9
    Amplicon Number 10
    Amplicon Number 11

Key:
DO: Template and Complement Strands don't overlap on amplicons as expected. Read Discarded. Unlikely to generate good 2D.
1...11: Amplicon numbers followed by current coverage counts satisfying highest stringency tests.
BF: Bad Files - these files contain more events than the -length threshold. They are unlikely to be of interest in an amplicon run and are skipped unless the -cautious flag is set.
TF: Total Files - this counts down the number of files left to process.
HF: Hairpin Found - this is the number of reads containing a hairpin - these are tested.
NH: No Hairpin - these reads are ignored as no hairpin has been found.

This program will output files to the specified output directory. Files which are already basecalled will be written to a subfolder called "Downloads". Raw files will just be written to the specified output directory.

Note that the example read set here has been preselected to contain 10 reads mapping to each amplicon and all the reads are 2D - thus the HF value is 110 and the DO value is 0.

It is also important to note that this script uses the hairpin_found flag in the read files. Since minKNOW version 0.51.1.62 (released February 22nd 2016) the writing of this flag to read files has been problematic. Thus reads generated with this version of minKNOW may appear to have far fewer 2D reads than they really do.


# ampliconSPLIT.py

ampliconSPLIT.py simulates read until on either raw or basecalled reads. It ignores the first 50 events of a read and matches the subsequent 250 events against a reference squiggle. Read until only ever processes template data and so only a template model is required.

This script will process reads into subfolders in the targetpath corresponding to each amplicon described in the -ids file - e.g /targetpath/{amplicon_number}. If the read is basecalled, it will be placed in a downloads subfolder - e.g /targetpath/{amplicon_number}/downloads .

To print the ampliconSPLIT.py help statement at the prompt ($) type:

    C:\path\to\RUscripts>python ampliconSPLIT.py -h

which will output:

    usage: ampliconSPLIT.py [-h] -fasta FASTA -ids IDS -w WATCHDIR -o TARGETPATH
                        -d DEPTH -procs PROCS -t TEMP_MODEL [-v] [-ver]

    ampliconSPLIT: A program designed to identify and group individual amplicons
    from minION reads prior to base calling. The depth setting limits the number
    of reads copied to each sub folder. Developed by Matt Loose @mattloose or
    matt.loose@nottingham.ac.uk for help! Args that start with '--' (eg.
    --reference_fasta_file) can also be set in a config file
    (/Users/mattloose/fixes/RUscripts/UPDATES/amp.config or ) by using .ini or
    .yaml-style syntax (eg. reference_fasta_file=value). If an arg is specified in
    more than one place, then command-line values override config file values
    which override defaults.

    optional arguments:
    -h, --help            show this help message and exit
    -fasta FASTA, --reference_fasta_file FASTA
                        The fasta format file for the reference sequence for
                        your organism.
    -ids IDS, --reference_amplicon_positions IDS
                        A file containing a list of amplicon positions defined
                        for the reference sequence. 1 amplicon per line in the
                        format fasta_sequence_name:start-stop e.g
                        EM_079517:27-1938
    -w WATCHDIR, --watch-dir WATCHDIR
                        The path to the folder containing the downloads
                        directory with fast5 reads to analyse - e.g.
                        C:\data\minion\downloads (for windows).
    -o TARGETPATH, --output-dir TARGETPATH
                        The path to the destination folder for the
                        preprocessed reads
    -d DEPTH, --depth DEPTH
                        The desired coverage depth for each amplicon. Note
                        this is unlikely to be achieved for each amplicon and
                        should probably be an overestimate of the minimum
                        coverage required.
    -procs PROCS, --proc_num PROCS
                        The number of processors to run this on.
    -t TEMP_MODEL, --template_model TEMP_MODEL
                        The appropriate template model file to use. This file
                        can be generated uing the getmodels.py script.
    -v, --verbose-true    Print detailed messages while processing files.
    -ver, --version       show program's version number and exit

An example run command using test data:

    C:\path\to\RUscripts>python ampliconSPLIT.py -fasta J02459.fasta -ids lambda_amplicons.txt -w RUtestset/ -o test -d 10 -procs 8 -t template_r7.3_e6_70bps_6mer_6.model

Which will output:

    processing the reference fasta.
    ID J02459
    length 48502
    FORWARD STRAND
    REVERSE STRAND
    Groking amplicons
    ******AMP DICTIONARY*******
    <type 'dict'>
    {1: 52, 2: 2065, 3: 4070, 4: 6059, 5: 8012, 6: 10008, 7: 12006, 8: 14011, 9: 16076, 10: 18022, 11: 20053}
    {'DO': 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 'BF': 0, 'HF': 0, 'NH': 0}
    We want to build a custom reference that is smaller than the original reference.
    First get a list of all the positions we will need to search
    Generating a custom fasta
    J02459
    processing the custom fasta
    Attempting to match reads and split into folders based on 250 events, excluding the first 50.
    Amplicon Read Counts
        Amplicon Number: 1 Reads: 4
        Amplicon Number: 2 Reads: 5
        Amplicon Number: 3 Reads: 4
        Amplicon Number: 4 Reads: 5
        Amplicon Number: 5 Reads: 6
        Amplicon Number: 6 Reads: 6
        Amplicon Number: 7 Reads: 5
        Amplicon Number: 8 Reads: 5
        Amplicon Number: 9 Reads: 5
        Amplicon Number: 10 Reads: 5
        Amplicon Number: 11 Reads: 5
    Copying Amplicon Data
    Amplicon Number 1
    Amplicon Number 2
    Amplicon Number 3
    Amplicon Number 4
    Amplicon Number 5
    Amplicon Number 6
    Amplicon Number 7
    Amplicon Number 8
    Amplicon Number 9
    Amplicon Number 10
    Amplicon Number 11

The depth parameter (-d) sets the number of reads that will be copied. Setting this to a value greater than the number of reads analysed will sort all reads in the dataset.


# test_gReadUntil.py

*** Note this script is located with the ReadUntil folder of this repository

This script runs entirely independently of the read until API and allows for simulation of selective sequencing of a genome. You can provide a pool of reads and the script will use the first 250 events to map a read and copy those mapping within the desired area to a specified folder.

Help is available by typing:

    C:\path\to\RUscripts\ReadUntil\python test_gReadUntil.py -h

which will output:

    usage: test_gReadUntil.py [-h] -fasta FASTA -targets [TARGETS [TARGETS ...]]
                          -procs PROCS -m TEMP_MODEL [-log LOGFILE] -w
                          WATCHDIR [-o OUTPUT_FOLDER] [-v] [-ver]

    real_read_until: A program providing read until with the Oxford Nanopore
    minION device. This program will ultimately be driven by minoTour to enable
    selective remote sequencing. This program is heavily based on original code
    generously provided by Oxford Nanopore Technologies.

    optional arguments:
    -h, --help            show this help message and exit
    -fasta FASTA, --reference_fasta_file FASTA
                        The fasta format file describing the reference
                        sequence for your organism.
    -targets [TARGETS [TARGETS ...]]
                        Positional IDs to enrich for in the form seqid:start-
                        stop . Can be space seperated eg: J02459:10000-15000
                        J02459:35000-40000
    -procs PROCS, --proc_num PROCS
                        The number of processors to run this on.
    -m TEMP_MODEL, --model TEMP_MODEL
                        The appropriate template model file to use
    -log LOGFILE, --log-file LOGFILE
                        The name of the log file that data will be written to
                        regarding the decision made by this program to process
                        read until.
    -w WATCHDIR, --watch-dir WATCHDIR
                        The path to the folder containing the downloads
                        directory with fast5 reads to analyse - e.g.
                        C:\data\minion\downloads (for windows).
    -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                        Path to a folder to symbolically place reads
                        representing match and not match.
    -v, --verbose-true    Print detailed messages while processing files.
    -ver, --version       show program's version number and exit

A typical command line to select reads mapping from 10-15kb in the lambda genome would be:

    C:\path\to\RUscripts\ReadUntil\python test_gReadUntil.py -fasta ../J02459.fasta -targets J02459:10000-15000 -procs 4  -m ../template_r7.3_e6_70bps_6mer_6.model -w ../RUtestset/ -o test

This would give the following output:

    ***********************************************************************************************
    **** This code will open a collection of reads and simulate read until on them. It will    ****
    **** copy reads into a secondary folder for subsequent processing by another analysis      ****
    **** package.                                                                              ****
    ***********************************************************************************************
    processing the reference fasta.
    J02459:10000-15000
    ID J02459
    length 48502
    1 J02459:10000-15000 10000 15000 J02459
    We want to extract this chunk J02459_1
    ID J02459_1
    length 5000
    FORWARD STRAND
    REVERSE STRAND
    ../RUtestset/llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_ch10_file49_strand.fast5 No Match
    ../RUtestset/llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_ch10_file69_strand.fast5 Sequence Found
    .
    .
    .
    ../RUtestset/llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_ch8_file58_strand.fast5 Sequence Found
    ../RUtestset/llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_ch9_file100_strand.fast5 No Match
    ../RUtestset/llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_ch8_file11_strand.fast5 No Match

Sequence found indicates that the read is derived from the desired region. No Match indicates that the read is from another region.


# Real Read Until

Important Note:

Running read until will influence the behaviour of your flow cell and change the output of your sequencing experiment. You are strongly advised to run simulations of read until prior to running on a live flow cell. The code as presented here is a demonstration of read until and one implementation. Users run this code entirely at their own risk.


We provide the two scripts used to generate the data presented in "Real time selective sequencing using nanopore technology".

These scripts absolutely require the read until API from Oxford Nanopore.

# Read Until API

This version of Read Until requires that the user runs a command window as Administrator.

The read until API changes depending on the version of minKNOW in use. The read until API for minKNOW versions 0.48.1.3 - 0.51.1.51 is available on request from Oxford Nanopore and the code provided here is compatible with Git Commit

    50ee65c1404870a2e04cb6dc211dcdeeddd3fd4f

and minKNOW versions <= 0.51.1.51

Note that this is not the latest commit. The latest version as of 30th March 2016 does not function as anticipated. It will work for simulation runs with ws_event_sample, but minKNOW itself does not function correctly with the API. This is actively being resolved by ONT. Thus, testing and simulation is possible with our code and later versions of minKNOW and ws_event_sampler but these scripts WILL NOT WORK AT THIS TIME with minKNOW > 0.51.1.51

# Installing the API

Copy the following files and folders (and there contents) into the RUscripts\ReadUntil folder:

    read_until.py
	event_sampler.thrift
    event_sampler/*
    example.py

# Running ws_event_sampler

Read until works by streaming data from minKNOW during a run to an external analysis pipeline and receiving instructions to reject reads as required. The program streaming the data is called ws_event_sampler.exe and this is found in the folder

    C:\grouper\binaries\bin

which is available after installing minKNOW.

Running ws_event_sampler is done as follows:

    C:\grouper\binaries>bin\ws_event_sampler.exe

*** NOTE: The command itself is bin\ws_event_sampler.exe so the command is executed from the folder binaries, not from within the folder bin. Executing the file from within the bin folder will result in the program crashing with a message about being unable to open a config file.

*** NOTE: ws_event_sampler provides no visible indication that it is running other than blocking the command line. If it fails, it will terminate - sometimes with an error message, but often not. Support for ws_event_sampler is available direct from ONT. We hope we have provided sufficient guidance for its use in the context of our scripts here.

ws_event_sampler can be run in simulation mode to facilitate the development of analysis scripts. When running in simulation mode, ws_event_sampler provides reads driven from a model file stored locally. This same model file must be used to analyse data when in simulation mode. When streaming live data direct from minKnow it is essential to use the correct model for the chemistry and pore type in use.

The read until API provides comprehensive instructions on how to run ws_event_sampler. Below we provide worked examples for our scripts.

To test that ws_event_sampler is functioning correctly, do the following:

1) Open a command window as administrator - to do this, right click on the cmd icon and select "Run as Administrator"

2) In this window navigate to C:\grouper\binaries

3) To see help for ws_event_sampler type:

    C:\grouper\binaries>bin\ws_event_sampler.exe -h

which will output:

    Usage: ws_event_sampler [Options]

    Common options:
      -h [ --help ]                         Produce help message
      -c [ --config ] arg (=./conf/global.conf)
                                            Application config file.

    Specific options:
      -p [ --port ] arg                     Port to listen for requests. If this is
                                            not specified it is taken from the
                                            ws_event_sampler_port field in the
                                            config file.
      -s [ --sim ]                          Start in standalone simulation mode.
      --sim-channels arg                    Number of channels simulation uses.
                                            E.g. use a low number when debugging to
                                            make the pushed data manageable. This
                                            overrides the channel_count setting in
                                            the conf file but only during
                                            simulation.
      --sim-model arg                       Path to model file to be used by
                                            simulation when generating events.
      --sim-fasta arg                       Path to fasta file from which to draw
                                            fragments to be used by simulation.
      --sim-delay arg (=0)                  Seconds delay before starting
                                            simulation. By default the simulation
                                            starts as soon as this server starts
                                            up.
      --sim-fragment-length arg             Random lengths are bounded, approximate
                                            normal distributions, picked by rolling
                                            lots of dice. You need to specify die
                                            sides, die count and a base amount. For
                                            example, the default setting of '5 250
                                            500' will generate reads with mean
                                            length 1000 and bounds [500,1500].
      --sim-log arg                         Filename to use for csv-format log of
                                            the simulation. Will overwrite if
                                            already exists.

To set up a simulation on a specific port type:

    C:\grouper\binaries>bin\ws_event_sampler.exe -s --sim-channels 100 --sim-fasta C:\path\to\RUscripts\J02459.fasta --sim-log log.txt

This will establish the data stream on port 12345, simulating 100 channels with a read distribution as described derived from the lambda sequence we have provided here. ws_event_sampler in simulation mode will write out a log file to log.txt enabling tracking of events.

In the example.py script, edit the line which states:

    host="ws://localhost:<port>"

to

    host="ws://localhost:12345"

Note - you can specify a different port value if you wish.

Then from another command window execute the examply.py script (note this is supplied by the ONT API)

    C:\path\to\RUscripts\ReadUntil>python example.py

This should output something like (depending on the version of the API in use):

    Client connection started. Beginning unblock loop...
    channel_name= 344  read_number= 12  events_in_sample= 200  first event: start= 1104931.0  mean= 71.1864189873
    channel_name= 345  read_number= 12  events_in_sample= 200  first event: start= 1111035.0  mean= 68.1394555562
    channel_name= 346  read_number= 12  events_in_sample= 200  first event: start= 1104312.0  mean= 67.1975140682
    channel_name= 347  read_number= 15  events_in_sample= 200  first event: start= 1203518.0  mean= 72.436242832
    channel_name= 340  read_number= 12  events_in_sample= 200  first event: start= 1104326.0  mean= 60.9232571465
    channel_name= 341  read_number= 12  events_in_sample= 200  first event: start= 1102520.0  mean= 69.496870232
    channel_name= 342  read_number= 12  events_in_sample= 200  first event: start= 1097286.0  mean= 65.2731245584
    channel_name= 343  read_number= 12  events_in_sample= 200  first event: start= 1102224.0  mean= 65.2680067692
    channel_name= 348  read_number= 12  events_in_sample= 200  first event: start= 1119138.0  mean= 77.3645669056
    channel_name= 349  read_number= 12  events_in_sample= 200  first event: start= 1092464.0  mean= 72.4693304753
    channel_name= 298  read_number= 12  events_in_sample= 200  first event: start= 1109404.0  mean= 66.3816659282
    channel_name= 299  read_number= 12  events_in_sample= 200  first event: start= 1103390.0  mean= 66.1452746381

If you see:

    Client connection started. Beginning unblock loop...
    Unhandled exception in thread started by <bound method ReadUntil.data_received_own_thread of <read_until.ReadUntil object at 0x00000000023AA
    FD0>>
    Traceback (most recent call last):
    File "C:\...\ReadUntil\read_until.py", line 110, in data_received_own_thread
    self.data_received(channels_update)
    File "example.py", line 46, in data_received
    for channel_id, data in channels.iteritems():
    AttributeError: 'NoneType' object has no attribute 'iteritems'
    Unhandled exception in thread started by <bound method ReadUntil.data_received_own_thread of <read_until.ReadUntil object at 0x00000000023AA
    FD0>>
    Traceback (most recent call last):
    File "C:\...\ReadUntil\read_until.py", line 110, in data_received_own_thread
    self.data_received(channels_update)
    File "example.py", line 46, in data_received
    for channel_id, data in channels.iteritems():
    AttributeError: 'NoneType' object has no attribute 'iteritems'
    Unhandled exception in thread started by <bound method ReadUntil.data_received_own_thread of <read_until.ReadUntil object at 0x00000000023AA
    FD0>>

Then you have likely got a clash between the API version we have specified and the version of minKNOW installed.


# gReadUntil.py

This script enables the selection of a specific region of a genome. We demonstrate this using Lambda.

This code optionally processes reads on a odd numbered channels - it allows all even numbered channels to be sequenced regardless of the target. This enables direct comparison of read until on one flow cell.

The default read type generated by ws_event_sampler is derived from E. coli. The current CPU bound version of DTW we are using will struggle on a genome of this size. Thus we suggest you run the simulator with (recalling to run this script in a seperate cmd window as an administrator):

    C:\grouper\binaries>bin\ws_event_sampler.exe -p 12345 -s --sim-channels 100 --sim-fragment-length 5 250 500 --sim-fasta C:\path\to\RUScripts\J02459.fasta --sim-log log.txt

You can obviously edit the number of channels and fragment lengths as you wish.

gReadUntil.py help statement:

    C:\path\to\RUScripts\ReadUntil\python gReadUntil.py -h

which outputs:

    usage: gReadUntil.py [-h] -fasta FASTA -targets [TARGETS [TARGETS ...]] -procs
                         PROCS -t TIME -m TEMP_MODEL [-ip IP] -p PORT
                         [-log LOGFILE] [-length LENGTH] [-skip] [-v] [-ver]

    gReadUntil.py: A program providing read until for genome sequences with the
    Oxford Nanopore minION device. This program will ultimately be driven by
    minoTour to enable selective remote sequencing. This program is partly based
    on original code generously provided by Oxford Nanopore Technologies.

    optional arguments:
      -h, --help            show this help message and exit
      -fasta FASTA, --reference_fasta_file FASTA
                            The fasta format file describing the reference
                            sequence for your organism.
      -targets [TARGETS [TARGETS ...]]
                            Positional IDs to enrich for in the form seqid:start-
                            stop . Can be space seperated eg: J02459:10000-15000
                            J02459:35000-40000
      -procs PROCS, --proc_num PROCS
                            The number of processors to run this on.
      -t TIME, --time TIME  This is an error catch for when we cannot keep up with
                            the rate of sequencing on the device. It takes a
                            finite amount of time to process through the all the
                            channels from the sequencer. If we cannot process
                            through the array quickly enough then we will 'fall
                            behind' and lose the ability to filter sequences.
                            Rather than do that we set a threshold after which we
                            allow the sequencing to complete naturally. The
                            default is 300 seconds which equates to 9kb of
                            sequencing at the standard rate.
      -m TEMP_MODEL, --model TEMP_MODEL
                            The appropriate template model file to use
      -ip IP, --ip-address IP
                            The IP address of the machine running minKNOW.
      -p PORT, --port PORT  The port that ws_event_sampler is running on.
      -log LOGFILE, --log-file LOGFILE
                            The name of the log file that data will be written to
                            regarding the decision made by this program to process
                            read until.
      -length LENGTH, --library-length LENGTH
                            Provide the average expected length of your library.
                            This offset will be applied to reads that are likely
                            to extend into your region of interest on either
                            strand.
      -skip, --skip_even    If set, this will allow all reads from even numbered
                            channels to be sequenced regardless of where they map.
                            This provides an internal control.
      -v, --verbose-true    Print detailed messages while processing files.
      -ver, --version       show program's version number and exit


An example run of gReadUntil.py:

    C:\path\to\RUScripts\ReadUntil\python gReadUntil.py -f ..\J02459.fasta -targets J02459:10000-20000 J02459:30000-45000 -procs 8 -t 100 -m C:\grouper\binaries\conf\synthesis\model.txt -ip 127.0.0.1 -p 12345 -length 1000

This configuration will map reads to the reference genome and select reads which map between 10-20 kb and 30-45 kb. We assume the library has an average length of 1kb (-length 1000) and so reads which start 500 bases 5' of the site of interest or 500 bases 3' (on the reverse strand) will be sequenced.

The script will output (note: You must choose Y or N to initiate read until):

    {'J02459': 48502}
    processing the reference fasta.
    ID J02459
    length 48502
    FORWARD STRAND
    REVERSE STRAND
    <type 'dict'>
    Failed to connect to minKNOW. Likely reasons include minKNOW not running, the wrong IP address for the minKNOW server or firewall issues. Two way control will not be possible with minKNOW. <urlopen error [Errno 61] Connection refused>
                                    `       
              ;`               ,;       
               :;;;;;;;;;;;;;;;,        
                  , .;;;;;, ,           
           @@@@@@@@ ;;;;;;; @@@@@@@@    
         @@@@@@@@@@# ;;;;; +@@@@@@@@@@  
        #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@
        .@@@@@    @@@@@@@@@@@    @@@@@:
         .@@@@`    @@@@@@@@@     @@@@,  
           '@@@@@@+ @@@@@@@ '@@@@@@+    
              ;@@@#  @@@@@. +@@@; inoTour read until routines.       
                   .;;;;;;;,            
                  ;;;.   .;;;`          
                  ;;       ;;`          
    Welcome to the .;;:     ,;;,          

    This script WILL implement read until.
    If you proceed it is at your own risk.

    Seriously - are you happy to proceed? Entering yes will make it your fault... [Y/n]
    ***********************************************************************************************
    ****        This version of the code will process reads regardless of channel.             ****
    ***********************************************************************************************
    1459369047.53
    Running Analysis
    Hanging around, waiting for the server...
    Running Analysis
    Hanging around, waiting for the server...
    Running Analysis
    Client connection started. Beginning unblock loop...
    2016-03-30 21:17:50 Unblocking  10
    2016-03-30 21:17:51 Unblocking  13
    2016-03-30 21:17:52 Unblocking  17
    2016-03-30 21:18:06 Unblocking  2
    2016-03-30 21:18:09 Unblocking  1

Here the code was executed prior to starting ws_event_sampler so the code waited for the server to start (Hanging around...). Once ws_event_sampler was started, the Connection is initiated and the code reports the number of channels it is rejecting reads from (Unblocking) at regular intervals.

The experiment can be terminated by Ctrl-C.

ws_event_sampler can then be terminated with Ctrl-C.

Inspection of the log.txt file from ws_event_sampler will allow validation of the code (an excerpt is shown):

    channel_name,seq_name,seq_index,seq_length,seq_head,read_number,read_start,read_event_count,read_reason
    .
    .
    .
    11,J02459,1019,1012,TGGTTGCCGACGGATG,4,217446,390,Unblock
    19,J02459,23886,1021,TTAAGTCTTCTTTCCC,4,220468,376,Unblock
    20,J02459,18434,988,GATAGCTGAAAACTGT,5,92238,1055,Natural
    22,J02459,3953,1023,TTCCCGGAATTACGCC,4,222261,364,Unblock
    24,J02459,5461,988,CAGATCACCGCAGCGG,4,216961,393,Unblock
    26,J02459,47375,1012,TCAGAATAAAACAATT,4,218657,383,Unblock
    27,J02459,32542,971,GGTTTTCATTGATGAT,5,92238,1038,Natural
    28,J02459,45131,994,TAATCGACCTTATTCC,4,222470,363,Unblock
    36,J02459,32223,979,GGAGTGATGTCGCGTT,5,97420,1046,Natural
    .
    .
    .
    4,J02459,9503,1020,GCCTTCCAGCCGGAGG,16,1035027,0,Stopped
    3,J02459,31374,1002,ACCAATTTCAGCCAGT,19,1136723,0,Stopped
    2,J02459,31534,969,GTCACCCACATGCTGT,20,1032682,0,Stopped
    1,J02459,12539,977,GCGGCGATGCTGACCG,16,1057685,0,Stopped

Note: The simulator only simulates reads from the forward strand. It is possible to force it to simulate both forward and reverse strands by supplying a fasta file containing both strands.

In the log file excerpt, Unblock = a rejected read, Natural = a read that has completed normally, Stopped = a read that finished prematurely due to ws_event_sampler being shut down.

By inspection of the above excerpt, all reads ending Naturally have start sites within the specified regions, whereas reads wich were rejected fall outside of the required regions.

# Running gReadUntil.py on a live run

Important Note:

Running read until will influence the behaviour of your flow cell and change the output of your sequencing experiment. You are strongly advised to run simulations of read until prior to running on a live flow cell. The code as presented here is a demonstration of read until and one implementation. Users run this code entirely at their own risk.


To run on live data, first change the configuration of ws_event_sampler:

    C:\grouper\binaries>bin\ws_event_sampler.exe -p 12345

This will now stream live data from whatever is being sequenced in minKNOW on port 12345

To run gReadUntil you MUST switch the model file to one appropriate for your chemistry and pore type. In our case:

    C:\path\to\RUScripts\ReadUntil\python gReadUntil.py -f ..\J02459.fasta -targets J02459:10000-20000 J02459:30000-45000 -procs 8 -t 100 -m ..\template_r7.3_e6_70bps_6mer_6.model -ip 127.0.0.1 -p 12345 -length 1000

Note: gReadUntil.py will write a message to the minKNOW messages window to tell you that we are remotely interacting with it.

Note: You should set the -length parameter to the expected average length of your library.

Note: As long as the appropriate ports are not blocked by your firewalls there is no reason why this code cannot be run on a separate computer - just configure the -ip statement accordingly.

Note: ws_event_sampler doesn't write out a log file when running on real data. gReadUntil can optionally output a log file in these cases with the -log option (e.g -log greaduntil.log). This will output the following information:

    Message:Channel,Read Number,Decision,Read Start Time,Ref ID,Distance,Orientation,Mapping Site
    INFO:24,2,REJ,23180.0,J02459,40.2154438853,F,3946
    INFO:25,2,REJ,23180.0,J02459,55.6327533467,F,6118
    INFO:26,2,SEQ,23180.0,J02459,76.1823545251,F,30637
    INFO:27,2,SEQ,23180.0,J02459,59.4201521716,F,34466
    INFO:28,2,REJ,23180.0,J02459,54.9151805152,F,47491
    INFO:20,2,SEQ,23180.0,J02459,52.824718186,F,30341
    INFO:29,2,REJ,23180.0,J02459,79.6781198595,F,20583
    INFO:21,2,REJ,23180.0,J02459,38.8591792471,F,7760
    INFO:23,2,REJ,23180.0,J02459,88.8619009677,F,20119
    INFO:22,2,REJ,23180.0,J02459,42.0277615893,F,5046
    INFO:4,2,REJ,23180.0,J02459,70.3307290757,F,26892
    INFO:8,2,REJ,23180.0,J02459,53.4937773307,F,5241
    INFO:59,2,SEQ,23180.0,J02459,54.0416642491,F,16105
    INFO:58,2,SEQ,23180.0,J02459,41.8579884151,F,12718
    INFO:55,2,SEQ,23180.0,J02459,40.0786231863,F,44787
    INFO:54,2,REJ,23180.0,J02459,57.7941719418,F,4125
    INFO:57,2,SEQ,23180.0,J02459,38.6372122551,F,10265

# aReadUntil.py

This script enables balanced sequencing of individual amplicons from a pool of amplicons. In its current form the script will try to balance the number of amplicons sequenced and so ensure uniform coverage. Alternatively, the script will enable the selection of individual amplicons, or even different coverage depths for each amplicon.

# Read Tracking

These methods are complicated by how the script tracks the final reads. We currently define three different optional counting points for a read:

1) read
    We track each channel individually. A read is assumed to have completed if another read is seen at that channel. This method assumes that every read start reported by read until will result in a true final sequenced read. By observation, we know that this assumption can fail. Presumably unseen minKNOW errors contribute to some of these failures. It is also possible that an observed read is not a true sequence or failed to generate sufficient data.
2) file
    Here we track files as they are written to disk. Files written to disk are paired up with entries tracked via read until and so matching data are conserved. When running a 2D library prep, a common failure are 1D reads - i.e a complement sequence is not detected. Thus it is possible that when counting files, 1D only reads might lower the yield.
3) 2d
    Here we track files in the same way, but further we inspect for the 2d flag within a file. This allows us to count the number of 2d potential reads (all prior to base calling). It should be noted that the current version of minKNOW at the time of writing - 0.51.1.62 - has some problems with event and hairpin detection, resulting in lower than expected 2d read counts.

# minKNOW Interaction

Because aReadUntil implements read tracking it needs to be able to 'see' reads as they are written to disk. This script also provides an implementation of the 'Run Until' concept - so it can optionally stop your sequencer when it believes a task is complete. It therefore communicates with minKNOW and sends messages to the Messages window on minKNOW. Because of this it expects minKNOW to be running and if the software cannot communicate with minKNOW it will stop with an error message. It makes this requirement even when running from the simulator.

# Running with metrichor/basecallers

aReadUntil interacts with read files. Metrichor moves files as they are basecalled from the folder they are in to a second folder. To avoid conflict, we have made aReadUntil move files as well. You should therefore point metrichor at this folder, not the usual output folder. The order of events is:

    1) minKnow writes file to folder (typically C:\data\reads)
    2) aReadUntil observes read, processes it and moves it to a 'done' folder ( C:\data\reads\done)
    3) metrichor reads files from the done folder and moves them to 'uploaded' and 'downloaded' as appropriate (C:\data\reads\done\uploaded and C:\data\reads\done\downloads)


To see the help message from aReadUntil type:

    C:\path\to\RUScripts\ReadUntil\>python aReadUntil.py -h

which will output:

    usage: aReadUntil.py [-h] -fasta FASTA [-c] -ids IDS [-d DEPTH]
                         [-cd CUSTOMDEPTH] -e DEPTHERROR -procs PROCS -t TIME -m
                         TEMP_MODEL -g GOAL [-precision] [-seq SPEED] [-i] [-s]
                         [-ip IP] -p PORT [-wt WRITETIME] -w WATCHDIR
                         [-log LOGFILE] [-v] [-v2] [-sim] [-ver]

    aReadUntil: A program providing read until with the Oxford Nanopore minION
    device. This program will ultimately be driven by minoTour to enable selective
    remote sequencing. This program is based on original code generously provided
    by Oxford Nanopore Technologies. Note that whilst some parameters can be set
    via a config file, the explicit parameters to stop a run (-s), prevent read
    until working (-i) and switch to presicion mode (-precison) can only be set
    via the command line. Args that start with '--' (eg. --reference-fasta-file)
    can also be set in a config file
    (/Users/mattloose/Dropbox/fixes/RUscripts/UPDATES/ReadUntil/aReadUntil.config
    or ) by using .ini or .yaml-style syntax (eg. reference-fasta-file=value). If
    an arg is specified in more than one place, then command-line values override
    config file values which override defaults.

    optional arguments:
      -h, --help            show this help message and exit
      -fasta FASTA, --reference-fasta-file FASTA
                            The fasta format file describing the reference
                            sequence for your organism.
      -c, --custom-genome   This will use a reduced search space genome to match
                            against.
      -ids IDS, --ids IDS   A file containing a list of amplicon positions defined
                            for the reference sequence. 1 amplicon per line in the
                            format fasta_sequence_name:start-stop e.g
                            EM_079517:27-1938
      -d DEPTH, --depth DEPTH
                            The desired coverage depth for each amplicon. Note
                            this is unlikely to be achieved for each amplicon and
                            should probably be an overestimate of the minimum
                            coverage required.
      -cd CUSTOMDEPTH, --custom-depth CUSTOMDEPTH
                            A comma seperated list of custom depths for each
                            amplicon. You must provide a coverage depth for each
                            amplicon in the order they are presented in the ids
                            file.
      -e DEPTHERROR, --error DEPTHERROR
                            Set an error range for coverage depth.
      -procs PROCS, --processor-number PROCS
                            The number of processors to run this on.
      -t TIME, --time TIME  This is an error catch for when we cannot keep up with
                            the rate of sequencing on the device. It takes a
                            finite amount of time to process through the all the
                            channels from the sequencer. If we cannot process
                            through the array quickly enough then we will 'fall
                            behind' and lose the ability to filter sequences.
                            Rather than do that we set a threshold after which we
                            allow the sequencing to complete naturally. The
                            default is 300 seconds which equates to 9kb of
                            sequencing at the standard rate.
      -m TEMP_MODEL, --tempalte-model TEMP_MODEL
                            The appropriate template model file to use
      -g GOAL, --goal GOAL  The measure by which reads will be counted - either
                            based on the presence of files ( -g file) or potential
                            2D files generated (-g 2d) or new reads generated ( -g
                            read )
      -precision            This option will attempt to obtain exactly the number
                            of reads required per amplicon. It is provided as a
                            novelty to illustrate the theoretical level of control
                            of the device. In reality it will slow down the time
                            taken to reach a specific goal due to the possibility
                            of reads failing and the delay in writing true reads
                            to disk.
      -seq SPEED, --seq-speed SPEED
                            This is the assumed sequencing speed. The default is
                            set at 30b/s (the speed of the simulator). This should
                            be configured to the appropriate value for your
                            chemistry.
      -i                    This will prevent read until from working but will
                            otherwise report what is happening in the sequencer.
      -s                    This will enable read until to stop your sequencing
                            when it is complete.
      -ip IP, --ip-address IP
                            The IP address of the minKNOW machine.
      -p PORT, --port PORT  The port that ws_event_sampler is running on.
      -wt WRITETIME, --write-time WRITETIME
                            If you are automatically stopping the minKNOW run, the
                            stop command will wait n seconds after the last read
                            has completed to ensure all files are written. Default
                            value is 15.
      -w WATCHDIR, --watch-dir WATCHDIR
                            The path to the folder containing the downloads
                            directory with fast5 reads to analyse - e.g.
                            C:\data\minion\downloads (for windows). This folder
                            must already exist. The script will not create it for
                            you. This is to prevent the wrong folder being
                            monitored for files which would disrupt read until.
      -log LOGFILE, --log-file LOGFILE
                            The name of the log file that data will be written to
                            regarding the decision made by this program to process
                            read until.
      -v, --verbose-detail  Print more detailed coverage info.
      -v2, --verbose-true   Print detailed messages while processing files.
      -sim, --sim           This action will write artificial fast5 files to a
                            folder for testing purposes.
      -ver, --version       show program's version number and exit

There are a number of important parameters here.

-i > Using this option, read until will actually be disabled. This allows the user to monitor what is happening on a sequencer via the read until api without ever actually rejecting reads. All other elements of the code (including squiggle matching) function. This allows for testing of the code without affecting sequencing.

-s > This option enables 'Run Until' - so the script will switch off your sequencer when it determines the goal has been reached.

-wt > Used in conjunction with 'Run Until' this option will allow the sequencer to continue running for a period of time after the goal has been reached. If using file or 2d tracking the default value of 15 seconds is sufficient. For read tracking it may be desirable to set this value to 180 or greater. This parameter can be tuned according to the rate at which minKNOW writes reads to disk.

-seq > The various chemistries used for nanopore sequencing operate at different speeds. This script tracks the number of reads currently being sequenced. To ensure that we don't end up with a read which never finishes being counted as 'currently sequenced' we calculate the anticipated length of the amplicons and with knowledge of the sequencing speed allow reads to 'time-out' if they have been considered active for too long. For the simulator, the default value of 30 b/s is correct, but this should be specified for the chemistry in use.  

-c  > This option will remove sequence from the center and 3'end of each strand of each amplicon from the reference sequence being matched too. This reduces the search space and so speeds up the rate of squiggle matching.

-precision > This option will attempt to limit the over sequencing of reads. It often adds time to a sequencing run rather than reducing it. This is a consequence of the observation that not all reads that begin to be sequenced will end up being written to disk for a variety of failure reasons. We do not suggest using this option at this time.

-sim  > This option will write out artificial fast5 files derived from the read until data stream to fully enable testing of the scripts. Reads are automatically written to the specified watchdir (set with the -w flag)


# Simulating an amplicon sequencing run.

To simulate amplicon sequencing we provide a workaround to enable the ws_event_sampler to stream individual amplicons. This requires a specially engineered fasta file, lambda_amplicons.fasta. This file contains amplicons in both forward and reverse complement orientation. Of the 11 amplicons, all odd numbered amplicons are present at 1x (1,3,5,7,9,11). Even numbered amplicons are present at different concentrations to simulate an uneven library.

To run ws_event_sampler with this file in a cmd window with administrator privileges enter:

    C:\grouper\binaries>bin\ws_event_sampler.exe -p 12345 -s --sim-channels 512 --sim-fragment-length 1 1 1900 --sim-fasta C:\path\to\ReadUntil\lambda_amplicons.fasta --sim-log amplicon-sim.log

Note: You can vary the number of channels as you feel appropriate for your computer. For the sim fragment length, we specify a length of 1900 bases per read. This forces the fragment to include the entire lenght of each amplicon with some flex around the start site, approximating what you might see on a real run.

To run an appropriate command for aReadUntil.py enter this command in another cmd window:

        C:\Users\plzmwl\Desktop\RUscripts\UPDATES\ReadUntil>python aReadUntil.py -fasta ..\J02459.fasta -ids ..\lambda_amplicons.txt -procs 8 -c -t 40 -m C:\grouper\binaries\conf\synthesis\model.txt -g 2d -seq 30 -ip 127.0.0.1 -p 12345 -w test  -s -sim -d 50 -e 0 -i


This command will use the lambda reference (-fasta ..\J02459.fasta) and the amplicon definitions (-ids ..\lambda_amplicons.txt) and match reads across 8 processors (-procs 8).

To optimise the experiment, we are tuning the reference sequence to omit regions that are unlikely to appear in the reads using the -c flag. This reduces the search time, but will misidentify reads that do not start within the first appoximately 400 bases of a specific amplicon.

The -t flag is set to 20 seconds. If a read has not been processed within 40 seconds of its start time it is automatically rejected. This is to ensure that a long queue of reads does not build up waiting for processing on slow systems. In the event of a timeout, the message:

    Read timeout

Will be written to the screen. This is a warning to indicate that the script is not keeping up with the run. It would be advisable to consider upping the number of processors available or reducing the scale of the experiment if possible. One or two Read Timeouts are reasonable and they will reduce as the flow cell decays, but this should be monitored. One way to increase the processing power available is to run this code from a second server. See the notes on this below.

For running with the simulator, the correct model must be used (currently as shown the model.txt file).

This run is looking at 2d reads (-g 2d) and assumes that sequencing is happening at approximately 30 b/s (-seq 30). We assume you are running on the same computer as minKNOW so the -ip is 127.0.0.1 and we have set ws_event_sampler running on port 12345 so -p 12345. The script is monitoring a folder called  "test" for the reads. As this is a simulation, the -sim flag is set and the script will write reads itself to the test folder.

The aim of this experiment is to sequence each amplicon at a depth of 50x (-d 50). The -e flag here is set to 0 and thus all amplicons must reach 50x. If -e is set to 5, then a value of 45x would be acceptable.

Finally for this first test we will not actually send any rejection messages (-i). This allows us to monitor the reads being produced without any interference from read until.


This command will run through various checks on parameters entered and request you confirm that you wish to run read until on your system.

    minoTour software is monitoring read until and will alert you when you have reached your target of 50x coverage on each amplicon.
    {'J02459': 48502}
    Mean amplicon length: 1929
    Autocalculated time threshold(s): 321
                                    `
              ;`               ,;
               :;;;;;;;;;;;;;;;,
                  , .;;;;;, ,
           @@@@@@@@ ;;;;;;; @@@@@@@@
         @@@@@@@@@@# ;;;;; +@@@@@@@@@@
        #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@
        .@@@@@    @@@@@@@@@@@    @@@@@:
         .@@@@`    @@@@@@@@@     @@@@,
           '@@@@@@+ @@@@@@@ '@@@@@@+
              ;@@@#  @@@@@. +@@@; inoTour read until routines.
                   .;;;;;;;,
                  ;;;.   .;;;`
                  ;;       ;;`
    Welcome to the .;;:     ,;;,

    This script will not implement read until. It will just report whatever is happening via the read unitl API.
    Are you happy to proceed? [Y/n]

This output tells you exactly what the script will do. You can agree or disagree. Disagreeing will exit the script.

The script will also send a message to the messages window of minKNOW:

    10:30:34: minoTour software is monitoring read until on this version of minKNOW and will alert you when you have reached your target of 50x coverage on each amplicon.

Agreeing will produce the following:

    Sample Name is:  sample_id
    1459697305.66
    Running Analysis
    2016-04-03 16:28:25 CACHED: 0 PROCESSED: 0
    Sun Apr  3 16:28:25 2016 : Obs: 0 Rej: Client connection started. Beginning unblock loop...
    0 Seq: 0 Done: 0 File: 0 2D: 0
    Obs Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    CurSDetails: {}
    2016-04-03 16:28:30 CACHED: 0 PROCESSED: 0
    Sun Apr  3 16:28:30 2016 : Obs: 22 Rej: 0 Seq: 22 Done: 0 File: 0 2D: 0
    Obs Details: {1: 1, 2: 0, 3: 8, 4: 3, 5: 1, 6: 4, 7: 3, 8: 1, 9: 0, 10: 1, 11: 1}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 1, 2: 0, 3: 8, 4: 3, 5: 1, 6: 4, 7: 3, 8: 1, 9: 0, 10: 1, 11: 1}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    CurSDetails: {1: 1, 3: 8, 4: 3, 5: 1, 6: 4, 7: 3, 8: 1, 10: 1, 11: 1}
    2016-04-03 16:28:35 CACHED: 233 PROCESSED: 0
    Sun Apr  3 16:28:35 2016 : Obs: 254 Rej: 0 Seq: 254 Done: 0 File: 0 2D: 0
    Obs Details: {1: 20, 2: 23, 3: 31, 4: 31, 5: 12, 6: 31, 7: 20, 8: 21, 9: 13, 10: 37, 11: 15}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 20, 2: 23, 3: 31, 4: 31, 5: 12, 6: 31, 7: 20, 8: 21, 9: 13, 10: 37, 11: 15}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    CurSDetails: {1: 20, 2: 23, 3: 31, 4: 31, 5: 12, 6: 31, 7: 20, 8: 21, 9: 13, 10: 37, 11: 15}

These text blocks will update every 5 seconds with a summary of what we are seeing via read until.

To briefly explain:

    2016-04-03 16:28:35 CACHED: 233 PROCESSED: 0

Date and time stamp. Cached and processed refers to the number of fast5 files that have been seen CACHED and then as they are analysed and matched to the original read until observation they move to PROCESSED.

    Sun Apr  3 16:28:35 2016 : Obs: 254 Rej: 0 Seq: 254 Done: 0 File: 0 2D: 0

This line repeats the time stamp. It then reports the number of read starts seen (Obs), the number of reads that have been rejected by the script (Rej), the number of reads sequenced (Seq). The 'Done', 'File' and '2D' values refer to total read counts at different stages of the process. Done reads are those where a subsequent read has been observed at that specific channel and so we assume the previous read completed. File are those reads that have been matched with a fast5 file written to disk. 2D identifies the 2D subset of those reads.

    Obs Details: {1: 20, 2: 23, 3: 31, 4: 31, 5: 12, 6: 31, 7: 20, 8: 21, 9: 13, 10: 37, 11: 15}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 20, 2: 23, 3: 31, 4: 31, 5: 12, 6: 31, 7: 20, 8: 21, 9: 13, 10: 37, 11: 15}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}


These lines provide a summary of the counts for each specific amplicon type in each of the measured read categories.

    CurSDetails: {1: 20, 2: 23, 3: 31, 4: 31, 5: 12, 6: 31, 7: 20, 8: 21, 9: 13, 10: 37, 11: 15}

The final line provides a count of the reads currently being sequenced according to the script.

After a period of time:

    2016-04-03 16:29:11 CACHED: 169 PROCESSED: 830
    Sun Apr  3 16:29:11 2016 : Obs: 999 Rej: 0 Seq: 998 Done: 487 File: 830 2D: 830
    Obs Details: {1: 64, 2: 117, 3: 94, 4: 134, 5: 60, 6: 135, 7: 70, 8: 82, 9: 74, 10: 110, 11: 59}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 64, 2: 117, 3: 94, 4: 134, 5: 60, 6: 134, 7: 70, 8: 82, 9: 74, 10: 110, 11: 59}
    DoneDetails: {1: 30, 2: 57, 3: 53, 4: 68, 5: 30, 6: 60, 7: 35, 8: 39, 9: 26, 10: 64, 11: 25}
    FileDetails: {1: 52, 2: 103, 3: 79, 4: 112, 5: 55, 6: 106, 7: 59, 8: 64, 9: 54, 10: 97, 11: 49}
    2D  Details: {1: 52, 2: 103, 3: 79, 4: 112, 5: 55, 6: 106, 7: 59, 8: 64, 9: 54, 10: 97, 11: 49}
    CurSDetails: {1: 34, 2: 60, 3: 41, 4: 66, 5: 30, 6: 75, 7: 35, 8: 43, 9: 48, 10: 46, 11: 34}
    You have reached your goal now. Sequencing should be stopped now!
    You have reached your goal now. Sequencing should be stopped now!
    You have reached your goal now. Sequencing should be stopped now!
    terminating sub-processes....


To end the run hit ctrl-C . (Note we have tried to catch all errors, but you may see spurious errors at this point. They will not affect the sequencing. We are continuing to work on this!)

The run is now complete and each amplicon has at least 50 x coverage of 2D reads. As read until has not been enabled here, no reads were rejected (Rej Details: all 0). This specific run took only 3 minutes 48 seconds (although numerous caveats have to be applied to this calculation - it is only a proxy for sequencing speed. For example, reads are only simulated by the event sampler in 1D thus each read would in reality take twice as long to generate.)

Now we will repeat this run but enable read until:

    C:\Users\plzmwl\Desktop\RUscripts\UPDATES\ReadUntil>python aReadUntil.py -fasta ..\J02459.fasta -ids ..\lambda_amplicons.txt -procs 8 -c -t 40 -m C:\grouper\binaries\conf\synthesis\model.txt -g 2d -seq 30 -ip 127.0.0.1 -p 12345 -w test  -s -sim -d 50 -e 0

Which will output:

    minoTour software is implementing read until and will send a stop sequencing message when the run is complete (defined by 50x coverage on each amplicon) AND SHUT DOWN minKNOW!!!!
    {'J02459': 48502}
    Mean amplicon length: 1929
    Autocalculated time threshold(s): 321
                                    `       
              ;`               ,;       
               :;;;;;;;;;;;;;;;,        
                  , .;;;;;, ,           
           @@@@@@@@ ;;;;;;; @@@@@@@@    
         @@@@@@@@@@# ;;;;; +@@@@@@@@@@  
        #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@
        .@@@@@    @@@@@@@@@@@    @@@@@:
         .@@@@`    @@@@@@@@@     @@@@,  
           '@@@@@@+ @@@@@@@ '@@@@@@+    
              ;@@@#  @@@@@. +@@@; inoTour read until routines.       
                   .;;;;;;;,            
                  ;;;.   .;;;`          
                  ;;       ;;`          
    Welcome to the .;;:     ,;;,          

    This script WILL implement read until.
    If you proceed it is at your own risk.


    Seriously - are you happy to proceed? Entering yes will make it your fault... [Y/n] Y
    Sample Name is:  sample_id
    1459697094.12
    Running Analysis
    2016-04-03 16:24:54 CACHED: 0 PROCESSED: 0
    Sun Apr  3 16:24:54 2016 : Obs: 0 Rej: 0 Seq: Client connection started. Beginning unblock loop...
    0 Done: 0 File: 0 2D: 0
    Obs Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    CurSDetails: {}



Initially this looks similar to the first run.

Over time these messages will appear:

    2016-04-03 16:25:26 CACHED: 36 PROCESSED: 512
    Sun Apr  3 16:25:26 2016 : Obs: 549 Rej: 23 Seq: 525 Done: 45 File: 512 2D: 512
    Obs Details: {1: 39, 2: 61, 3: 52, 4: 57, 5: 39, 6: 78, 7: 43, 8: 56, 9: 38, 10: 60, 11: 27}
    Rej Details: {1: 0, 2: 4, 3: 0, 4: 5, 5: 0, 6: 10, 7: 0, 8: 2, 9: 0, 10: 2, 11: 0}
    Seq Details: {1: 39, 2: 56, 3: 52, 4: 52, 5: 39, 6: 68, 7: 43, 8: 54, 9: 38, 10: 58, 11: 27}
    DoneDetails: {1: 0, 2: 11, 3: 4, 4: 4, 5: 4, 6: 4, 7: 3, 8: 4, 9: 3, 10: 5, 11: 3}
    FileDetails: {1: 37, 2: 57, 3: 48, 4: 52, 5: 38, 6: 68, 7: 42, 8: 54, 9: 33, 10: 58, 11: 25}
    2D  Details: {1: 37, 2: 57, 3: 48, 4: 52, 5: 38, 6: 68, 7: 42, 8: 54, 9: 33, 10: 58, 11: 25}
    CurSDetails: {1: 39, 2: 49, 3: 52, 4: 50, 5: 39, 6: 68, 7: 41, 8: 52, 9: 38, 10: 57, 11: 27}
    Unblocking  18
    Unblocking  20
    Unblocking  15
    Unblocking  18
    Unblocking  23
    Unblocking  16

The "Unblocking" numbers refer to the number of pores being unblocked at 1 second intervals. The other point to note is that the "Rej Details" line is now being populated. Observant readers will note that this number starts to increase before 2D reads have reached the 50x threshold. This is because aReadUntil takes in to account the reads it is expecting to see as well as those already written to disk so as not to overshoot.

Finally you will see this message:

    2016-04-03 16:25:45 CACHED: 344 PROCESSED: 744
    Sun Apr  3 16:25:45 2016 : Obs: 1094 Rej: 412 Seq: 681 Done: 512 File: 744 2D: 610
    Obs Details: {1: 79, 2: 126, 3: 89, 4: 125, 5: 90, 6: 145, 7: 69, 8: 117, 9: 73, 10: 117, 11: 64}
    Rej Details: {1: 13, 2: 69, 3: 21, 4: 73, 5: 22, 6: 77, 7: 6, 8: 63, 9: 9, 10: 59, 11: 0}
    Seq Details: {1: 66, 2: 56, 3: 68, 4: 52, 5: 68, 6: 68, 7: 63, 8: 54, 9: 64, 10: 58, 11: 64}
    DoneDetails: {1: 37, 2: 57, 3: 48, 4: 52, 5: 38, 6: 68, 7: 42, 8: 54, 9: 33, 10: 58, 11: 25}
    FileDetails: {1: 53, 2: 85, 3: 64, 4: 79, 5: 59, 6: 98, 7: 51, 8: 78, 9: 50, 10: 83, 11: 44}
    2D  Details: {1: 53, 2: 57, 3: 64, 4: 52, 5: 59, 6: 68, 7: 51, 8: 54, 9: 50, 10: 58, 11: 44}
    CurSDetails: {1: 56, 2: 33, 3: 56, 4: 33, 5: 52, 6: 50, 7: 46, 8: 35, 9: 54, 10: 42, 11: 55}
    Unblocking  13
    Unblocking  23
    Unblocking  17
    Sequencing complete. Now waiting for 15 seconds to ensure reads are written correctly to disk before stopping minKNOW.
    2016-04-03 16:25:59 CACHED: 262 PROCESSED: 1088
    Sun Apr  3 16:25:59 2016 : Obs: 1380 Rej: 694 Seq: 687 Done: 512 File: 1088 2D: 681
    Obs Details: {1: 97, 2: 157, 3: 117, 4: 158, 5: 110, 6: 184, 7: 93, 8: 145, 9: 93, 10: 143, 11: 85}
    Rej Details: {1: 31, 2: 100, 3: 49, 4: 106, 5: 42, 6: 116, 7: 30, 8: 91, 9: 29, 10: 85, 11: 15}
    Seq Details: {1: 66, 2: 56, 3: 68, 4: 52, 5: 68, 6: 68, 7: 63, 8: 54, 9: 64, 10: 58, 11: 70}
    DoneDetails: {1: 37, 2: 57, 3: 48, 4: 52, 5: 38, 6: 68, 7: 42, 8: 54, 9: 33, 10: 58, 11: 25}
    FileDetails: {1: 79, 2: 125, 3: 89, 4: 124, 5: 90, 6: 144, 7: 69, 8: 116, 9: 73, 10: 116, 11: 63}
    2D  Details: {1: 66, 2: 57, 3: 68, 4: 52, 5: 68, 6: 68, 7: 63, 8: 54, 9: 64, 10: 58, 11: 63}
    CurSDetails: {1: 55, 2: 33, 3: 56, 4: 33, 5: 52, 6: 48, 7: 45, 8: 34, 9: 54, 10: 41, 11: 61}
    The minoTours work is done. Sequencing will be stopped now!

Note that each amplicon has at least 50x coverage.

On these specific examples, the amplicons are already reasonably balanced and so no significant time saving is seen. If we try and induce a different patter, greater speed ups are seen. To illustrate we will use the custom depth option to sequence a 'ladder' pattern, going from 0x coverage on amplicon 1 to 500x coverage on amplicon 11.

To do this, type the following command:

    C:\Users\plzmwl\Desktop\RUscripts\UPDATES\ReadUntil>python aReadUntil.py -fasta ..\J02459.fasta -ids ..\lambda_amplicons.txt -procs 8 -c -t 40 -m C:\grouper\binaries\conf\synthesis\model.txt -g 2d -seq 30 -ip 127.0.0.1 -p 12345 -w test  -s -sim -cd 0,50,100,150,200,250,300,350,400,450,500 -e 0

We have substituted the -d flag for -cd.

This will output:

    minoTour software is implementing read until and will send a stop sequencing message when the run is complete (defined by 0,50,100,150,200,250,300,350,400,450,500x coverage on each amplicon) AND SHUT DOWN minKNOW!!!!
    {'J02459': 48502}
    Mean amplicon length: 1929
    Autocalculated time threshold(s): 321
                                    `       
              ;`               ,;       
               :;;;;;;;;;;;;;;;,        
                  , .;;;;;, ,           
           @@@@@@@@ ;;;;;;; @@@@@@@@    
         @@@@@@@@@@# ;;;;; +@@@@@@@@@@  
        #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@
        .@@@@@    @@@@@@@@@@@    @@@@@:
         .@@@@`    @@@@@@@@@     @@@@,  
           '@@@@@@+ @@@@@@@ '@@@@@@+    
              ;@@@#  @@@@@. +@@@; inoTour read until routines.       
                   .;;;;;;;,            
                  ;;;.   .;;;`          
                  ;;       ;;`          
    Welcome to the .;;:     ,;;,          

    This script WILL implement read until.
    If you proceed it is at your own risk.


    Seriously - are you happy to proceed? Entering yes will make it your fault... [Y/n] Y
    Sample Name is:  sample_id
    1459698615.17
    Running Analysis
    2016-04-03 16:50:15 CACHED: 0 PROCESSED: 0
    Sun Apr  3 16:50:15 2016 : Obs: 0 Rej: 0 Seq: 0 Done: 0 File: 0 2D: 0
    Obs Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Client connection started. Beginning unblock loop...
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    CurSDetails: {}

Part way through the run the selctive rejection can be seen (look at the Rej: Details):

    2016-04-03 16:42:07 CACHED: 113 PROCESSED: 4632
    Sun Apr  3 16:42:07 2016 : Obs: 4766 Rej: 1327 Seq: 3438 Done: 3078 File: 4632 2D: 3368
    Obs Details: {1: 350, 2: 552, 3: 389, 4: 526, 5: 334, 6: 656, 7: 319, 8: 423, 9: 424, 10: 471, 11: 325}
    Rej Details: {1: 349, 2: 442, 3: 182, 4: 211, 5: 0, 6: 146, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 0, 2: 111, 3: 207, 4: 314, 5: 334, 6: 511, 7: 319, 8: 423, 9: 424, 10: 470, 11: 325}
    DoneDetails: {1: 0, 2: 111, 3: 207, 4: 314, 5: 279, 6: 511, 7: 271, 8: 365, 9: 355, 10: 388, 11: 277}
    FileDetails: {1: 343, 2: 539, 3: 384, 4: 513, 5: 324, 6: 628, 7: 312, 8: 414, 9: 408, 10: 451, 11: 316}
    2D  Details: {1: 0, 2: 111, 3: 207, 4: 314, 5: 324, 6: 511, 7: 312, 8: 414, 9: 408, 10: 451, 11: 316}
    CurSDetails: {4: 2, 5: 69, 6: 23, 7: 71, 8: 91, 9: 87, 10: 105, 11: 64}
    Unblocking  16

And at the end of the run:

    Sequencing complete. Now waiting for 15 seconds to ensure reads are written correctly to disk before stopping minKNOW.
    2016-04-03 16:57:11 CACHED: 630 PROCESSED: 7231
    Sun Apr  3 16:57:11 2016 : Obs: 7860 Rej: 4959 Seq: 2900 Done: 2849 File: 7231 2D: 2883
    Obs Details: {1: 550, 2: 924, 3: 704, 4: 861, 5: 583, 6: 1069, 7: 516, 8: 706, 9: 645, 10: 760, 11: 546}
    Rej Details: {1: 550, 2: 853, 3: 591, 4: 699, 5: 362, 6: 810, 7: 208, 8: 346, 9: 237, 10: 287, 11: 20}
    Seq Details: {1: 0, 2: 71, 3: 111, 4: 161, 5: 221, 6: 259, 7: 309, 8: 360, 9: 408, 10: 474, 11: 526}
    DoneDetails: {1: 0, 2: 71, 3: 111, 4: 161, 5: 221, 6: 259, 7: 309, 8: 360, 9: 408, 10: 474, 11: 475}
    FileDetails: {1: 502, 2: 841, 3: 651, 4: 788, 5: 544, 6: 976, 7: 475, 8: 638, 9: 604, 10: 703, 11: 509}
    2D  Details: {1: 0, 2: 71, 3: 111, 4: 161, 5: 221, 6: 259, 7: 309, 8: 360, 9: 408, 10: 474, 11: 509}
    CurSDetails: {5: 3, 7: 69, 8: 32, 9: 71, 10: 97, 11: 240}
    The minoTours work is done. Sequencing will be stopped now!
    terminating sub-processes....
    Killed: 9

And without read until:

    python aReadUntil.py -fasta ../J02459.fasta -ids ../lambda_amplicons.txt -cd 0,50,100,150,200,250,300,350,400,450,500 -e 0 -procs 8 -t 20 -m /Applications/MinKNOW.app/Contents/Resources/conf/synthesis/model.txt -g 2d -ip 127.0.0.1 -w test -p 12345 -seq 30 -sim -c -s -i
    minoTour software is monitoring read until and will alert you when you have reached your target of 0,50,100,150,200,250,300,350,400,450,500x coverage on each amplicon.
    {'J02459': 48502}
    Mean amplicon length: 1929
    Autocalculated time threshold(s): 321
                                        `       
                  ;`               ,;       
                   :;;;;;;;;;;;;;;;,        
                      , .;;;;;, ,           
               @@@@@@@@ ;;;;;;; @@@@@@@@    
             @@@@@@@@@@# ;;;;; +@@@@@@@@@@  
            #@@@@@@`@@@@@ .;. @@@@@.@@@@@@@
            .@@@@@    @@@@@@@@@@@    @@@@@:
             .@@@@`    @@@@@@@@@     @@@@,  
               '@@@@@@+ @@@@@@@ '@@@@@@+    
                  ;@@@#  @@@@@. +@@@; inoTour read until routines.       
                       .;;;;;;;,            
                      ;;;.   .;;;`          
                      ;;       ;;`          
      Welcome to the .;;:     ,;;,          

    This script will not implement read until. It will just report whatever is happening via the read unitl API.
    Are you happy to proceed? [Y/n] Y
    Sample Name is:  sample_id
    1459699101.41
    Running Analysis
    2016-04-03 16:58:21 CACHED: 0 PROCESSED: 0
    Sun Apr  3 16:58:21 2016 : Obs: Client connection started. Beginning unblock loop...
    0 Rej: 0 Seq: 0 Done: 0 File: 0 2D: 0
    Obs Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    DoneDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    FileDetails: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    2D  Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    CurSDetails: {}

    2016-04-03 17:07:37 CACHED: 244 PROCESSED: 6626
    Sun Apr  3 17:07:37 2016 : Obs: 6872 Rej: 0 Seq: 6873 Done: 6365 File: 6626 2D: 6626
    Obs Details: {1: 508, 2: 759, 3: 556, 4: 755, 5: 494, 6: 944, 7: 494, 8: 613, 9: 579, 10: 664, 11: 506}
    Rej Details: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
    Seq Details: {1: 508, 2: 760, 3: 556, 4: 755, 5: 494, 6: 944, 7: 494, 8: 613, 9: 580, 10: 663, 11: 506}
    DoneDetails: {1: 466, 2: 697, 3: 512, 4: 696, 5: 463, 6: 871, 7: 458, 8: 570, 9: 539, 10: 620, 11: 473}
    FileDetails: {1: 488, 2: 728, 3: 536, 4: 725, 5: 474, 6: 909, 7: 479, 8: 594, 9: 556, 10: 645, 11: 492}
    2D  Details: {1: 488, 2: 728, 3: 536, 4: 725, 5: 474, 6: 909, 7: 479, 8: 594, 9: 556, 10: 645, 11: 492}
    CurSDetails: {1: 42, 2: 64, 3: 44, 4: 59, 5: 32, 6: 73, 7: 36, 8: 43, 9: 41, 10: 44, 11: 34}
    You have reached your goal now. Sequencing should be stopped now!
    You have reached your goal now. Sequencing should be stopped now!
    You have reached your goal now. Sequencing should be stopped now!
    terminating sub-processes....
    Killed: 9


# Running live read until amplicon sequencing.

Important Note:

Running read until will influence the behaviour of your flow cell and change the output of your sequencing experiment. You are strongly advised to run simulations of read until prior to running on a live flow cell. The code as presented here is a demonstration of read until and one implementation. Users run this code entirely at their own risk.

To run amplicon balancing on live data, first change the configuration of ws_event_sampler:

    C:\grouper\binaries>bin\ws_event_sampler.exe -p 12345

This will now stream live data from whatever is being sequenced in minKNOW on port 12345.

For the aReadUntil.py script it is essential that you:

1). Use the correct model file for your chemistry.
2). Point aReadUntil at the correct folder for seeing the raw reads as output by minKNOW (typically C:\data\reads). aReadUntil.py does not create this folder itself. It must already exist. If aReadUntil cannot see the reads, it cannot keep track of reads actually written to disk.
3). aReadUntil will move reads it has seen to a subfolder called 'done' within the reads folder (so typically C:\data\reads\done) and it is this folder that you should point metrichor too.

You can then run any of the script options above to try aReadUntil.

COMMON FAILURES.

1) ws_event_sampler crashes: If ws_event_sampler crashes, the read until scripts will no longer 'see' the events and thus the sequencer will sequence all reads regardless of where they map too. The scripts alert you to a dropped connection with the message "Hanging around waiting for server". You should check for this periodically. We hope that ws_event_sampler will become more stable over time.
2) Processor timeouts: If the load on the CPU is too high, reads will 'timeout'. This will result in a loss of selective sequencing.
3) Missing amplicons: If an amplicon is entirely missing from a sample, read until will never complete. It is wise to carefully observe a run and check this hasn't occurred.

To monitor a run without implementing read until, alerting when you have 40x coverage of each amplicon:

    C:\path\to\RUscripts\ReadUntil>python aReadUntil.py -fasta ..\J02459.fasta -ids ..\lambda_amplicons.txt -procs 8 -c -t 40 -m ..\template_r7.3_e6_70bps_6mer_6.model -g 2d -seq 30 -ip 127.0.0.1 -p 12345 -w test  -s -sim -d 40 -e 0 -i

To implement read until, alerting when you have 40x coverage of each amplicon:

    C:\path\to\RUscripts\ReadUntil>python aReadUntil.py -fasta ..\J02459.fasta -ids ..\lambda_amplicons.txt -procs 8 -c -t 40 -m ..\template_r7.3_e6_70bps_6mer_6.model -g 2d -seq 30 -ip 127.0.0.1 -p 12345 -w test  -s -sim -d 40 -e 0

#Running read until from a remote machine.

Running aReadUntil.py from a remote machine is perfectly possible. Ensure that you replace the ip address for the target machine (-ip). The most critical issue is that the aReadUntil.py script must be able to see the reads that are written to disk. You therefore need a way to synchronise reads from the machine running minKNOW to the target server.

One possible method is to use rsync in an infinitely looping shell script.

#Getting help with read until.

We are happy to help where we can. Please feel free to contact via twitter ({at}mattloose) or email matt.loose{@}nottingham.ac.uk . We will try to update this read me file whenever significant changes are made to the API or minKNOW versions that might affect read until. For issues with the API itself or the ws_event_sampler program, users will most likely be redirected to ONT.
