# AmpBalanceP
AmpBalance Parallel

Very simple little program to try and split up amplicons.

This should work on the default version of python as registered on an ONT PC. It does require the mlpy library to be installed:

pip install mlpy-3.5.0-cp27-none-win_amd64.whl

After that you should just be able to run:

python ampbalance_v2P.py -w /path/to/reads -o /path/to/outfolder -d 50 -p 8

Note that the full options are available with -h

Many options are now encoded in the ampW.config file for simplicity. 

usage: ampbalance_v2P.py [-h] -fasta FASTA -ids IDS -w WATCHDIR -o TARGETPATH
                         -d DEPTH -procs PROCS -l LENGTH [-v]

AmpBalance now outputs lines as follows:

{'DO': 1456, 1: 41, 2: 30, 3: 77, 4: 49, 5: 36, 6: 647, 7: 408, 8: 112, 9: 598, 10: 38, 11: 19, 'BF': 2, 'TF': 301, 'HF': 3505, 'NH': 1203}

{'DO': 1456, 1: 41, 2: 30, 3: 77, 4: 49, 5: 36, 6: 647, 7: 409, 8: 112, 9: 598, 10: 38, 11: 19, 'BF': 2, 'TF': 300, 'HF': 3506, 'NH': 1203}

{'DO': 1456, 1: 41, 2: 30, 3: 77, 4: 49, 5: 36, 6: 647, 7: 409, 8: 112, 9: 598, 10: 38, 11: 19, 'BF': 2, 'TF': 299, 'HF': 3506, 'NH': 1204}

{'DO': 1456, 1: 41, 2: 30, 3: 77, 4: 49, 5: 36, 6: 647, 7: 409, 8: 112, 9: 598, 10: 38, 11: 19, 'BF': 2, 'TF': 298, 'HF': 3506, 'NH': 1205}

{'DO': 1457, 1: 41, 2: 30, 3: 77, 4: 49, 5: 36, 6: 647, 7: 409, 8: 112, 9: 598, 10: 38, 11: 19, 'BF': 2, 'TF': 297, 'HF': 3506, 'NH': 1206}

{'DO': 1457, 1: 41, 2: 30, 3: 77, 4: 49, 5: 36, 6: 647, 7: 409, 8: 112, 9: 598, 10: 38, 11: 19, 'BF': 2, 'TF': 296, 'HF': 3506, 'NH': 1206}

Where:
    DO: Template and Complement Strands don't overlap on amplicons as expected. Read Discarded. Unlikely to generate good 2D.
    
    1...11: Amplicon numbers followed by current coverage counts satisfying highest stringency tests.
    
    BF: Bad Files - these irritations cause windows to crash due to funky memory errors so get skipped automagically.
    
    TF: Total Files - this counts down the number of files left to process.
    
    HF: Hairpin Found - this is the number of reads containing a hairpin - these are tested.
    
    NH: No Hairpin - these reads are ignored as no hairpin has been found.
