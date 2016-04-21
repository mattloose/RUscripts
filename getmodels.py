#!/usr/bin/env python
import sys,os
import h5py
import configargparse
sys.path.append("ReadUntil")
from ruutils import query_yes_no



__version__ = "1.0"
__date__ = "29th March 2016"


def get_model_location(hdf,strand,args):
    '''
    Function to identify the presence of a Model field in an hdf file.
    This should find the template or complement model file if they are present.
    '''
    for element in hdf:
        for element2 in hdf[element]:
            for element3 in hdf[element][element2]:
                try:
                    for element4 in hdf[element][element2][element3]:
                        if any("Model" in s for s in [element,element2,element3,element4]):
                            if any(strand in s for s in [element,element2,element3,element4]):
                                return element,element2,element3,element4
                except Exception, value:
                    if (args.debug is True):
                        print "*** DEBUG MESSAGE ***"
                        print "Looking for model."
                        print "Error : " , value
                        print "*** END DEBUG ***"
                    print hdf[element][element2][element3][element4]
                    print "Most likely a " + strand + " model hasn't been found in this read."

def get_model_type(hdf,strand,args):
    '''
    Function to identify the model type in an hdf file.
    '''
    for element in hdf:
        for element2 in hdf[element]:
            for element3 in hdf[element][element2]:
                try:
                    for element4 in hdf[element][element2][element3]:
                        if (any("general" in s for s in [element,element2,element3,element4])):
                            findval = "Basecall"
                            if str(hdf[element][element2][element3][element4]).lower().find(findval.lower()) != -1:
                                return hdf[element][element2][element3][element4].attrs['model_type']
                except Exception, value:
                    if (args.debug is True):
                        print "*** DEBUG MESSAGE ***"
                        print "Looking for model type information."
                        print "Error : " , value
                        print "*** END DEBUG ***"
                    print "Most likely a " + strand + " model definition hasn't been found in this read."




if __name__ == "__main__":

    parser = configargparse.ArgParser(
    description='getmodels: A program to extract model files from Nanopore reads. Ideally you should provide a read that has given a 2D read - i.e comes from the pass folder on a 2D basecalling run. For Read Until purposes, extraction of just a template model is sufficient. Model files will be appended with a number to identify the kmer length of the model.'
    )
    parser.add('-r', '--read', type=str, dest='read', help="Provide a read file to extract the current model from.")
    parser.add('-v', '--verbose-true', action='store_true', help="Print detailed messages while processing files to aid in debugging.", default=False, dest='debug')
    parser.add_argument('-ver', '--version', action='version',version=('%(prog)s version={version} date={date}').format(version=__version__,date=__date__))
    args = parser.parse_args()



    if (args.read is None):
        print "You must pass a single read file to this script using the -r option."
        print "*** This code will now exit. It has not written out any model files."
        sys.exit()

    print "Processing " +  args.read

    try:
        hdf = h5py.File(args.read, 'r')
    except:
        print "Sorry - this file appears to be corrupt. Please try another."
        print "*** This code will now exit. It has not written out any model files."
        sys.exit()



    modelstrand=("template","complement")

    for model in modelstrand:
        modeltype=get_model_type(hdf,model,args)
        #sys.exit()
        print "Looking for " + model + " model."

        try:
            temp1,temp2,temp3,temp4 = get_model_location(hdf,model,args)
            try:
                ####Check folderdict.py line 42
                lookuparray=dict()
                for entry in range(0,len(hdf[temp1][temp2][temp3][temp4].dtype.descr)):
                    lookuparray[hdf[temp1][temp2][temp3][temp4].dtype.descr[entry][0]]=entry
                kmerlen = str(len(hdf[temp1][temp2][temp3][temp4][0][lookuparray['kmer']]))
                if os.path.isfile(model+"_"+ modeltype + "_" + kmerlen +".model"):
                    if not query_yes_no("It looks like this model file already exists.\nOverwriting this file should cause no problems as long as your source read file contains the most recent model for your chemistry and pore type OR is the correct type for the files you are working with.\nAre you happy to proceed?"):
                        os._exit(1)
                file = open(model+"_"+ modeltype + "_" + kmerlen +".model", "w")
                file.write("kmer\tlevel_mean\tlevel_stdv\n")
                for kmerrow in range(0,len(hdf[temp1][temp2][temp3][temp4])):
                    writestring = str(hdf[temp1][temp2][temp3][temp4][kmerrow][lookuparray['kmer']])+"\t"+str(hdf[temp1][temp2][temp3][temp4][kmerrow][lookuparray['level_mean']])+"\t"+str(hdf[temp1][temp2][temp3][temp4][kmerrow][lookuparray['level_stdv']])+"\n"
                    file.write(writestring)
                file.close()
                print "***: Model Found"
                print model + " model file write completed."
                print "Kmer length is:",kmerlen
                print "File format is:"
                print "Kmer\tMean\tStandard Dev"
            except:
                print "An error has occurred writing the "+model+ " file to disk."
        except:
            print "No model file found for " + model + ". If this isn't intentional then please try another read."
    hdf.close()

    sys.exit()
