#!/usr/bin/env python
import sys
import h5py
import configargparse

def get_model_location(hdf,strand):
    #print filepath

    for element in hdf:
        #print element
        for element2 in hdf[element]:
            #print element2
            for element3 in hdf[element][element2]:
                #print element3
                try:
                    for element4 in hdf[element][element2][element3]:
                        #if "model" in (element,element2,element3,element4):
                        if any("Model" in s for s in [element,element2,element3,element4]):
                            if any(strand in s for s in [element,element2,element3,element4]):
                                return element,element2,element3,element4
                except:
                    print ""




if __name__ == "__main__":

    parser = configargparse.ArgParser(description='getmodels: A program to extract model files from Nanopore reads. You must provide a read that has given a 2D read - i.e comes from the pass folder.')
    parser.add('-read', '--read', type=str, dest='read', required=True, default=None, help="Provide a read file to extract the current model from.")
    args = parser.parse_args()
    hdf = h5py.File(args.read, 'r')
    temp1,temp2,temp3,temp4 = get_model_location(hdf,"template")
    comp1,comp2,comp3,comp4 = get_model_location(hdf,"complement")

    print temp1,temp2,temp3,temp4
    print comp1,comp2,comp3,comp4

    file = open("template.model", "w")
    for thing in hdf[temp1][temp2][temp3][temp4]:
        writestring = str(thing[0])+"\t"+str(thing[1])+"\t"+str(thing[2])+"\n"
        file.write(writestring)
        #print thing[0],thing[1],thing[2]
    file.close()

    file = open("complement.model", "w")
    for thing in hdf[comp1][comp2][comp3][comp4]:
        writestring = str(thing[0])+"\t"+str(thing[1])+"\t"+str(thing[2])+"\n"
        file.write(writestring)
        #print thing[0],thing[1],thing[2]
    file.close()

    hdf.close()

    print "File Write Completed."
    print "Kmer length is:",len(str(thing[0]))
    print "File format is:"
    print "Kmer\tMean\tStandard Dev"
    sys.exit()
