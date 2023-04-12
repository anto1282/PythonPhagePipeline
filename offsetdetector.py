#!/usr/bin/python3

import sys


def offsetDetector(read1,read2):
    maxASCII = None
    minASCII = None
    #print("Finding phred-offset")
    reads = [read1,read2]
    for read in reads:
        with open(read) as file:
            phredflag = False
            for line in file:
                if line.startswith("+"):
                    phredflag = True
                elif line.startswith("@"):
                    phredflag = False
                elif phredflag == True:
                    
                    for char in line.strip():
                        if maxASCII is None or char > maxASCII:
                            maxASCII = char
                        if minASCII is None or char < minASCII:
                            minASCII = char
                    
                    if minASCII < "@":
                        phredoffset = "33"
                        #print("Phred-Offset:", phredoffset)
                        return phredoffset
                    
                    if maxASCII > "K":
                        phredoffset = "64"
                        #print(phredoffset)
                        return phredoffset
                    
print(offsetDetector(sys.argv[1],sys.argv[2]))