#!/usr/bin/env python3

import sys


#Script that saves predicted viral contigs in one file and the predicted non-viral contigs 
cutoff = float(sys.argv[3])

#try: 
print("Running DeepVirExtractor") 
if cutoff > 1 or cutoff < 0:
    print("Cutoff must be between 0 and 1")
    raise ValueError

virusnames = set()

contigfile = sys.argv[2]
virusfile = "predicted_viruses.fasta"
nonvirusfile = "non_viral_assemblies.fasta"
predfile = sys.argv[1]

linecount = 0
print("\nCreating set of viral entries...")
with open(predfile,'r') as file: 
    for line in file:
        if linecount > 0:
            if float(line.split()[2]) > cutoff:
                virusnames.add(line.split()[0])
            
        linecount += 1

print(linecount - 1, "total entries.")
print(len(virusnames), "fasta entries passed cutoff.")

virusoutfile = open(virusfile,"w")
nonvirusoutfile = open(nonvirusfile,'w')

with open(contigfile, 'r') as file:
    virusflag = False
    nonvirusflag = False
    seqcount = 0
    for line in file:
        #print(virusflag)
        if line.startswith(">"):
            virusflag = False             
            nonvirusflag = False
        
        if virusflag == True:
            virusoutfile.write(line)
        elif nonvirusflag == True:
            nonvirusoutfile.write(line)
        elif line.startswith(">") and line[1:].strip().split()[0] in virusnames:
            virusflag = True
            virusoutfile.write(line)
            seqcount += 1
        else:
            nonvirusflag = True
            nonvirusoutfile.write(line)
        
        
    print(seqcount, "sequence entries written to output file:", virusfile)
    print("The remaining entries was written to the non-viral file:", nonvirusfile)
    
    
virusoutfile.close()
nonvirusoutfile.close()
""""
except IndexError as error:
    print("Command argument missing:", error)
    sys.exit()
except ValueError as error:
    print("Command argument missing:", error)
    sys.exit()
"""
