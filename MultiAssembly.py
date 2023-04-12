#!/usr/bin/python3
import subprocess, sys

output_dir = sys.argv[2]
read1 = sys.argv[3]
read2 = sys.argv[4]
phred_offset = sys.argv[5]
directory = sys.argv[6]


#!/usr/bin/python3
import subprocess, re

def SPADES(read1,read2,directory,spades_tag,phred_offset):
    output_dir = "assembly" 
    
    if "Assembly" not in spades_tag:
        print("Running spades.py")
        subprocess.run(["conda", "run", "-n", "ASSEMBLY", "spades.py", "-o", output_dir, "-1", read1, "-2", read2, "--meta","--phred-offset",phred_offset], cwd = directory)
        print("Spades.py finished. \n")
    
    return output_dir


def SubSampling(read1,read2,directory,sampleRate,sampleSeed, SkipTag): #Subsampling using Reformat.sh 
    
    read1WithNoFileFormat = re.search(r'(\w+)\.fastq',read1).groups()[0]
    read2WithNoFileFormat = re.search(r'(\w+)\.fastq',read2).groups()[0]
    read1Trimmed = read1WithNoFileFormat + "_trimmed.fastq"
    read2Trimmed = read2WithNoFileFormat + "_trimmed.fastq"
    if SkipTag == None:
        print("Subsampling reads using reformat.sh with samplerate =", sampleRate, "and sampleseed =", sampleSeed, "\n")
        subprocess.run(["conda","run","-n", "QC","reformat.sh","in=" + read1, "in2=" + read2, "out=" + read1Trimmed, "out2=" + read2Trimmed,"samplerate=" + str(sampleRate),"sampleseed=" + str(sampleSeed),"overwrite=true"], cwd = directory)
    print("Finished subsampling")
    return read1Trimmed, read2Trimmed


def MultiAssembly(read1, read2, directory, phred_offset, sampleRate, nrofassemblies, SkipTag):
    
    maxN50 = None
    maxseed = None
    #Initialisation
    contiglengthcutoff = 500
    if "Assembly" not in SkipTag:
        Skip = None
    
    if "Assembly" not in SkipTag:
        coverage = None

        print("\nFinding the optimal subsampling rate")
        print("Starting with samplingrate:",sampleRate,"\n")
        while coverage == None or coverage < 20 or coverage > 100:
            sampleSeed = nrofassemblies


            read1Trimmed, read2Trimmed = SubSampling(read1,read2,directory,sampleRate,sampleSeed, Skip)
            assemblydirectory = SPADES(read1Trimmed,read2Trimmed,directory, SkipTag, phred_offset)
            #trimmedContigs = contigTrimming(directory, assemblydirectory + "/contigs.fasta",contiglengthcutoff)
            coverage, coveragestatsfile = coverageFinderAverage(read1Trimmed,read2Trimmed,directory,trimmedContigs)
            subprocess.run(["rm","-rf",assemblydirectory], cwd = directory)
            subprocess.run(["rm",trimmedContigs], cwd = directory)
            subprocess.run(["rm",read1Trimmed], cwd = directory)
            subprocess.run(["rm",read2Trimmed], cwd = directory)
            #subprocess.run(["rm",coveragestatsfile], cwd = directory)
                

            if coverage < 20:
                sampleRate = sampleRate * (25 / coverage) #Aiming for x25 coverage
            elif coverage > 100:
                sampleRate = sampleRate * (80 / coverage) #Aiming for x80 coverage
            

            if sampleRate > 1:
                sampleRate = 1
                print(coverage)
                print(sampleRate)
                break
            print(coverage)
            print(sampleRate)
    
        print("Found best sampleRate:", sampleRate)
        for sampleSeed in range(1,int(nrofassemblies)+1):
            print("Assembling using SPADES.py", sampleSeed ,"out of", nrofassemblies, "...")
            read1Trimmed, read2Trimmed = SubSampling(read1,read2,directory,sampleRate,sampleSeed, Skip)
            assemblydirectory = SPADES(read1Trimmed,read2Trimmed,directory, SkipTag, phred_offset)
            N50val = N50(directory,assemblydirectory)
            if maxN50 is None or N50val > maxN50:
                maxseed = sampleSeed
                maxN50 = N50val
            print("\nMax N50:", maxN50)
            print("Best seed:", maxseed, "\n")
            subprocess.run(["rm","-rf",assemblydirectory],cwd = directory)
            subprocess.run(["rm",read1Trimmed],cwd = directory)
            subprocess.run(["rm",read2Trimmed],cwd = directory)
        

    print("Running assembly for the best subsampling seed...\n")
    read1Trimmed, read2Trimmed = SubSampling(read1,read2,directory,sampleRate, maxseed, Skip)
    assemblydirectory = SPADES(read1Trimmed,read2Trimmed,directory, SkipTag, phred_offset)
    trimmedContigs = contigTrimming(directory, assemblydirectory + "/contigs.fasta",contiglengthcutoff)
    coverage, coveragestatsfile = coverageFinderAverage(read1Trimmed,read2Trimmed,directory,trimmedContigs)
    print(coverage)
    print(sampleRate)
    print("Finished assembly for the best subsampling seed:", maxseed,"\n")
    
    return assemblydirectory, read1Trimmed, read2Trimmed


def N50(directory,assemblydirectory): #Calculating N50 using stats.sh from BBtools
    filename = "N50assemblystats.txt"
    subprocess.run(["conda","run","-n","QC","stats.sh","in=" + assemblydirectory + "/contigs.fasta",">",filename],cwd = directory)
    
    with open(filename,'r') as file:
        for line in file:
            if line.startswith("Main genome contig N/L50:"):
                nl50 = line.split()[4]
                reg_obj = re.search(r'(\w+)\/',nl50)
                N50 = float(reg_obj.groups()[0])
                print(N50)
    subprocess.run(["rm",filename],cwd = directory)
    print("N50 =", N50)
    return N50


def coverageFinderAverage(read1,read2,directory,contigfilepath):
    print("Finding maximum coverage among assemblies which are larger than half of the size of the largest contig.")
    contigs = contigfilepath
    coveragestats = "coveragestats.txt"
    subprocess.run(["conda","run","-n","QC","bbmap.sh","ref=" + contigs,"in=" + read1,"in2=" + read2,"out=coverage_mapping.sam","nodisk=t","fast=t","covstats="+coveragestats],cwd = directory)
    print("Finished")
    linecount = 0
    sumcoverage = 0
    longestcontiglength = None
    with open(coveragestats,'r') as covfile:
        for line in covfile:
            linesplit = line.split()
            if linecount == 1:
                print(line)
                longestcontiglength = float(linesplit[2])
                sumcoverage = float(linesplit[1])
            elif linecount > 1:
                if float(linesplit[2]) > longestcontiglength * 0.7:
                    print(line)
                    sumcoverage += float(linesplit[1])
                else:
                    break
                
            linecount += 1      
    averagecoverage = sumcoverage / linecount - 1
    return averagecoverage, coveragestats



