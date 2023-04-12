#!/usr/bin/python3
import subprocess, re, os, glob

def offsetDetector(read1,read2,directory):
    maxASCII = None
    minASCII = None
    print("Finding phred-offset")
    reads = [read1,read2]
    for read in reads:
        with open(directory + "/" + read) as file:
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
                        print("Phred-Offset:", phredoffset)
                        return phredoffset
                    
                    if maxASCII > "K":
                        phredoffset = "64"
                        print(phredoffset)
                        return phredoffset

def SPADES(read1,read2,directory,spades_tag,phred_offset):
    output_dir = "assembly" 
    
    if "Assembly" not in spades_tag:
        print("Running spades.py")
        subprocess.run(["conda", "run", "-n", "ASSEMBLY", "spades.py", "-o", output_dir, "-1", read1, "-2", read2, "--meta","--phred-offset",phred_offset], cwd = directory)
        print("Spades.py finished. \n")
    
    return output_dir

 subprocess.run(reformat.sh in=${r1} in2=${r2} out=${subsampled_reads[0]}_#5_subsampled.fastq out2=${subsampled_reads[1]}_#5_subsampled.fastq samplerate=${samplerate} sampleseed=5)
   
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
            trimmedContigs = contigTrimming(directory, assemblydirectory + "/contigs.fasta",contiglengthcutoff)
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


def contigTrimming(directory,Contigs_fasta, minLength=200):
    Contigs_trimmed = "contigs_trimmed.fasta"

    subprocess.run(["conda", "run", "-n","QC","reformat.sh","in="+Contigs_fasta, "out=" + Contigs_trimmed, "minlength="+str(minLength),"overwrite=True"], cwd = directory)
    return Contigs_trimmed


def DeepVirFinder(pathtoDeepVirFinder,assemblydirectory,threads,inFile,PredTag):
    DVPDir = "DeepVirPredictions"
    filename = glob.glob(DVPDir + "/contigs*")
    if PredTag != "skip":
        print("Running DeepVirFinder")
        
        if not os.path.exists(DVPDir): 
            subprocess.run(["mkdir","../" + DVPDir],cwd = assemblydirectory)
        
        subprocess.run(["conda","run","-n","VIRFINDER","python", pathtoDeepVirFinder + "/dvf.py", "-i", inFile,"-o","../" + DVPDir,"-c", str(threads)],cwd = assemblydirectory)
    else:
        print("DeepVirFinder was skipped")
    resultpath = filename[0]
    return resultpath


# TODO Map reads to contigs to reveal coverage of potential phages

'''def coverageFinder(read1,read2,directory,contigfilepath): #Finds coverage of the largest / first contig in a contigs.fasta file
    print("Finding coverage of assemblies")
    contigs = contigfilepath
    coveragestats = "coveragestats.txt"
    subprocess.run(["conda","run","-n","QC","bbmap.sh","ref=" + contigs,"in=" + read1,"in2=" + read2,"out=coverage_mapping.sam","nodisk=t","fast=t","covstats="+coveragestats],cwd = directory)
    
    linecount = 0
    with open(coveragestats,'r') as covfile:
        for line in covfile:
            if linecount == 1:
                coverage = float(line.split()[1])
                break
            linecount += 1

    with open(coveragestats,'r') as covfile:
        maxcoverage = 0
        for line in covfile:
            if linecount != 0:
                if line.split()[]
                elif float(line.split()[1]) > maxcoverage:
                    maxcoverage = float(line.split()[1])
                
            linecount += 1      
        coverage = maxcoverage
    return coverage, coveragestats
'''


   #Check average coverage, adjust samplerate to improve coverage
   #Check coverage before multi assembly and after


def coverageFinder(read1,read2,directory,contigfilepath): #Finds coverage of the largest / first contig in a contigs.fasta file
    print("Finding coverage of assemblies")
    contigs = contigfilepath
    coveragestats = "coveragestats.txt"
    subprocess.run(["conda","run","-n","QC","bbmap.sh","ref=" + contigs,"in=" + read1,"in2=" + read2,"out=coverage_mapping.sam","nodisk=t","fast=t","covstats="+coveragestats],cwd = directory)
    
    linecount = 0
    with open(coveragestats,'r') as covfile:
        for line in covfile:
            if linecount == 1:
                coverage = float(line.split()[1])
                break
            linecount += 1
    return coverage, coveragestats

def coverageFinderMax(read1,read2,directory,contigfilepath):
    print("Finding maximum coverage among assemblies which are larger than half of the size of the largest contig.")
    contigs = contigfilepath
    coveragestats = "coveragestats.txt"
    subprocess.run(["conda","run","-n","QC","bbmap.sh","ref=" + contigs,"in=" + read1,"in2=" + read2,"out=coverage_mapping.sam","nodisk=t","fast=t","covstats="+coveragestats],cwd = directory)
    
    linecount = 0
    maxcoverage = 0
    longestcontiglength = None
    with open(coveragestats,'r') as covfile:
        for line in covfile:
            linesplit = line.split()
            if linecount == 1:
                print(line)
                longestcontiglength = float(linesplit[2])
                maxcoverage = float(linesplit[1])
            elif linecount > 1:
                if float(linesplit[2]) > longestcontiglength / 2 and float(linesplit[1]) > maxcoverage:
                    print(line)
                    maxcoverage = float(linesplit[1])
                
            linecount += 1      
        
    return maxcoverage, coveragestats


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


def coverageFinderMaxWeighted(read1,read2,directory,contigfilepath):
    print("Finding maximum coverage among assemblies which are larger than half of the size of the largest contig.")
    contigs = contigfilepath
    coveragestats = "coveragestats.txt"
    subprocess.run(["conda","run","-n","QC","bbmap.sh","ref=" + contigs,"in=" + read1,"in2=" + read2,"out=coverage_mapping.sam","nodisk=t","fast=t","covstats="+coveragestats],cwd = directory)
    
    linecount = 0
    maxcoverage = 0
    coverageLengthProduct = 0
    with open(coveragestats,'r') as covfile:
        for line in covfile:
            linesplit = line.split()
            coverage = linesplit[1]
            length = linesplit[2]
            if linecount > 0:
                product = coverage * length
                if product > coverageLengthProduct:
                    print(line)
                    maxcoverage = coverage
                    coverageLengthProduct = product
                
            linecount += 1      
        
    return maxcoverage, coveragestats


def PHAROKKA(directory, viralcontigs,threads):

    print("Running pharokka.py")
    print("Using:", threads, "threads.")
    pathToDB = "../PHAROKKADB"
    subprocess.run(["conda", "run", "-n", "PHAROKKA", "pharokka.py","-i", viralcontigs, "-o", "pharokka","-f","-t",str(threads),"-d",pathToDB, "-g","prodigal","--meta"],cwd = directory)

    print("Pharokka.py finished running.")