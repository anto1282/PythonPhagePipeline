#!/usr/bin/python3

import subprocess, sys, os, multiprocessing, glob
from argparse import ArgumentParser
import Assembly
import TaxRemover
import Kraken2
from DeepVirExtractor import DeepVirExtractor

#Command line arguments

parser = ArgumentParser(prog = 'PhAnTom.py',description="Help me!")
parser.add_argument("-i","--input", action="store", dest = "sraAccNr", help ="Input valid accesion number from Sequence Read Archive.")
parser.add_argument("--flag",'-f', action = "store", dest = "whatSPADES", default = "--meta", help = "Runs spades with the specified flag. (skip skips spades)")
parser.add_argument("--threads","-t",action = "store", dest = "threads", type = int, default = multiprocessing.cpu_count(),help = "Specifies the total number of threads used in various of the programs in the pipeline.")
parser.add_argument("-r","--ref", action="store", dest="refFile", help ="Input fasta file for filtering out", default=False)
parser.add_argument("-p","-pred",action="store",dest = "virpredflag", default = "dontskip", help = "Write skip to skip deepvirfinder.")
parser.add_argument("-a","--assemblies",action="store",dest = "nrofassemblies",default = "1", help ="Determines the number of assemblies and subsamplings to be performed (Default is 1)")
parser.add_argument("-s", "--skip", action="store", dest = "skip",nargs="+",default= "1", help="Please input the steps you want to skip: \n Wget \n Trimming \n Assembly \n DeepVirFinder \n Pharokka") 


args = parser.parse_args()
#print("Hello", args.skip)

#Creates results directory for the pipeline results
def directory_maker(sraAccNr):
    parent_directory =  "../Results" + sraAccNr
    if not os.path.exists(parent_directory): 
        subprocess.run(["mkdir",parent_directory])
    return parent_directory

#Downloads sequence reads from SRA using fasterq-dump
def sra_get(sraAccNr,directory):
    if not os.path.exists(directory + "/" +sraAccNr + "_1.fastq") and not os.path.exists(directory + "/" + sraAccNr + "_2.fastq"):
        print("Downloading", sraAccNr, "from SRA using fasterq-dump")
        subprocess.run(["fasterq-dump", sraAccNr, "--progress", "--split-files"], cwd = directory)
    else:
        print("Reads already present in directory! Continuing...")
    return sraAccNr + "_1.fastq", sraAccNr + "_2.fastq"


def wget_wget(sraAccNr,directory): #Replacement for fasterq-dump
    if not os.path.exists(directory + "/" +sraAccNr + "_1.fastq") and not os.path.exists(directory + "/" + sraAccNr + "_2.fastq"):
        print("Downloading", sraAccNr, "from SRA using wget")
        link = "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fasta?acc=" + sraAccNr
        filename = sraAccNr + ".gz"
        subprocess.run(["wget", link,"-O",filename], cwd = directory)
        subprocess.run(["gzip","-d",filename],cwd = directory)
    else:
        print("Reads already present in directory! Continuing...")
    return sraAccNr + "_1.fastq", sraAccNr + "_2.fastq"


#TODO Spørg Bent/Thomas/Ole om programmet skal kunne ændre i hvor meget den trimmer af (Kvalitet)
def trimming(read1, read2, directory, refFile):
    

    read1_trimmed = "trimmed/" + read1
    read2_trimmed = "trimmed/" + read2

    if not os.path.exists(directory + "/trimmed"):
        subprocess.run(["mkdir","trimmed"], cwd = directory) 

    adaptoutput1 = read1 + ".noadap"
    adaptoutput2 = read2 + ".noadap"

    subprocess.run(["mamba", "run", "-n", "QC","AdapterRemoval","--file1", read1,  "--file2", read2, "--output1", adaptoutput1, "--output2", adaptoutput2], cwd =directory)
    print("AdapterRemoval finished.")

    if args.refFile:
        subprocess.run(["mamba", "run", "-n", "QC","bbduk.sh","-in=" + adaptoutput1,  "-in2=" + adaptoutput2, "-out=" + read1_trimmed, "-out2=" + read2_trimmed, "ref=" + refFile , "trimq=25", "qtrim=w","forcetrimleft=15" ,"overwrite=true"], cwd =directory)
        print("Trim finished.")
    else:
        subprocess.run(["mamba", "run", "-n", "QC","bbduk.sh","-in=%s" % adaptoutput1,  "-in2=%s" % adaptoutput2, "-out=%s" % read1_trimmed, "-out2=%s" % read2_trimmed, "trimq=25", "qtrim=w","forcetrimleft=15" ,"overwrite=true"], cwd =directory)
        print("Trim finished.")
    
    #Cleaning up
    subprocess.run(["rm",adaptoutput1],cwd = directory)
    subprocess.run(["rm",adaptoutput2],cwd = directory)

    return read1_trimmed, read2_trimmed

def main():
    #Initialization

    threads = args.threads

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()

    #Implement steps to stop pipeline if arguments are missing

    sraAccNr = args.sraAccNr 
    
    parent_directory = directory_maker(sraAccNr)

    #read1, read2 = sra_get(sraAccNr,parent_directory)
    read1, read2 = sra_get(sraAccNr,parent_directory) #Testing wget instead of fasterq
   
   
    phredOffset = Assembly.offsetDetector(read1,read2,parent_directory)

    refFile = args.refFile
    if "Trimming" not in args.skip:
        read1Trimmed, read2Trimmed = trimming(read1, read2, parent_directory, refFile)
        Kraken2.Kraken(parent_directory,read1Trimmed,read2Trimmed, "../KrakenDB")
    
        read1Trimmed, read2Trimmed = TaxRemover.EuRemover(parent_directory,read1Trimmed, read2Trimmed, sraAccNr)
    else:
        read1Trimmed, read2Trimmed = "trimmed/" + read1, "trimmed/" + read2
        print("Trimming was skipped")

    print(read1Trimmed,read2Trimmed)
    assemblydirectory, read1Trimmed, read2Trimmed = Assembly.MultiAssembly(read1Trimmed,read2Trimmed,parent_directory,phredOffset,0.3,args.nrofassemblies, args.skip)

    
    #Maybe this is not necessary?
    Contigs_Trimmed = Assembly.contigTrimming(assemblydirectory, "contigs.fasta", minLength=500) #Filters off too short contigs

    pathToDeepVirFinder = "../../DeepVirFinder"
    if "DeepVirFinder" not in args.skip:
        predfile = Assembly.DeepVirFinder(pathToDeepVirFinder, assemblydirectory,threads, Contigs_Trimmed, args.virpredflag)
    else:
        predfile = glob.glob("DeepVirPredictions" + "/contigs*")[0]
        print("DeepVirFinder was skipped")

    viralcontigs,nonviralcontigs = DeepVirExtractor(predfile,assemblydirectory,parent_directory,0.95)
    
    if "Pharokka" not in args.skip:
        Assembly.PHAROKKA(parent_directory, viralcontigs, threads)
    else:
        print("Pharokka was skipped")
    
    Kraken2.KrakenContig(parent_directory,nonviralcontigs, "../KrakenDB")



main()

### The end ###
print("\nThanks for using the PhAnTomic's pipeline :-P")