#!/usr/bin/python3
import subprocess, sys

read1 = sys.argv[1]
read2 = sys.argv[2]
samplerate = sys.argv[3]
sampleseed = sys.argv[4]
covorn50 = sys.argv[5]

def SubSampling(read1,read2,sampleRate,sampleSeed,covorn50): #Subsampling using Reformat.sh 
    if covorn50 == "coverage":
        sampleRate = float(sampleRate)
        #i = i / 100
        out1 = "subs#cov"+ str(sampleRate) + "_read1.fastq"
        out2 = "subs#cov"+ str(sampleRate) + "_read2.fastq"
        subprocess.run(["reformat.sh","in=" + read1, "in2=" + read2, "out=" + out1, "out2=" + out2,"samplerate=" + str(sampleRate),"sampleseed=" + str(sampleSeed),"overwrite=true"])
        #print((int(i * 100)),end="")
    elif covorn50 == "n50":
        sampleRate = int(sampleRate) / 100
        sampleSeed = (int(sampleSeed))
        out1 = "subs#n50_"+ str(sampleSeed) + "_read1.fastq"
        out2 = "subs#n50_"+ str(sampleSeed) + "_read2.fastq"
        subprocess.run(["reformat.sh","in=" + read1, "in2=" + read2, "out=" + out1, "out2=" + out2,"samplerate=" + str(sampleRate),"sampleseed=" + str(sampleSeed),"overwrite=true"])


SubSampling(read1,read2,samplerate,sampleseed,covorn50)