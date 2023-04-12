#!/usr/bin/python3
import subprocess, sys, re

read1 = sys.argv[1]
read2 = sys.argv[2]
contig = sys.argv[3]


def coverageFinderAverage(read1,read2,contigfilepath):
    contigs = contigfilepath
    coveragestats = contigfilepath + "_coveragestats.txt"
    
    subprocess.run(["bbmap.sh","ref=" + contigs,"in=" + read1,"in2=" + read2,"out=coverage_mapping.sam","nodisk=t","fast=t","covstats="+coveragestats])
    
    linecount = 0
    sumcoverage = 0
    longestcontiglength = None
    with open(coveragestats,'r') as covfile:
        for line in covfile:
            linesplit = line.split()
            if linecount == 1:
                
                longestcontiglength = float(linesplit[2])
                sumcoverage = float(linesplit[1])
            elif linecount > 1:
                if float(linesplit[2]) > longestcontiglength * 0.7:
                    
                    sumcoverage += float(linesplit[1])
                else:
                    break
                
            linecount += 1      
    averagecoverage = sumcoverage / (linecount - 1)
    subprocess.run(["rm",coveragestats])

    if averagecoverage > 20 and 100 > averagecoverage:
        reg_obj = re.search(r'subs#cov([0-1]\.[0-9]+)_read1\.fastq',read1)
        return int(float(reg_obj.group(1).strip())*100)
    else: return 0

print(coverageFinderAverage(read1,read2,contig),end = "")