import os
import gzip
import argparse
import numpy
import pysam
import sys

n = 22 
chr_Num = range(1,n+1)

#Opening Tabix 
def open_tabix (chr_Num):
        
	for item in chr_Num:
                chr_Pos = 'chr' + str(item)
                myTabix = pysam.TabixFile('/net/akey/vol1/scratch/1KGenomes_VCF/ALL.%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.bgz' %(chr_Pos))
                header = myTabix.header
                tabix_List.append(myTabix)

	        for line in header:
        	        #Skips line that aren't needed
			if line.startswith('##'):
				continue
	                splt_Header = line.split('\t')
        	        inds.append(splt_Header[9:])
        return inds

#Function for finding 
def common_inds (inds, selectInds):
	for i in range(len(inds)):
		temp_List = []
		for j in range(len(inds[i])):
			for k in range(len(selectInds)):
				if inds[i][j] == selectInds[k]:
					temp_List.append(j)
		all_Inds.append(temp_List)
	return all_Inds

def tabix_regions (lines, pad):
      	#gets chr, start, and stop
	splt_Line = lines.split('\t')
	chr_Pos = splt_Line[0]
	start_Pos = int(splt_Line[1])
        stop_Pos = int(splt_Line[2])
	ind_Pos = splt_Line[3]
	len_Pos = splt_Line[5]
	num_Inds = splt_Line[6]	
	
	beg_Pad = start_Pos - pad
	end_Pad = stop_Pos + pad

	if (chr_Pos == 'X') or (chr_Pos == 'Y'):
		return

	#Gets number of chr
	wordList = list(chr_Pos)
        numList = wordList[3:]
        chrList = ''.join(numList)

        #Gets regions
        myRegion = tabix_List[int(chrList)-1].fetch(chrList, beg_Pad, end_Pad)
	#padRegion = tabix_List[int(chrList)-1].fetch(chrList, beg_Pad, end_Pad)
	return chr_Pos, start_Pos, stop_Pos, ind_Pos, len_Pos, num_Inds, myRegion

#Skips lines that aren't needed and header lines
def new_array(myRegion, all_Inds):
	pi_indList = []
	for line in myRegion:
		if line.startswith('##') or line.startswith('#'): 
			continue

        #Cast numList as a numpy.array
		mainLines = line.split('\t')
		ref_List = mainLines[3]
		alt_List = mainLines[4]
		numList = mainLines[9:]
        	num_array = numpy.array(numList)
        	comp_array = numpy.array(all_Inds[0]) 
		if len(ref_List) > 1 and len(alt_List) > 1:
                	continue #probably these returns should just be continues

        #Subset with a list of the indices that correspond to the selectInds 
		pi_indList.append(num_array[comp_array])
			 
	
	return pi_indList

def pi_variables(pi_indList):
        #Split line by ':'
        pi_calcList = []

        for item in pi_indList:
		#for each site create a new list of 0s and 1s
		temp_calcList = []
		for h in item:
	       		pi_numList = h.split(':')
                	piNum = pi_numList[0]

                #Split line by '|' 
                	piNumSplit = piNum.split('|') #gives string of numbers

                #Changes values to integers and adds them to a list
                	x = int(piNumSplit[0])
                	y = int(piNumSplit[1])
                	temp_calcList.append(x)
                	temp_calcList.append(y)
		#append that new list to pi_calcList
		pi_calcList.append(temp_calcList)		

	return pi_calcList
      
def pi_calculation (pi_calcList):	

	zeroList = []
	oneList = []

	totalZero = []
	totalOne = []

	piTotal = 0.0

	#Find number of 0 and 1 in list (add another loop) 
        for item in pi_calcList:
		zeroList = []
		oneList = []
		for allele in item:
                	if allele == 0:
                        	zero = item
                        	zeroList.append(zero)
                	else:
                        	one = item
                        	oneList.append(one)

        	#Get variables to make equations        
       		numZero =  float(len(zeroList))
        	numOne = float(len(oneList))

	        #Make more variables for equations
        	p = numZero/(numZero + numOne)
        	n = float(len(item))
 
        	#Calculate pi
        	pi = (n/(n-1))*(2.0*p*(1-p))
		piTotal += pi
	return piTotal
	#needs to be divided by the total number of bases in the region 
					#I would return piTotal
					#back in the main part of the script, divide piTotal by the total number of sites
					#i.e. stop-start+1 
 
parser = argparse.ArgumentParser(description = 'Calculate pi for specific individuals')
parser.add_argument('-i',type=str,default='',help='List of individuals')
parser.add_argument('-r',type=str,default='',help='Bedfile Regions')
parser.add_argument('-p',type=int,default=0,help='Add Pad')
args = parser.parse_args()

inds = []
VCF_Inds = []
all_Inds = []
selectInds = []
tabix_List = []
comList = []
indexList = []
varTotal = 0

#open the list of individuals to use to calculate pi
openIndividuals = open(args.i,"r") 

#go through that list, storing the names of the inds in a list
for line in openIndividuals:
        #store each individual in a list
        selectInds.append(line[:-1])

regions = open(args.r,"r")
	
#skip the first line
regions.readline()

#Opens tabix files
inds = open_tabix(chr_Num)

#Finds position of common individuals
all_Inds = common_inds(inds, selectInds)
	
#Gets tabix regions

pad = args.p

for line in regions:
	lines = line.strip()

	chr_Pos, start_Pos, stop_Pos,inds_Pos, len_Pos, num_Inds, myRegion = tabix_regions(lines, pad)

#	beg_Pad = start_Pos - pad
#	end_Pad = stop_Pos + pad

#	print padRegion
#	raw_input()

	if myRegion is None:
		sys.stderr.write("myRegion is None\n")
		continue

	#Skips lines that aren't needed and header lines
	pi_indList = new_array(myRegion, all_Inds) #should process lines from a VCF

	if pi_indList is None:
		sys.stderr.write("pi_indList is None\n")
		continue


	#Gets pi variables
	pi_calcList = pi_variables(pi_indList)

	#Calculate pi (n)
	piTotal = pi_calculation(pi_calcList)

	final_Pi = piTotal/((stop_Pos - start_Pos + 1)*2)
	print chr_Pos, start_Pos, stop_Pos, inds_Pos, len_Pos, num_Inds, final_Pi
#	print "%s\t%f"%(line.strip(), pi)
#	raw_input()
