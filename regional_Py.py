import os
import gzip
import argparse
import numpy
import pysam

#Opening Tabix 
def open_tabix (chr_Num):
        #Gets a list of chr numbers
        n = 22
        chr_Num = range(1,n+1)

#Opens tabix files
        for item in chr_Num:
                chr_Pos = 'chr' + str(item)
                myTabix = pysam.TabixFile('/net/akey/vol1/scratch/1KGenomes_VCF/ALL.%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.bgz' %(chr_Pos))
                header = myTabix.header
                tabix_List.append(myTabix)
        return tabix_List

def VCF (header):
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

def tabix_regions (regions):
	for line in regions:
  	     	#gets chr, start, and stop
		splt_Line = line.split('\t')
		chr_Pos = splt_Line[0]
		start_Pos = int(splt_Line[1])
        	stop_Pos = int(splt_Line[2])

		#Gets number of chr
		wordList = list(chr_Pos)
        	numList = wordList[3:]
        	chrList = ''.join(numList)

        	#Gets regions
        	myRegion = tabix_List[int(chrList)-1].fetch(chrList, start_Pos, stop_Pos)
	return myRegion

def common_list (myRegion, inds, selectInds):
	for line in myRegion:
		for i in range(len(inds)):
			for j in range(len(selectInds)):
                	        if inds[i] == selectInds[j]:
                 	               comList.append(i)
	return comList

#Skips lines that aren't needed and header lines
def numpy_array (line, comList):
	if line.startswith('##') or line.startswith('#'): 
		continue

        #Cast numList as a numpy.array
        mainLines = line.split('\t')
        ref_List = mainLines[3]
        alt_List = mainLines[4]
        numList = mainLines[9:]
        num_array = numpy.array(numList)
        comp_array = numpy.array(comList)

        if len(ref_List) > 1 and len(alt_List) > 1:
                continue

        #Subset with a list of the indices that correspond to the selectInds    
        pi_indList = num_array[comp_array]

        #Gets ref and alt columnns
        ref_List = mainLines[3]
        alt_List = mainLines[4]

        #Skips lines that have more than one base pair
        if len(ref_List) > 1 and len(alt_List) > 1:
                continue
	return pi_indList

def pi_variables(pi_indList):
        #Split line by ':'
        pi_calcList = []
        zeroList = []
        oneList = []

        for item in pi_indList:
                pi_numList = item.split(':')
                piNum = pi_numList[0]

                #Split line by '|' 
                piNumSplit = piNum.split('|') #gives string of numbers

                #Changes values to integers and adds them to a list
                x = int(piNumSplit[0])
                y = int(piNumSplit[1])
                pi_calcList.append(x)
                pi_calcList.append(y)

	return pi_calcList
      
def pi_calculation (pi_calcList):	
	#Find number of 0 and 1 in list 
        for item in pi_calcList:
                if item == 0:
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
        n = len(pi_calcList)
 
        #Calculate pi
        pi = (n/(n-1))*(2*p*(1-p))
	return pi   
 
parser = argparse.ArgumentParser(description = 'Calculate pi for specific individuals')
parser.add_argument('-i',type=str,default='',help='List of individuals')
parser.add_argument('-r',type=str,default='',help='Bedfile Regions')
args = parser.parse_args()

inds = []
VCF_Inds = []
all_Inds = []
selectInds = []
tabix_List = []
comList = []
indexList = []
pi_calcList = []
zeroList = []
oneList = []
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
tabix_List = open_tabix(chr_Num)

# Puts all individuals from VCF and tabix files into lists
inds = VCF(header)

#Finds position of common individuals
all_Inds = common_inds(inds, selectInds)
	
#Gets tabix regions
myRegion = tabix_regions(regions)

#Find common individuals 	
comList = common_list(myRegion, inds, selectInds)

#Skips lines that aren't needed and header lines
pi_indList = numpy_array(line, comList)

#Gets pi variables
pi_calcList = pi_variables(pi_indList)

#Calculate pi (n)
pi = pi_calculation(pi_calcList)

print pi
raw_input()

                #Calculate pi (n^2)     
#               for i in range(len(pi_calcList)-1):
#                       for j in range(i+1,len(pi_calcList)):
#                               if pi_calcList[i] != pi_calcList[j]:
#                                       varTotal += 1
                #Calculate mean pi      
#               length = len(pi_calcList)       
#               calc = float((length*(length-1))/2)
#               meanPi = varTotal/calc
