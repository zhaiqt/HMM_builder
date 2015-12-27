#import sys
import argparse
import math
import json
# generate a dictionary contain the HMM model

parser = argparse.ArgumentParser( prog='build HMM model',description="using fastafiles from IMGT", epilog='python HMMbuilder.py -i inputfile -o outputfile')
parser.add_argument ('-i','--input',help='Input File Name', default="./data/mFR1.txt")
parser.add_argument('-n','--size', help="string size", default=6)
parser.add_argument('-t','--tail',help='3 prime end',action='store_false')
parser.print_help()
args=parser.parse_args()

##show values##
#print (parser.print_help())
print ("input file: %s" % args.input)
print ("string size: %s" % args.size)
print ("Is it 5 prime end? %s" % args.tail)

#------parse the input file, generate key value pair of each sequences {Key: (sequence, length of sequence)}
def read_ProFastaFiles(fastaFile):
	Infile1 = open(fastaFile, 'r')
	fastaName=''
	fastaSeq=''
	allDict={}
	aminoAcids=set('ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv')
	
	#record the length of each sequence
	#lenSeq=()

	for line in Infile1:
		line = line.strip().strip('\n')
		if line and line.startswith(">"):
			line=line.lstrip('>')
			#process previous entry
			allDict[fastaName]= fastaSeq.upper()
			#begin a new entry
			fastaName=line
			fastaSeq=''
		elif line =='':  #/check whether it is the last line
			allDict[fastaName]= fastaSeq
			#lenSeq.append(len(fastaSeq))
		elif  set(line) <= aminoAcids and len(line)>10: 
			fastaSeq = fastaSeq + line
	Infile1.close()
	return allDict;

##---------QC fasta files, remove the entry which is too short to too long----------#
def remove_outlier(fastaDic, length_ruler):
	#record the length of each sequence
	cleaned_dict={}
	for seqName in fastaDic:
		if len(fastaDic[seqName]) ==length_ruler:
			cleaned_dict[seqName]=fastaDic[seqName]
	return cleaned_dict

##--------- define median -----------------------------------------#
# the input the fasta dictionary
def median(fastaDic):
	lenSeq=[]
	for seqName in fastaDic:
		lenSeq.append(len(fastaDic[seqName]))

	#check whether only numbers in the list
	for i in lenSeq:
		if not isinstance(i,(int,long,float,complex)):
			return "Error: the list is not composed of pure number"
	#find the median
	lenSeq=sorted(lenSeq)
	media_pos=int(len(lenSeq)/2)
	print "The median length of the fasta is " + str(lenSeq[media_pos-1])+ '\n'
	return lenSeq[media_pos-1]

##-----------generate the AA count matrix---------------------------------------#
def get_PWM_table(trunFastaDic,string_size, is5prime):

	for entryName in trunFastaDic:
		seq=trunFastaDic[entryName]
		if string_size> len(seq):
			print "Error:the request matrix string size is bigger than the fasta "
			return 

	# count the amino acid in each position
	count_matrix={}
	delta=0.000001 
	if is5prime:	
		for i in range(string_size):
			pos_matrix={}
			aminoAcids='ARNDCEQGHILKMFPSTWYV'
			for aminoAcid in aminoAcids: # initiate the 0 count, without it, error occur
				pos_matrix[aminoAcid]=0
			for entryName in trunFastaDic:# count absolute times
				pos_matrix[trunFastaDic[entryName][i]] += 1
			for aminoAcid in pos_matrix:
				pos_matrix[aminoAcid]=math.log( float(pos_matrix[aminoAcid])/len(trunFastaDic)/0.05 +delta)
			count_matrix[i]=pos_matrix
	else:
		for i in range(len(seq)-string_size,len(seq)):
			pos_matrix={}
			aminoAcids='ARNDCEQGHILKMFPSTWYV'	
			for aminoAcid in aminoAcids: # initiate the 0 count, without it, error occur
				pos_matrix[aminoAcid]=0
			for entryName in trunFastaDic:# count absolute times
				pos_matrix[trunFastaDic[entryName][i]] += 1
			for aminoAcid in pos_matrix:
				#pos_matrix[aminoAcid]=float(pos_matrix[aminoAcid])/len(trunFastaDic)/0.05
				pos_matrix[aminoAcid]=math.log( float(pos_matrix[aminoAcid])/len(trunFastaDic)/0.05 +delta)
			count_matrix[i]=pos_matrix
			#count_matrix[i+string_size-len(seq)]=pos_matrix	
	#print the resulting dictionary for checkup
	# for key in count_matrix:
	# 	print "position " + str(key)
	# 	print count_matrix[key]
	#print count_matrix

	return count_matrix




############################## main ##############################################
allFasta = read_ProFastaFiles(args.input)
print str(len(allFasta))+' sequences have been read.'

#remove key value pair , whose sequence lenth is too short
medianLenth=median(allFasta)
cleanedFasta=remove_outlier(allFasta,medianLenth)
print str(len(cleanedFasta))+' sequences are used for generate PMW table.'

# genreate PWMtable, write it to Json format
PWMtable=get_PWM_table(cleanedFasta,args.size, args.tail)
if args.tail :
	PWMtable_name=args.input.rstrip(".txt")+"_"+str(args.size)+"_head_PMW_table.json"
else:
	PWMtable_name=args.input.rstrip(".txt")+"_"+str(args.size)+"_tail_PMW_table.json"
with open(PWMtable_name,'w') as fp:
	json.dump(PWMtable,fp)
print "json file has been writen to the file:" + PWMtable_name

# get PMWtable dictionary from Json format
with open(PWMtable_name,'r') as fp:
	read_PMW_dict=json.load(fp)
for key in read_PMW_dict:
	print "position " + str(key)
	print read_PMW_dict[key]
#print read_PMW_dict
