#import sys
import argparse
# generate a dictionary contain the HMM model

parser = argparse.ArgumentParser( prog='build HMM model',description="using fastafiles from IMGT", epilog='python HMMbuilder.py -i inputfile -o outputfile')
parser.add_argument ('-i','--input',help='Input File Name', default="./data/mFR1.txt")
#parser.add_argument('-5','--head', help="trim 5' region", default='store_true')
#parser.add_argument('-o','--output',help='Output File Name',required="true")
parser.print_help()
args=parser.parse_args()

##show values##
#print (parser.print_help())
print ("input file: %s" % args.input)


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
def generate_AA_matrix(trunFastaDic,string_size, is5prime):

	for entryName in trunFastaDic:
		seq=trunFastaDic[entryName]
		if string_size> len(seq):
			print "Error:the request matrix string size is bigger than the fasta "
			return 

	count_matrix={}
	if is5prime:	
		for i in range(string_size-1):
			pos_matrix={}
			aminoAcids='ARNDCEQGHILKMFPSTWYV'
			for aminoAcid in aminoAcids:
				pos_matrix[aminoAcid]=0

			for entryName in trunFastaDic:
				pos_matrix[trunFastaDic[entryName][i]] += 1

			count_matrix[i]=pos_matrix
	else:
		for i in range(len(seq)-string_size,len(seq)):
			pos_matrix={}
			aminoAcids='ARNDCEQGHILKMFPSTWYV'
			for aminoAcid in aminoAcids:
				pos_matrix[aminoAcid]=0

			for entryName in trunFastaDic:
				pos_matrix[trunFastaDic[entryName][i]] += 1

			count_matrix[i+string_size-len(seq)]=pos_matrix			

	# print the resulting dictionary for checkup
	for key in count_matrix:
		print "position " + str(key)
		print count_matrix[key]

	return count_matrix


########## prepare input and output files ###############################
Infilename1 = args.input
Infile1 = open(Infilename1, 'r')
	
Outfilename1=args.input.rstrip(".txt")+"HMM.txt"
Outfile1 = open(Outfilename1, 'w')

Outfilename2=args.input.rstrip(".txt")+"remove_outlier.txt"
Outfile2 = open(Outfilename2, 'w')

############################## main ##############################################
allFasta={}
allFasta = read_ProFastaFiles(args.input)

for key in allFasta:
	Outfile1.write(str(key) + "\t" + allFasta[key] + '\n')
	if not allFasta[key].upper().startswith('N'):
		# generate trimmed fasta file and dictinary	
		Outfile1.write(str(key) + "\t" + allFasta[key] + '\n')
print str(len(allFasta))+' sequences have been writen into ' + Outfilename1 + '\n'

#remove key value pair , whose sequence lenth is too short

medianLenth=median(allFasta)
cleanedFasta=remove_outlier(allFasta,medianLenth)

AA_matrix=generate_AA_matrix(cleanedFasta,8, False)

for key in cleanedFasta:
	Outfile1.write(str(key) + "\t" + cleanedFasta[key] + '\n')
	if not cleanedFasta[key].upper().startswith('N'):
		# generate trimmed fasta file and dictinary	
		Outfile2.write(str(key) + "\t" + cleanedFasta[key] + '\n')
print str(len(cleanedFasta))+' sequences without short ones have been writen into ' + Outfilename2 + '\n'

Outfile1.close()
Outfile2.close()
