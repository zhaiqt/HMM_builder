import argparse
import json

#parser = argparse.ArgumentParser( prog='calculate PMW, and find the most match string',description="", epilog='python get_PMW_table.py -i inputfile')
#parser.add_argument ('-i','--input',help='Input File Name', default="./data/mFR1.txt")

#---------------------------calculate score of given sequence------------------#
def find_high_PMW_score (PMW_dict,fragment):
	input_seq=fragment
	motif_size=len(PMW_dict)
	find_flag=False
	while not find_flag:
		stop= min(len(input_seq),motif_size)+1
		frag_window=input_seq[0:stop]
		score=0

		# print "input_seq is " + input_seq
		for i in range(motif_size):
			score= score+ PMW_dict[str(i)][input_seq[i]]
			#tmp_table=PMW_dict[str(i)]['E']
		#print frag_window +' : ' + str(score)

		if score <3: 
			input_seq=input_seq[1:]
		else:
			find_flag=True
	# print frag_window
	# print fragment.find(frag_window)
	return frag_window
#--------------------------------find all the boundaries of given sequence-----------#
def annotate_Ab_CDR (Ab_seq,species):
	PMW_FR1head_name='./data/'+species+'/FR1_head_PMW_'+species+'.json'
	PMW_FR1tail_name='./data/'+species+'/FR1_tail_PMW_'+species+'.json'
	PMW_FR2head_name='./data/'+species+'/FR2_head_PMW_'+species+'.json'
	PMW_FR2tail_name='./data/'+species+'/FR2_tail_PMW_'+species+'.json'
	PMW_FR3head_name='./data/'+species+'/FR3_head_PMW_'+species+'.json'
	PMW_FR3tail_name='./data/'+species+'/FR3_tail_PMW_'+species+'.json'
	PMW_FR4head_name='./data/'+species+'/FR4_head_PMW_'+species+'.json'

	
	with open(PMW_FR1head_name,'r') as fp:
		PMW_FR1head=json.load(fp)
	with open(PMW_FR1tail_name,'r') as fp:
		PMW_FR1tail=json.load(fp)
	with open(PMW_FR2head_name,'r') as fp:
		PMW_FR2head=json.load(fp)
	with open(PMW_FR2tail_name,'r') as fp:
		PMW_FR2tail=json.load(fp)
	with open(PMW_FR3head_name,'r') as fp:
		PMW_FR3head=json.load(fp)
	with open(PMW_FR3tail_name,'r') as fp:
		PMW_FR3tail=json.load(fp)
	with open(PMW_FR4head_name,'r') as fp:
		PMW_FR4head=json.load(fp)

	# test_string= "QVQLQQSGAELARPGASVKMSCKASGYTFTSYTMHWVKQRPGQGLEWIGYINPSSGYTKYNQKFKDKATLTADKSSSTAYMQLSSLTSEDSAVYYCAR"
	#[FR1head_pos,FR1tail_pos,FR2head_pos,FR2tail_pos,FR3head_pos,FR3tail_pos,FR4head_pos]
	try:
		FR1head=find_high_PMW_score(PMW_FR1head,Ab_seq)
		FR1head_pos=Ab_seq.find(FR1head)
		print "FR1head sequence starts at " + str(FR1head_pos) +','+ FR1head
	except:
		print "FR1head can't be found in %s." % Ab_seq

	try:
		FR1tail=find_high_PMW_score(PMW_FR1tail,Ab_seq[FR1head_pos:])
		FR1tail_pos=Ab_seq.find(FR1tail)+len(FR1tail)-1
		print "FR1head sequence ends at " + str(FR1tail_pos) +','+ FR1tail
	except:
		print "FR1tail can't be found in %s." % Ab_seq

	try:
		FR2head=find_high_PMW_score(PMW_FR2head,Ab_seq[FR1tail_pos:])
		FR2head_pos=Ab_seq.find(FR2head)
		print "FR2head sequence starts at " + str(FR2head_pos) +','+ FR2head
	except:
		print "FR2head can't be found in %s." % Ab_seq

	try:
		FR2tail=find_high_PMW_score(PMW_FR2tail,Ab_seq[FR2head_pos:])
		FR2tail_pos=Ab_seq.find(FR2tail)+len(FR2tail)-1
		print "FR2tail sequence ends at " + str(FR2tail_pos) +','+ FR2tail
	except:
		print "FR2tail can't be found in %s." % Ab_seq

	try:
		FR3head=find_high_PMW_score(PMW_FR3head,Ab_seq[FR2tail_pos:])
		FR3head_pos=Ab_seq.find(FR3head)
		print "FR3head sequence starts at " + str(FR3head_pos) +','+ FR3head
	except:
		print "FR3head can't be found in %s." % Ab_seq

	try:
		FR3tail=find_high_PMW_score(PMW_FR3tail,Ab_seq[FR3head_pos:])
		FR3tail_pos=Ab_seq.find(FR3tail)+len(FR3tail)-1
		print "FR3tail sequence ends at " + str(FR3tail_pos) +','+ FR3tail
	except:
		print "FR3tail can't be found in %s." % Ab_seq

	try:
		FR4head=find_high_PMW_score(PMW_FR4head,Ab_seq[FR3tail_pos:])
		FR4head_pos=Ab_seq.find(FR4head)
		print "FR4head sequence starts at " + str(FR4head_pos) +','+ FR4head
	except:
		print "FR4head can't be found in %s." % Ab_seq

	FRM1=Ab_seq[FR1head_pos:FR1tail_pos+1]
	CDR1=Ab_seq[FR1tail_pos+1:FR2head_pos]
	FRM2=Ab_seq[FR2head_pos:FR2tail_pos+1]
	CDR2=Ab_seq[FR2tail_pos+1:FR3head_pos]
	FRM3=Ab_seq[FR3head_pos:FR3tail_pos+1]
	CDR3=Ab_seq[FR3tail_pos+1:FR4head_pos]
	# print FRM1 
	# print CDR1
	# print FRM2
	# print CDR2
	# print FRM3
	# print CDR3
	# boundary_list=[FR1head_pos,FR1tail_pos,FR2head_pos,FR2tail_pos,FR3head_pos,FR3tail_pos,FR4head_pos]
	frag_list=[FRM1,CDR1,FRM2,CDR2,FRM3,CDR3]
	print frag_list
	#print boundary_list
	return frag_list
#	for key in PMW_FR4head:
		# print "position " + str(key)
		# print PMW_FR4head[key]
##########  main ###############

# generate a dictionary contain the HMM model

# parser = argparse.ArgumentParser( prog='annotate_Ab_CDR',description="annote_Ab_CDR, return sequence of each CDR and Framwork", epilog='python cal_PMW_score.py -i inputfile -s species')
# parser.add_argument ('-s','--species',help='species', default="mouse")
#test_string= "QVQLQQSGAELARPGASVKMSCKASGYTFTSYTMHWVKQRPGQGLEWIGYINPSSGYTKYNQKFKDKATLTADKSSSTAYMQLSSLTSEDSAVYYCAR"
species="mouse"
test_string="EVKLEESGAELVRPGASVTLSCAASGYTFTDFEIHWVKQPPVGGLEWIGTLDPETGGTAYNQNFKGRATLTADKSSSTAYMELRSLTSEDSAVYYCTRWGKKFYYYGTSYAMDYWGQGTSVTVSSAKTTPPSVYPLAPGXXATNSMVTLGCLVKGYFPEPVTVTWNSGSLSGGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKIVPRD"
annotate_Ab_CDR (test_string,species)

