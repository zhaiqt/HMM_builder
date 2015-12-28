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
		print frag_window +' : ' + str(score)

		if score <3: 
			input_seq=input_seq[1:]
		else:
			find_flag=True
	# print frag_window
	# print fragment.find(frag_window)
	return frag_window
#--------------------------------find all the boundaries of given sequence-----------#
def annotate_Ab_CDR (Ab_seq,species):
	PMW_mFR1head_name='./data/'+species+'/FR1_head_PMW_'+species+'.json'
	PMW_mFR1tail_name='./data/'+species+'/FR1_tail_PMW_'+species+'.json'
	PMW_mFR2head_name='./data/'+species+'/FR2_head_PMW_'+species+'.json'
	PMW_mFR2tail_name='./data/'+species+'/FR2_tail_PMW_'+species+'.json'
	PMW_mFR3head_name='./data/'+species+'/FR3_head_PMW_'+species+'.json'
	PMW_mFR3tail_name='./data/'+species+'/FR3_tail_PMW_'+species+'.json'
	PMW_mFR4head_name='./data/'+species+'/FR4_head_PMW_'+species+'.json'

	
	with open(PMW_mFR1head_name,'r') as fp:
		PMW_mFR1head=json.load(fp)
	with open(PMW_mFR1tail_name,'r') as fp:
		PMW_mFR1tail=json.load(fp)
	with open(PMW_mFR2head_name,'r') as fp:
		PMW_mFR2head=json.load(fp)
	with open(PMW_mFR2tail_name,'r') as fp:
		PMW_mFR2tail=json.load(fp)
	with open(PMW_mFR3head_name,'r') as fp:
		PMW_mFR3head=json.load(fp)
	with open(PMW_mFR3tail_name,'r') as fp:
		PMW_mFR3tail=json.load(fp)
	with open(PMW_mFR4head_name,'r') as fp:
		PMW_mFR4head=json.load(fp)

	# test_string= "QVQLQQSGAELARPGASVKMSCKASGYTFTSYTMHWVKQRPGQGLEWIGYINPSSGYTKYNQKFKDKATLTADKSSSTAYMQLSSLTSEDSAVYYCAR"
	try:
		mFR1head=find_high_PMW_score(PMW_mFR1head,Ab_seq)
		mFR1head_pos=Ab_seq.find(mFR1head)
		print "mFR1head sequence starts at " + str(mFR1head_pos) +','+ mFR1head
	except:
		print "mFR1head can't be found in %s." % Ab_seq

	try:
		mFR1tail=find_high_PMW_score(PMW_mFR1tail,Ab_seq[mFR1head_pos:])
		mFR1tail_pos=Ab_seq.find(mFR1tail)+len(mFR1tail)-1
		print "mFR1head sequence ends at " + str(mFR1tail_pos) +','+ mFR1tail
	except:
		print "mFR1tail can't be found in %s." % Ab_seq

	try:
		mFR2head=find_high_PMW_score(PMW_mFR2head,Ab_seq[mFR1tail_pos:])
		mFR2head_pos=Ab_seq.find(mFR2head)
		print "mFR2head sequence starts at " + str(mFR2head_pos) +','+ mFR2head
	except:
		print "mFR2head can't be found in %s." % Ab_seq

	try:
		mFR2tail=find_high_PMW_score(PMW_mFR2tail,Ab_seq[mFR2head_pos:])
		mFR2tail_pos=Ab_seq.find(mFR2tail)+len(mFR2tail)-1
		print "mFR1head sequence ends at " + str(mFR2tail_pos) +','+ mFR2tail
	except:
		print "mFR2tail can't be found in %s." % Ab_seq

	try:
		mFR3head=find_high_PMW_score(PMW_mFR3head,Ab_seq[mFR2tail_pos:])
		mFR3head_pos=Ab_seq.find(mFR3head)
		print "mFR3head sequence starts at " + str(mFR3head_pos) +','+ mFR2head
	except:
		print "mFR3head can't be found in %s." % Ab_seq

	try:
		mFR3tail=find_high_PMW_score(PMW_mFR3tail,Ab_seq[mFR3head_pos:])
		mFR3tail_pos=Ab_seq.find(mFR3tail)+len(mFR3tail)-1
		print "mFR3head sequence ends at " + str(mFR3tail_pos) +','+ mFR3tail
	except:
		print "mFR3tail can't be found in %s." % Ab_seq

	try:
		mFR4head=find_high_PMW_score(PMW_mFR4head,Ab_seq[mFR3tail_pos:])
		mFR4head_pos=Ab_seq.find(mFR4head)
		print "mFR4head sequence starts at " + str(mFR4head_pos) +','+ mFR4head
	except:
		print "mFR4head can't be found in %s." % Ab_seq

	return
#	for key in PMW_mFR4head:
		# print "position " + str(key)
		# print PMW_mFR4head[key]
##########  main ###############

# generate a dictionary contain the HMM model

# parser = argparse.ArgumentParser( prog='annotate_Ab_CDR',description="annote_Ab_CDR, return sequence of each CDR and Framwork", epilog='python cal_PMW_score.py -i inputfile -s species')
# parser.add_argument ('-s','--species',help='species', default="mouse")
#test_string= "QVQLQQSGAELARPGASVKMSCKASGYTFTSYTMHWVKQRPGQGLEWIGYINPSSGYTKYNQKFKDKATLTADKSSSTAYMQLSSLTSEDSAVYYCAR"
species="mouse"
test_string="EVKLEESGAELVRPGASVTLSCAASGYTFTDFEIHWVKQPPVGGLEWIGTLDPETGGTAYNQNFKGRATLTADKSSSTAYMELRSLTSEDSAVYYCTRWGKKFYYYGTSYAMDYWGQGTSVTVSSAKTTPPSVYPLAPGXXATNSMVTLGCLVKGYFPEPVTVTWNSGSLSGGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKIVPRD"
annotate_Ab_CDR (test_string,species)


