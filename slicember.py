#requires samtools, bwa, velvet (idba_ud/Ray/SPAdes)
import sys, random, string, time, re, os, collections, shutil, getopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append('/24-2/home/hamid/Tools/suffix_tree-2.1/build/lib')  # not good
import suffix_tree

# help functions
def ReverseComplement(seq):
	seq1 = 'ATCGTAGCatcgtagcNNNNN'
	seq_dict = { seq1[i]:seq1[i+4] for i in range(17) if i < 4 or 8<=i<12 or i==16 }
	return "".join([seq_dict[base] for base in reversed(seq)])

def TrimRpttvEnds(seq):
	chck_len = 20
	pattern, rep_patterns = re.compile(r'(.+?)\1+'), []
	rep_patterns.extend(pattern.findall(seq[0:chck_len]) or [seq[0:chck_len]])
	rep_patterns = list(set(rep_patterns))
	rep_patterns.sort(key = len)
	for pat in rep_patterns[::-1]:
		if seq.startswith(pat):
			while seq.startswith(pat):
				seq = seq[len(pat):]
			seq = pat + seq
			break
	l = len(seq)
	pattern, rep_patterns = re.compile(r'(.+?)\1+'), []
	rep_patterns.extend(pattern.findall(seq[l-chck_len:l]) or [seq[l-chck_len:l]])
	rep_patterns = list(set(rep_patterns))
	rep_patterns.sort(key = len)
	for pat in rep_patterns[::-1]:
		if seq.endswith(pat):
			while seq.endswith(pat):
				seq = seq[:-len(pat)]
			seq = seq + pat
			break
	return seq

def total_len(l):
	total_len = 0 
	for record in l:
		total_len = total_len + len(record.seq)
	return total_len

def get_longest_overlap(s1, s2):
	i = len(s2)
	while not s1.endswith(s2[0:i]):
		i -= 1
	j = len(s1)
	while not s2.endswith(s1[0:j]):
		j -= 1
	if i > j:
		return s2[0:i]
	else:
		return s1[0:j]

def merge(s1,s2, lcps):
	splt_s1 = s1.split(lcps,1)
	splt_s2 = s2.split(lcps,1)
	if len(splt_s2[0]) == 0 and len(splt_s1[0]) > 0:
		return splt_s1[0]+lcps+splt_s2[1]
	elif len(splt_s2[0]) > 0 and len(splt_s1[0]) == 0:
		return splt_s2[0]+lcps+splt_s1[1]
	elif len(splt_s2[0]) == 0 and len(splt_s1[0]) ==0:
		if lcps in splt_s1[1]:
			return lcps + merge(splt_s1[1],s2,lcps)
		elif lcps in splt_s2[1]:
			return lcps + merge(s1,splt_s2[1],lcps)


def num_bridge_reads(ns1,ns2,orphan_reads): 
	s1_ors = [line for line in orphan_reads if line[1]==ns1]
	s2_ors = [line for line in orphan_reads if line[1]==ns2]
	s1_comm = [itm for itm in s1_ors if itm[0] in (x[0] for x in s2_ors)]
	s2_comm = [itm for itm in s2_ors if itm[0] in (x[0] for x in s1_ors)]
	comm_ors = s1_comm + s2_comm
	skt = set(tuple(i) for i in comm_ors)
	return (len(comm_ors)-(2*(len(comm_ors)-len(skt))))/2

# main functions
def do_slicing(working_dir, input_reads, ref_est_len, des_cov):
	print "Step 0: Slicing reads:"

	# find number of required slices and average read length
	num_reads = 0
	total_read_len = 0
	handle = open(input_reads, "rU")
	for record in SeqIO.parse(handle, "fastq") :
		num_reads = num_reads + 1
		total_read_len = total_read_len + len(record.seq)
	ave_read_len = total_read_len / num_reads
	num_slc = ((num_reads * ave_read_len) / (ref_est_len * des_cov))

	# make the directories
	if not os.path.exists(working_dir + '/input_data'):
		os.makedirs(working_dir + '/input_data')
	if not os.path.exists(working_dir + '/it1'):
		os.makedirs(working_dir + '/it1')
	if not os.path.exists(working_dir + '/it1/conserved_regions'):
		os.makedirs(working_dir + '/it1/conserved_regions')
	for i in range(num_slc):
		if not os.path.exists(working_dir + '/it1/slice' + str(i)):
			os.makedirs(working_dir + '/it1/slice' + str(i))

	##os.system('''awk '{OFS="\/t"; getline seq1; getline sep1; getline qual1; getline header2; getline seq2; getline sep2; getline qual2; print $0,seq1,sep1,qual1,header2,seq2,sep2,qual2}' ''' + input_reads +' | shuf | awk ''{OFS="\/n"; print $1,$2,$3,$4,$5,$6,$7,$8}'' > ' + working_dir + '/input_data/reads.shuffled.fq')
	os.system('cp ' + input_reads + ' ' + working_dir + '/input_data/reads.shuffled.fq');

	os.system('split -d -l ' + str(((ref_est_len * des_cov)/(ave_read_len*2))*8) + ' ' +working_dir + '/input_data/reads.shuffled.fq ' + working_dir + '/it1/reads_slc');
	for i in range(10):
		os.system('mv ' + working_dir + '/it1/reads_slc0' + str(i) + ' ' + working_dir + '/it1/reads_slc' + str(i));
	for i in range(num_slc):
		os.system('mv ' + working_dir + '/it1/reads_slc' + str(i) + ' ' + working_dir + '/it1/slice' + str(i) + '/reads_slc' + str(i) + '.fastq');
	os.system('rm ' + working_dir + '/it1/reads_slc*');

	print "	" + str(num_slc) + " slices were created"
	return num_slc

def assemble_slices(working_dir, num_slc, curr_it, k_mer=69):
	print "it " + str(curr_it) + " :"
	print "	Step 1: assembling " + str(num_slc) + " slices..."

	for i in range(num_slc):
		
		#IDBA_ud
		#os.system('~/Tools/idba-1.1.1/bin/fq2fa --paired --filter ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fastq ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fa');
		#os.system('~/Tools/idba-1.1.1/bin/idba_ud -r ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fa --mink 25  --maxk 65 --step 20  --num_threads 10 -o ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i));
		
		# Velvet
		os.system('velveth ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + ' ' + str(k_mer) + ' -short -fastq ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fastq')
		os.system('velvetg ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + ' -cov_cutoff 10 -exp_cov auto -min_contig_lgth 200')

		#SPADes
		#os.system('~/Tools/SPAdes/SPAdes-3.1.0/bin/spades.py --careful -k 25,45,65 --12 ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fastq  -o ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i));

		#Ray
		#os.system('mpiexec -n 80 Ray -k ' + str(k_mer) + ' -i ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fastq  -o ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/ray');

		time.sleep(5)

		#idba
		#os.system("sed -e '/^>/! s/[[:lower:]]/N/g' " + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contig.fa > ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contigs_slc' + str(i) + '.fa');
		
		#velvet
		os.system("sed -e '/^>/! s/[[:lower:]]/N/g' " + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contigs.fa > ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contigs_slc' + str(i) + '.fa');

		#spades
		#os.system("sed -e '/^>/! s/[[:lower:]]/N/g' " + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contigs.fasta > ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contigs_slc' + str(i) + '.fa');

		#Ray
		#os.system("sed -e '/^>/! s/[[:lower:]]/N/g' " + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/r.Contigs.fasta > ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/contigs_slc' + str(i) + '.fa');

def find_cnsrvd_rgns(working_dir, num_slc, curr_it, min_cr_len, des_cr_sbsts):
	print "		Trying: " + str(des_cr_sbsts) + " out of " + str(num_slc) + " slices. Min length is " + str(min_cr_len)
	
	# current iteration directory
	cur_it_dir = working_dir + '/it' + str(curr_it)

	cnsrvd_rgns = []
	mrgd_cntgs_lst = []
	tmp_mrgd_cntgs_lst =[]
	tmp_rvrs_cmplmnt_mrgd_cntgs_lst = []
	fcs = ['A','T','C','G','a','t','c','g','N','$',"'",'"']

	# build a suffix tree from all the assemblies and their reverse compliments
	for curr_slc in range(num_slc):
		mrgd__cntgs = ''
		rvrs_cmplmnt_mrgd_cntgs = ''
		seperator_char = '%'

		#idba
		#handle = open(cur_it_dir + '/slice' + str(curr_slc)+ '/contig.fa', 'r')

		#velvet
		handle = open(cur_it_dir + '/slice' + str(curr_slc)+ '/contigs.fa', 'r')

		#spades
		#handle = open(cur_it_dir + '/slice' + str(curr_slc)+ '/contigs.fasta', 'r')

		#Ray
		#handle = open(cur_it_dir + '/slice' + str(curr_slc)+ '/ray.Contigs.fasta', 'r')

		seperator_char1 = 'A'
		while seperator_char1 in fcs:
			seperator_char1 = random.choice(string.printable[0:94])
		fcs.append(seperator_char1)
		seperator_char2 = 'A'
		while seperator_char2 in fcs:
			seperator_char2 = random.choice(string.printable[0:94])
		fcs.append(seperator_char2)

		for record in SeqIO.parse(handle, "fasta"):
			mrgd__cntgs = mrgd__cntgs + seperator_char1 + str(record.seq)
			rvrs_cmplmnt_mrgd_cntgs = rvrs_cmplmnt_mrgd_cntgs + seperator_char2 + ReverseComplement(str(record.seq))
		tmp_mrgd_cntgs_lst.append(mrgd__cntgs)
		tmp_rvrs_cmplmnt_mrgd_cntgs_lst.append(rvrs_cmplmnt_mrgd_cntgs)
	mrgd_cntgs_lst = tmp_mrgd_cntgs_lst + tmp_rvrs_cmplmnt_mrgd_cntgs_lst
	st = suffix_tree.GeneralisedSuffixTree(mrgd_cntgs_lst)
	##

	# find segments conserved in at least k assemblies out of n
	for shared in st.k_of_n_sharedSubstrings(minimumLength=min_cr_len,num_strngs=num_slc, des_num_subsets=des_cr_sbsts):
		for seq,start,stop in shared:
			shrd_seq = TrimRpttvEnds(mrgd_cntgs_lst[seq][start:stop])
			rc_shrd_seq = ReverseComplement(shrd_seq)

			should_remove = [cnsrvd_rgns.index(s) for s in cnsrvd_rgns if (s in shrd_seq or s in rc_shrd_seq)]
			cnsrvd_rgns = [i for j, i in enumerate(cnsrvd_rgns) if j not in should_remove]
			if not any(shrd_seq in cr or rc_shrd_seq in cr for cr in cnsrvd_rgns):
				cnsrvd_rgns.append(shrd_seq)
	##

	# write the conserved segments to a file
	counter = 0
	cr_records = []
	for item in cnsrvd_rgns:
		counter = counter + 1
		curr_id = 'it' + str(curr_it) + '_org_cr_' + str(counter)
		record = SeqRecord(Seq(item),id=curr_id)
		cr_records.append(record)

	if len(cr_records)>0:
		output_handle = open(cur_it_dir + '/conserved_regions/org_conserved_regions_it' + str(curr_it) + '.fasta', 'w')
		SeqIO.write(cr_records, output_handle, "fasta")
		output_handle.close()
		##

		time.sleep(5)
		# find orphan reads along the new conserved segments
		os.system('bwa index -p ' + cur_it_dir + '/conserved_regions/curr_it_crs ' + cur_it_dir + '/conserved_regions/org_conserved_regions_it' + str(curr_it) + '.fasta');
		
		time.sleep(5)
		os.system('/24-2/home/hamid/Tools/bwa/bwa mem -a -t 4 -w 1 -B 1000 -O 1000 -E 1000 -U 10000 -T 55 ' + cur_it_dir + '/conserved_regions/curr_it_crs -p ' + cur_it_dir + '/slice0/reads_slc0.fastq >> ' + cur_it_dir + '/conserved_regions/mpd_to_curr_it_crs.sam');
		
		time.sleep(5)
		for curr_slc in range(2,num_slc):
			time.sleep(5)
			os.system('/24-2/home/hamid/Tools/bwa/bwa mem -a -t 4 -w 1 -B 1000 -O 1000 -E 1000 -U 10000 -T 55 ' + cur_it_dir + '/conserved_regions/curr_it_crs -p ' + cur_it_dir + '/slice' + str(curr_slc) + '/reads_slc' + str(curr_slc) + ".fastq | grep -v '^@PG' | grep -v '^@SQ' >> " + cur_it_dir + '/conserved_regions/mpd_to_curr_it_crs.sam');

		time.sleep(25)
		os.system('samtools view -h -S -q 55 ' + cur_it_dir + '''/conserved_regions/mpd_to_curr_it_crs.sam > ''' + cur_it_dir + '/conserved_regions/fltrd_mpd_to_curr_it_crs.sam');
		time.sleep(25)
		os.system('samtools view -S -f 0x0040 ' + cur_it_dir + '''/conserved_regions/fltrd_mpd_to_curr_it_crs.sam | awk '$3 != "*" {print $0}' | awk '$0 ~"NM:i:0"' | awk '$0 ~"XS:i:0"' | awk '$7 != "=" {print $1,$3,1}' > ''' + cur_it_dir + '/conserved_regions/all_ci_orph_rds.txt');
		time.sleep(25)
		os.system('samtools view -S -f 0x0080 ' + cur_it_dir + '''/conserved_regions/fltrd_mpd_to_curr_it_crs.sam | awk '$3 != "*" {print $0}' | awk '$0 ~"NM:i:0"' | awk '$0 ~"XS:i:0"' | awk '$7 != "=" {print $1,$3,2}' >> ''' + cur_it_dir + '/conserved_regions/all_ci_orph_rds.txt');
		##

	return cr_records

def mrg_cnsrvd_rgns(working_dir, curr_it, curr_it_cnsr, min_ovl_ln_0=200, min_ovl_ln_1=100, min_num_brdgs_1=8, min_ovl_ln_2=50, min_num_brdgs_2=80):
	print "		merging the CRs ... "

	cur_it_dir = working_dir + '/it' + str(curr_it)
	prev_it_dir = working_dir + '/it' + str(curr_it-1)

	# add all conserved regions found so far to new conserved regions
	if curr_it > 1:
		prev_it_cnsr = []
		handle = open(prev_it_dir + "/conserved_regions/mrgd_conserved_regions.fasta", "rU")
		for record in SeqIO.parse(handle, "fasta") :
			prev_it_cnsr.append(record)
		handle.close()
		cnsr = curr_it_cnsr + prev_it_cnsr
	else:
		cnsr = curr_it_cnsr

	# orphan reads maped to new conserved segments and orphan reads mapped to prev iteration conserved region
	orphan_reads = [line.split() for line in open(cur_it_dir + '/conserved_regions/all_ci_orph_rds.txt')] 
	if curr_it > 1:
		orphan_reads = orphan_reads + [line.split() for line in open(prev_it_dir + '/conserved_regions/orphan_reads_fnl.txt')]


	counter = 0	
	bfr_cr_ln, aftr_cr_ln = len(cnsr), 0
	while bfr_cr_ln != aftr_cr_ln:
		bfr_cr_ln = len(cnsr)
		for curr_item in cnsr:
			for another_item in cnsr:
				if curr_item != another_item:

					if str(another_item.seq) in str(curr_item.seq):
						for line in orphan_reads:
							if line[1]==another_item.id:
								line[1]=curr_item.id
						cnsr.remove(another_item)
					elif str(curr_item.seq) in str(another_item.seq):
						for line in orphan_reads:
							if line[1]==curr_item.id:
								line[1]=another_item.id
						cnsr.remove(curr_item)
						break
					elif ReverseComplement(str(another_item.seq)) in str(curr_item.seq):
						for line in orphan_reads:
							if line[1]==another_item.id:
								line[1]=curr_item.id
						cnsr.remove(another_item )
					elif ReverseComplement(str(curr_item.seq)) in str(another_item.seq):
						for line in orphan_reads:
							if line[1]==curr_item.id:
								line[1]=another_item.id
						cnsr.remove(curr_item)
						break
					else: 
						lngovrlp = get_longest_overlap(str(curr_item.seq),str(another_item.seq))
						lngovrlprc = get_longest_overlap(ReverseComplement(str(curr_item.seq)),str(another_item.seq))

	 					if (len(lngovrlp) > min_ovl_ln_0) or (len(lngovrlp) > min_ovl_ln_1 and num_bridge_reads(curr_item.id,another_item.id,orphan_reads) > min_num_brdgs_1) or (len(lngovrlp) > min_ovl_ln_2 and num_bridge_reads(curr_item.id,another_item.id,orphan_reads) > min_num_brdgs_2):
							counter = counter + 1
							new_id = 'it'+ str(curr_it) + '_mrgd_cr_' + str(counter)
							mrgd_record = SeqRecord(Seq(merge(str(curr_item.seq),str(another_item.seq),lngovrlp)),id=new_id)
							cnsr.append(mrgd_record)

							for line in orphan_reads:
								if line[1]==curr_item.id or line[1]==another_item.id:
									line[1]=mrgd_record.id
							cnsr.remove(curr_item)
							cnsr.remove(another_item)
							break
							
						elif (len(lngovrlprc) > min_ovl_ln_0) or (len(lngovrlprc) > min_ovl_ln_1 and num_bridge_reads(curr_item.id,another_item.id,orphan_reads) > min_num_brdgs_1) or (len(lngovrlp) > min_ovl_ln_2 and num_bridge_reads(curr_item.id,another_item.id,orphan_reads) > min_num_brdgs_2):
							counter = counter + 1
							new_id = 'it'+ str(curr_it) + '_mrgd_cr_' + str(counter)
							mrgd_record = SeqRecord(Seq(merge(ReverseComplement(str(curr_item.seq)),str(another_item.seq),lngovrlprc)),id=new_id)
							cnsr.append(mrgd_record)

							for line in orphan_reads:	
								if line[1]==curr_item.id or line[1]==another_item.id:
									line[1]=mrgd_record.id
							cnsr.remove(curr_item)
							cnsr.remove(another_item)
							break
		aftr_cr_ln = len(cnsr)

	return cnsr


def find_playful_reads(working_dir, num_slc, curr_it, guard_bps=50):
	print "	Step 3: Finding reads which should go to the next iteration..."

	cur_it_dir = working_dir + '/it' + str(curr_it)

	cr_records = []
	for seq_record in SeqIO.parse(cur_it_dir + '/conserved_regions/mrgd_conserved_regions.fasta', "fasta"):
		record = SeqRecord(seq_record.seq[guard_bps:len(seq_record.seq) - guard_bps],id=seq_record.id)
		cr_records.append(record)
	output_handle = open(cur_it_dir + '/conserved_regions/shrtn_mrgd_conserved_regions.fasta', "w")
	SeqIO.write(cr_records, output_handle, "fasta")
	output_handle.close()

	time.sleep(5)
	os.system('bwa index -p ' + cur_it_dir + '/conserved_regions/cnsrvd_frgmnts ' + cur_it_dir + '/conserved_regions/shrtn_mrgd_conserved_regions.fasta');
	
	time.sleep(5)
	for curr_slc in range(num_slc):
		os.system('/24-2/home/hamid/Tools/bwa/bwa mem -a -t 4 -w 1 -B 1000 -O 1000 -E 1000 -U 10000 -T 55 ' + cur_it_dir + '/conserved_regions/cnsrvd_frgmnts -p ' + cur_it_dir + '/slice' + str(curr_slc) + '/reads_slc' + str(curr_slc) + '.fastq > ' + cur_it_dir + '/slice' + str(curr_slc) + '/algnd_to_cnsrvd_frgmnts.sam');
		time.sleep(5)
		os.system('samtools view -h -S -q 55 '+ cur_it_dir + '/slice' + str(curr_slc) + '/algnd_to_cnsrvd_frgmnts.sam > ' + cur_it_dir + '/slice' + str(curr_slc) + '/fltrd_algnd_to_cnsrvd_frgmnts.sam');


		time.sleep(5)
		if os.path.isfile(cur_it_dir + '/slice' + str(curr_slc) + '/fltrd_algnd_to_cnsrvd_frgmnts.sam'):
			os.system("grep '=' " + cur_it_dir + '/slice' + str(curr_slc) + "/fltrd_algnd_to_cnsrvd_frgmnts.sam | awk '{print $1}' > " + cur_it_dir + '/slice' + str(curr_slc) + '/algnd_rnames.txt');
		else:
			os.system('samtools view -h -S -q 55 '+ cur_it_dir + '/slice' + str(curr_slc) + '/algnd_to_cnsrvd_frgmnts.sam > ' + cur_it_dir + '/slice' + str(curr_slc) + '/fltrd_algnd_to_cnsrvd_frgmnts.sam');
			time.sleep(5)
			os.system("grep '=' " + cur_it_dir + '/slice' + str(curr_slc) + "/fltrd_algnd_to_cnsrvd_frgmnts.sam | awk '{print $1}' > " + cur_it_dir + '/slice' + str(curr_slc) + '/algnd_rnames.txt');
		
		time.sleep(5)
		os.system('samtools view -S -f 0x0040 ' + cur_it_dir + "/slice" + str(curr_slc) + '''/fltrd_algnd_to_cnsrvd_frgmnts.sam | awk '$0 ~"NM:i:0"' | awk '$0 ~"XS:i:0"' | awk '$7 != "=" {print $1,$3,1}' >> ''' + cur_it_dir + "/conserved_regions/orphan_reads_fnl.txt");
		time.sleep(5)
		os.system('samtools view -S -f 0x0080 ' + cur_it_dir + "/slice" + str(curr_slc) + '''/fltrd_algnd_to_cnsrvd_frgmnts.sam | awk '$0 ~"NM:i:0"' | awk '$0 ~"XS:i:0"' | awk '$7 != "=" {print $1,$3,2}' >> ''' + cur_it_dir + "/conserved_regions/orphan_reads_fnl.txt");

		gr_names = []
		gpr_names = []

		time.sleep(5)
		os.system('grep ''^@16'' ' + cur_it_dir + '/slice' + str(curr_slc) + '/reads_slc' + str(curr_slc) + '.fastq | sed ''s/@//g'' > ' + cur_it_dir + '/slice' + str(curr_slc) + '/ar_headers_slc' + str(curr_slc) + '.txt');

		rm = open(cur_it_dir + '/slice' + str(curr_slc) + '/algnd_rnames.txt', 'r')
		for line in rm:
			line = line.strip('\n')
			gr_names.append(line)

		both_pairs = [x for x, y in collections.Counter(gr_names).items() if y > 1]

		for item in both_pairs:
			gpr_names.append(item+'/1')
			gpr_names.append(item+'/2')

		with open(cur_it_dir + '/slice' + str(curr_slc) + '/ar_headers_slc' + str(curr_slc) + '.txt', 'r') as all_headers:
			arns = all_headers.read()
			all_rnames = arns.split('\n')
		while '' in all_rnames:
			all_rnames.remove('')
		bd_names = list(set(all_rnames) - set(gpr_names))
		bd_names.sort()
		all_reads_index = SeqIO.index(cur_it_dir + '/slice' + str(curr_slc) + '/reads_slc' + str(curr_slc) + '.fastq', "fastq")
		bd_reads_index = (all_reads_index[r] for r in bd_names)
		SeqIO.write(bd_reads_index, cur_it_dir + '/slice' + str(curr_slc) + '/playful_reads_slc' + str(curr_slc) + '.fastq', "fastq") 


def prepare_next_it(working_dir, num_slc, curr_it):
	next_it = curr_it + 1

	for i in range(num_slc):
		if not os.path.exists(working_dir + '/it' + str(next_it) + '/slice' + str(i)):
			os.makedirs(working_dir + '/it' + str(next_it) + '/slice' + str(i))

	if not os.path.exists(working_dir + '/it' + str(next_it) + '/conserved_regions'):
		os.makedirs(working_dir + '/it' + str(next_it) + '/conserved_regions')

	for i in range(num_slc):
		os.system('cp ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/playful_reads_slc' + str(i) + '.fastq ' + ' ' + working_dir + '/it' + str(next_it) + '/slice' + str(i) + '/reads_slc' + str(i) + '.fastq')

	for i in range(num_slc):
		os.system('rm  ' + working_dir + '/it' + str(curr_it)  + '/conserved_regions/mpd_to_curr_it_crs.sam')

	for i in range(num_slc):
		os.system('rm  ' + working_dir + '/it' + str(curr_it)  + '/conserved_regions/fltrd_mpd_to_curr_it_crs.sam')

	for i in range(num_slc):
		os.system('rm ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/algnd_to_cnsrvd_frgmnts.sam')

	for i in range(num_slc):
		os.system('rm ' + working_dir + '/it' + str(curr_it) + '/slice' + str(i) + '/Sequences')




def main(argv):

	wd = '' #working directory
	irs = '' #input reads (interleved paired-end reads)
	rel = 0  #reference estimated length
	dc = 0 #desired coverage of each slice
	init_mcl = 0

	try:
		opts, args = getopt.getopt(argv,"r:i:c:o:n:",["reflen=","ifile=","descov=","odir=","initminlen="])
	except getopt.GetoptError:
		print 'slicembler.py -r <reference estimated length> -i <input read set> -c <coverage of each slice> -n <initial FOS minimum length> -o <output directory>'
		sys.exit(2)
	for opt, arg in opts:
		print opt
		print arg
		if opt == '-h':
			print 'slicembler.py -r <reference estimated length> -i <input read set> -c <coverage of each slice> -n <initial FOS minimum length> -o <output directory>'
			sys.exit()
		elif opt in ("-o", "--odir"):
			wd = arg
		elif opt in ("-r", "--reflen"):
			rel = arg
		elif opt in ("-i", "--ifile"):
			irs = arg
		elif opt in ("-c", "--descov"):
			dc = arg
		elif opt in ("-n", "--initminlen"):
			init_mcl = arg

	#### other important parameters: k_mer -> assemble_slices

	#### initialization
	ci = 1
	ended = False
	min_cr = 200
	ttl_len_prev_cnsrv_rgns = 0
	num_prev_cnsrv_rgns = 0
	mrgd_cnsrv_rgns = []
	if not os.path.exists(wd):
		os.makedirs(wd)

	#### Step 0: slicing the reads based on the desired coverage for each back
	num_slc = do_slicing(working_dir=wd, input_reads=irs, ref_est_len=rel, des_cov=dc)

	#### the main loop. Run the pipeline and make new iterations untill there no new CR is generated or the new CRs are very small
	while not(ended) and (ttl_len_prev_cnsrv_rgns < (1 * rel)): 

		#### assemble each of the slices
		assemble_slices(working_dir=wd, num_slc=num_slc, curr_it=ci, k_mer=69)

		dcs = num_slc
		mcl = init_mcl

		while not(ended) and (len(mrgd_cnsrv_rgns) == num_prev_cnsrv_rgns) and (total_len(mrgd_cnsrv_rgns) == ttl_len_prev_cnsrv_rgns):
			print "	Step 2: Finding conserved regions..."

			# find conserved regions 
			curr_orig_cnsrv_rgns = find_cnsrvd_rgns(working_dir=wd, curr_it=ci, min_cr_len=mcl, num_slc=num_slc, des_cr_sbsts=dcs)
			
			if len(curr_orig_cnsrv_rgns)>0:
				# merging 
				mrgd_cnsrv_rgns = mrg_cnsrvd_rgns(working_dir=wd, curr_it=ci, curr_it_cnsr=curr_orig_cnsrv_rgns , min_ovl_ln_0=100, min_ovl_ln_1=50, min_num_brdgs_1=(dc * num_slc)/1000, min_ovl_ln_2=20, min_num_brdgs_2=(dc * num_slc)/100)

			# relax conditions
			if ((dcs-1) > (num_slc/2)):
				dcs = dcs - 1
			elif (mcl/2) >= min_cr:
				mcl = mcl / 2
				dcs = num_slc
			else:
				ended = True


		ttl_len_prev_cnsrv_rgns = total_len(mrgd_cnsrv_rgns)
		num_prev_cnsrv_rgns = len(mrgd_cnsrv_rgns)

		# write the final set of conserved regions for this iteration
		output_handle = open(wd + '/it' + str(ci)+ '/conserved_regions/mrgd_conserved_regions.fasta', 'w')
		SeqIO.write(mrgd_cnsrv_rgns, output_handle, "fasta")
		output_handle.close()

		if not(ended):
			# find pair reads not mapped to (all) conserveg regions
			find_playful_reads( working_dir=wd, num_slc=num_slc, curr_it=ci, guard_bps=50) 
			# create directories and copy playful reads
			prepare_next_it( working_dir=wd, num_slc=num_slc, curr_it=ci)
			ci = ci + 1

		print "	There are " + str(len(mrgd_cnsrv_rgns)) + "CRs in total at the end of iteration " + str(ci-1)

	print "The pipeline ended after " + str(ci) + " iterations"

if __name__ == "__main__":
   main(sys.argv[1:])
