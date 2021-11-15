import numpy as np 
import pandas as pd
from ..core import ReporterScreen
from _supporting_fn import *

def _base_edit_to_from(start_base: chr = "A"):
	base_map = {"A":"G", "C":"T"}
	return(base_map[start_base])

class GuideEditCounter():
	def __init__(self, **kwargs):
		self.R1_filename = kwargs["R1"]
		self.R2_filename = kwargs["R2"]
		self.base_edited_from = kwargs["edited_base"]
		self.base_edited_to = _base_edit_to_from(self.base_edited_from)
		self.min_average_read_quality = kwargs["min_average_read_quality"]
		self.min_single_bp_quality = kwargs["min_single_bp_quality"]
		self.name = kwargs["name"]

		self.qstart_R1 = kwargs["qstart_R1"]
		self.qend_R1 = kwargs["qend_R1"]
		self.qstart_R2 = kwargs["qstart_R2"]
		self.qend_R2 = kwargs["qend_R2"]

		self.sgRNA_file = kwargs["sgRNA_file"]
		self._set_sgRNA_df()

    	self.screen = ReporterScreen(
    		X = np.zeros((len(self.guides_info_df), 1)), 
    		obs = self.guides_info_df, 
    		var = pd.DataFrame(index = self.name)
    		)

		self.count_reporter_edits = kwargs["count_reporter_edits"]
		if self.count_reporter_edits:
			self.screen.guides["edits"] = dict()
			self.screen.layers["edits"] = np.zeros_like(self.screen.X)
		self.count_edited_alleles = kwargs["count_alleles"]
		if self.count_edited_alleles:
			self.screen.guides["edited_alleles"] = dict()
		self.count_only_bcmatched = NotImplemented

		self.guide_start_sequence = kwargs["guide_start_seq"]
		self.guide_bc = kwargs["guide_bc"]
		if self.guide_bc:
			self.guide_bc_len = kwargs["guide_bc_len"]

		self.offset = kwargs["offset"]
		if self.count_reporter_edits:
			self.output_prefix = kwargs['output_prefix']
			self.reporter_length = kwargs['reporter_length']

		self.n_total_reads = _read_count_match(self.R1_filename, self.R2_filename)
		self.database_id = _get_database_id(self.name)

		self.keep_intermediate = kwargs[keep_intermediate]
		self.output_dir = os.path.join(os.path.abspath(args.output_folder), 
			'CRISPRessoCount_on_%s' % self.database_id)
        self.log_filename=self._jp('CRISPRessoCount_RUNNING_LOG.txt')
        self._write_start_log()



    def _set_sgRNA_df(self):
	    with open(self.sgRNA_filename) as infile:
	        sgRNA_df = pd.read_csv(infile)
	        if not ('name' in sgRNA_df.columns and 'sequence' in sgRNA_df.columns):
	            raise InputFileError("Input gRNA info file doesn't have the column 'name' or 'sequence'.")
	        if self.count_only_bcmatched and 'barcode' not in sgRNA_df.columns:
	        	raise InputFileError("Input gRNA info file doesn't have the column 'barcode'.")
	        sgRNA_df = sgRNA_df.set_index("name")
	    self.guides_info_df = sgRNA_df
	    self.guide_lengths = sgRNA_df.sequence.map(lambda s: len(s)).unique()


    def check_filter_fastq(self):
    	self.filtered_R1_filename = self._jp(os.path.basename(self.R1_filename).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz')
    	self.filtered_R2_filename = self._jp(os.path.basename(self.R2_filename).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz')
    	if path.exists(filtered_R1_filename) and path.exists(filtered_R2_filename):
            warn("Using preexisting filtered file")
        else:
        	self._check_names_filter_fastq()

    def get_counts(self):
        infile_R1 = _get_fastq_handle(self.R1_filename)
    	infile_R2 = _get_fastq_handle(self.R2_filename)
	
        self.nomatch_R1_filename = R1_filename.replace(".fastq", "_nomatch.fastq")
        self.nomatch_R2_filename = R2_filename.replace(".fastq", "_nomatch.fastq")
        self.semimatch_R1_filename = R1_filename.replace(".fastq", "_semimatch.fastq")
        self.semimatch_R2_filename = R2_filename.replace(".fastq", "_semimatch.fastq")

    	if not self.count_only_bcmatched:
    		# count X
    		self._get_guide_counts_bcmatch()
		else: # count both bc matched & unmatched guide counts
			self._get_guide_counts_bcmatch_semimatch()


    def _gRNA_eq(self, guide, observed):
	    if len(observed) != len(guide): return(False)
	    for observed_nt, guide_nt in zip(observed, guide):
	        if (observed_nt == guide_nt) or ((guide_nt, observed_nt) == (self.edit_start, self.edit_end)):
	            continue
	        else:
	            return(False)
	    return(True)


	def _get_guide_counts_bcmatch(self):
		NotImplemented


	def _get_guide_counts_bcmatch_semimatch(self):
		self.screen.layers["X"] = np.zeros_like((self.screen.X))
		R1_iter, R2_iter = self._get_fastq_iterators()

		outfile_R1_nomatch, outfile_R2_nomatch = self._get_fastq_handle("nomatch")
		outfile_R1_semimatch, outfile_R2_semimatch = self._get_fastq_handle("semimatch")
		outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc = self._get_fastq_handle("duplicate_wo_barcode")
		outfile_R1_dup, outfile_R2_dup = self._get_fastq_handle("duplicate")

		for i, (r1, r2) in tqdm(enumerate(zip(R1_iter, R2_iter)), total = self.n_reads_after_filtering):
	        _, R1_seq, _ = r1
	        _, R2_seq, _ = r2
	        _fastq_iter_to_text

	        bc_match, semimatch = self._match_read_to_sgRNA_bcmatch_semimatch(R1_seq, R2_seq)
	        
	        if len(bc_match) == 0:
	        	if len(semimatch) == 0:# no guide match
	        		_write_paired_end_reads(r1, r2, outfile_R1_nomatch, outfile_R2_nomatch)
	        	elif len(semimatch) >= 2:# Duplicate match if w/o barcode	        		
	        		_write_paired_end_reads(r1, r2, outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc)
	        	else: # guide match with no barcode match	        		
	        		self.screen["X_bcmatch"][semimatch[0], 0] += 1
    		elif len(bc_match) >= 2:
    			_write_paired_end_reads(r1, r2, outfile_R1_dup, outfile_R2_dup)
    		else: # unique barcode match
    			matched_guide_idx = bc_match[0]
    			self.screen.X[matched_guide_idx, 0] += 1

    			if self.count_reporter_edits:
    				self.count_reporter_edits()
    					# TBD: what if reporter seq doesn't match barcode & guide?
    					ref_reporter_seq = self.screen.guides.Reporter[matched_guide_idx]
    					read_reporter_seq = self._get_reporter_seq(R1_seq, R2_seq)

    					allele = _get_edited_allele(
    						ref_reporter_seq, 
    						read_reporter_seq, 
    						self.screen.guides.offset[matched_guide_idx]
    						)

    					if self.count_edited_alleles:
	    					self._write_allele(matched_guide_idx, allele)
    					self._write_edits(matched_guide_idx, allele)

	        
	def _write_allele(self, guide_idx: int, allele: Allele):
		if not allele in self.secreen.guides[guide_idx, :]["edited_alleles"].keys():
			self.secreen.guides[guide_idx, :]["edited_alleles"][alelle] = 1
		else:
			self.secreen.guides[guide_idx, :]["edited_alleles"][alelle] += 1


	def _write_edits(self, guide_idx: int, allele: Allele):
		for edit in allele.edits:
			if not edit in self.screen.guides[guide_idx, :]["edits"].keys():
				self.screen.guides[guide_idx, :]["edits"][edit] = 1
			else:
				self.screen.guides[guide_idx, :]["edits"][edit] += 1
			if edit.rel_pos == self.guides.target_pos:
				self.screen.edits[guide_idx, 0] += 1


    def get_reporter_seq(self, R1_seq, R2_seq):
    	"""This can be edited by user based on the read construct."""
    	reporter_seq = revcomp(R2_seq[self.guide_bc_len:(self.guide_bc_len + self.reporter_length)])
    	return(reporter_seq)


	def _match_read_to_sgRNA_bcmatch_semimatch(self, R1_seq, R2_seq):
		# This should be adjusted for each experimental recipes.
		bc_match = []
		semimatch = []
	    for i, (guide_name, row) in enumerate(self.guides_info_df.iterrows()):
	        if self._gRNA_eq(row["sequence"], seq) :
	            if ref_barcode == guide_bc: 
	            	bc_match.append(i)
	                self.screen.layers[bcmatch_layer][i, 0] += 1
	            else:
	            	semimatch.append(i)
	            	self.screen.layers[semimatch_layer][i, 0] += 1
    	return((bc_match, semimatch))
	    

	def _get_guide_position_seq_of_read(self, seq):
		guide_start_idx = self._get_guide_start_idx(seq)
		if guide_start_idx == -1: return None

		return([seq[guide_start_idx:(guide_start_idx + guide_length)]
			for guide_length in self.guide_lengths])

	def _get_guide_start_idx(self, seq):
		start_seq_idx = seq.find(self.guide_start_sequence)
		if start_seq_idx == -1: return -1
		return(start_seq_idx + len(self.guide_start_sequence))

	def get_gRNA_barcode(self, R1_seq, R2_seq):
		# This can be adjusted for different construct design.
    	return(revcomp(R2_seq[:guide_bc_len]))

    def _get_fastq_handle(self, out_type: Literal["semimatch", "nomatch"] = None):
		elif intermediate == "semimatch":
			R1_handle = _get_fastq_handle(self.semimatch_R1_filename)
			R2_handle = _get_fastq_handle(self.semimatch_R2_filename)
		elif intermediate == "nomatch":
			R1_handle = _get_fastq_handle(self.nomatch_R1_filename)
			R2_handle = _get_fastq_handle(self.nomatch_R2_filename)
		else:
			raise ValueError("Invalid value '{}' for argument out_type".format(out_type))
		return((R1_handle, R2_handle))


	def _get_fastq_iterators(self, intermediate: Literal["semimatch", "nomatch"] = None):
		R1_handle = _get_fastq_handle(self.R1_filename)
		R2_handle = _get_fastq_handle(self.R2_filename)

		R1_iterator = FastqGeneralIterator(R1_handle)
		R2_iterator = FastqGeneralIterator(R2_handle)

		return((R1_iterator, R2_iterator))


    def _get_seq_records(self):
        R1_handle = _get_fastq_handle(self.R1_filename)
        R2_handle = _get_fastq_handle(self.R2_filename)
        R1 = list(SeqIO.parse(R1_handle, "fastq"))
        R2 = list(SeqIO.parse(R2_handle, "fastq"))
        R1_handle.close()
        R2_handle.close()
        return((R1, R2))


    def _check_names_filter_fastq(self, filter_by_qual = False):
        if self.min_bp_quality > 0 or self.min_single_bp_quality > 0:
            info('Filtering reads with average bp quality < {} \
                and single bp quality < {} ... in up to position {} of R1 and {} of R2'.format(
                	self.min_average_read_quality,
                    self.min_single_bp_quality))

        if self.qend_R1 > 0 or self.qend_R2 > 0:
            info("In the filtering, bases up to position {} of R1 and {} of R2 are considered.".format(qend_R1, qend_R2))

        R1, R2 = self._get_seq_records()

        _check_readname_match(R1, R2)
        if filter_by_qual:
        	self.n_reads_after_filtering = _filter_read_quality(R1, R2)
        	if self.n_reads_after_filtering == 0:
        		raise NoReadsAfterQualityFiltering('No reads in input or no reads survived the average or single bp quality filtering.')
        	else:
        		info('Number of reads in input:%d\tNumber of reads after filtering:%d\n' % (N_READS_INPUT, N_READS_AFTER_PREPROCESSING))
    	else:
    		self.n_reads_after_filtering = self.n_total_reads



    def _filter_read_quality(self, R1 = None, R2 = None) -> int:
	    R1_filtered = gzip.open(self.filtered_R1_filename, 'w+')
        R2_filtered = gzip.open(self.filtered_R2_filename, 'w+')
        
        if R1 is None or R2 is None:
        	R1, R2 = _get_seq_records()

        n_reads_after_filtering = 0
        for i in range(len(R1)):
            R1_record = R1[i]
            R2_record = R2[i]

            if filter_by_qual:
                R1_quality_pass = _read_is_good_quality(R1_record, self.min_bp_quality, self.min_single_bp_quality, self.qend_R1)
                R2_quality_pass = _read_is_good_quality(R2_record, self.min_bp_quality, self.min_single_bp_quality, self.qend_R2)

                if R1_quality_pass and R2_quality_pass:
                	n_reads_after_filtering += 1
                    R1_filtered.write(record.format('fastq'))
                    R2_filtered.write(record.format('fastq'))
        return(n_reads_after_filtering)


    def _write_start_log(self):
    	with open(self.log_filename, "w+") as outfile:
    		outfile.write('[Command used]:\nCRISPRessoCount %s\n\n[Execution log]:\n' % ' '.join(sys.argv))


	def _jp(self, filename):
		return(os.path.join(self.output_dir,filename))

	def _get_database_name(self):
	    get_name_from_fasta = lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','').replace("_R1",'')
    
	    if not self.name:
	        database_id ='%s' % get_name_from_fasta(self.R1_filename)    
	    else:
	        database_id = self.name
	    return(database_id)