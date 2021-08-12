# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 13:17:59 2017

@author: Luca Pinello
Edited by Jayoung Ryu on Tue Aug 10 11:01 2021
"""

import os
import gzip
import argparse
import sys
import gzip
import subprocess as sb
from collections import defaultdict
import unicodedata
import re

import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info


_ROOT = os.path.abspath(os.path.dirname(__file__))


####Support functions###

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)
 


def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoCount!' % library_name)
                sys.exit(1)

def slugify(value): #adapted from the Django project
    
    value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '_', value).strip())
    value = unicode(re.sub('[-\s]+', '-', value))
    
    return str(value)

def get_fastq_handle(fastq_filename):
    if fastq_filename.endswith('.gz'):
            fastq_handle=gzip.open(fastq_filename)
    else:
            fastq_handle=open(fastq_filename)
    return(fastq_handle)

def read_is_good_quality(record: SeqIO.record, 
        min_bp_quality = 0, 
        min_single_bp_quality = 0, 
        qend = -1):
    mean_quality_pass = np.array(record.letter_annotations["phred_quality"])[:qend_R1].mean() >= min_bp_quality
    min_quality_pass = np.array(record.letter_annotations["phred_quality"])[:qend_R1].min()>=min_single_bp_quality
    return(mean_quality_pass and min_quality_pass)


def check_read_names_and_filter_quality(R1_filename, 
        R2_filename, 
        output_R1_filename=None, 
        output_R2_filename = None,
        filter_by_qual = False,
        min_bp_quality=20,
        min_single_bp_quality=0, 
        qend_R1 = -1, 
        qend_R2 = -1):

        R1_handle = get_fastq_handle(R1_filename)
        R2_handle = get_fastq_handle(R2_filename)

        if not output_R1_filename:
                output_R1_filename = R1_filename.replace('.fastq', '').replace('.gz', '') + '_filtered.fastq.gz'
                output_R2_filename = R2_filename.replace('.fastq', '').replace('.gz', '') + '_filtered.fastq.gz'

        try:
            if filter_by_qual:
                R1_filtered = gzip.open(output_R1_filename, 'w+')
                R2_filtered = gzip.open(output_R2_filename, 'w+')
            R1 = list(SeqIO.parse(R1_handle, "fastq"))
            R2 = list(SeqIO.parse(R2_handle, "fastq"))

            if len(R1) != len(R2):
                raise ValueError("The number of reads in R1 and R2 file does not match.")

            for i in range(len(R1)):
                R1_record = R1[i]
                R2_record = R2[i]

                if R1_record.name != R2_record.name : 
                    raise InputFASTQError("R1 and R2 read discordance in read {} and {}".format(R1_record.name, R2_record.name))
                if filter_by_qual:
                    R1_quality_pass = read_is_good_quality(R1_record, min_bp_quality, min_single_bp_quality, qend_R1)
                    R2_quality_pass = read_is_good_quality(R2_record, min_bp_quality, min_single_bp_quality, qend_R2)

                    if R1_quality_pass and R2_quality_pass:
                        R1_filtered.write(record.format('fastq'))
                        R2_filtered.write(record.format('fastq'))
        except:
                raise Exception('Error handling the fastq_filtered_outfile')

        if filter_by_qual: 
            return((output_R1_filename, output_R2_filename))
        else:
            return((R1_filename, R2_filename))



def find_wrong_nt(sequence):
    return(list(set(sequence.upper()).difference(set(['A','T','C','G','N']))))


def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     return(int(float(p.communicate()[0])/4.0))

def mask_sequence_positions(seq: str, pos: np.array) -> str:
    return("".join([seq[i] if i in pos else "N" for i in range(len(seq))]))
        


pd=check_library('pandas')
np=check_library('numpy')
Bio=check_library('Bio')
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from itertools import pairwise

###CLASS FOR MASKED SEQUENCE FOR GRNA COMPARISON####
class masked_sequence(Bio.Seq.Seq):
    def __init__(self, data, mask_base = None):
        super().__init__(data)
        if mask_base:
            self._mask = [self._data[i] == mask_base for i in len(self._data)]

    def __eq__(self, other):
        if isinstance(other, masked_sequence):
            query_nomask = np.nonzero([not (a or b) for a,b in pairwise(self._mask, other._mask)])
            return(all([a[i] == b[i] for x in query_nomask]))

            

###EXCEPTIONS############################
class NTException(Exception):
    pass

class InputFASTQError(Exception):
    pass
class NoReadsAfterQualityFiltering(Exception):
    pass
#########################################

def main():
    try:
        print '  \n~~~CRISPRessoCount~~~'
        print '-Utility to perform sgRNA enumeration from deep sequencing data-'
        print r'''
              )                                             )
             (           ________________________          (
            __)__       | __   __            ___ |        __)__
         C\|     \      |/  ` /  \ |  | |\ |  |  |     C\|     \
           \     /      |\__, \__/ \__/ | \|  |  |       \     /
            \___/       |________________________|        \___/
        '''
        
        
        print'\n[Luca Pinello 2017, send bugs, suggestions or *green coffee* to lucapinello AT gmail DOT com]\n\n',
        
        
        parser = argparse.ArgumentParser(description='CRISPRessoCount parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('--R1', type=str,  help='fastq file for read 1', required=True,default='Fastq filename' )
        parser.add_argument('--R2', type=str,  help='fastq file for read 2, sorted as the same name order as in --R1 file.', required=True, default='Fastq filename' )
        parser.add_argument('-b', '--edited_base', type = str, help = 'For base editors, the base that should be ignored when matching the gRNA sequence')
        parser.add_argument('--guide_start', type = int, help = "Starting position of guide in R1", type = str, required = True, default = 27)
        parser.add_argument('--guide_end', type = int, help = "Ending positon of guide in R1 (not inclusive)", type = str, required = True, default = 47)

        #optional
        parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
        parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
        parser.add_argument('-f','--sgRNA_file', type=str,  help='''sgRNA description file. The format requires three columns: gRNA, Reporter, gRNA_barcode.''')
        parser.add_argument('-n','--name',  help='Output name', default='')
        parser.add_argument('-o','--output_folder',  help='', default='')
        parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
        parser.add_argument('--qstart_R1', help='Start position of the read when filtering for quality score of the read 1', type = int, default = 0)
        parser.add_argument('--qend_R1', help = 'End position of the read when filtering for quality score of the read 1', type = int, default = -1)
        parser.add_argument('--qstart_R2', help = 'Same as qstart_R1, for read 2 fastq file')
        parser.add_argument('--qend_R2', help = 'Same as qstart_R2, for read 2 fastq file')
        parser.add_argument('--guide_bc_len', help = 'Guide barcode sequence length at the beginning of the R2', type = str, default = 0)

        args = parser.parse_args()
        
        
        #check files
        check_file(args.R1)
        check_file(args.R2)

        # Get read length
        fastq_handle = get_fastq_handle(args.R1)
        for record in SeqIO.parse(fastq_handle, "fastq"):
            read_length = len(record)
            break

        if args.sgRNA_file:
            check_file(args.sgRNA_file)
        
        # Edited base should be one of A/C/T/G
        if args.edited_base.upper() not in ['A', 'C', "T", "G"]:
            raise ValueError("The edited base should be one of A/C/T/G, {} provided.".format(args.edited_base))
        else:
            edited_base = args.edited_base.upper()
        info('Using specified edited base: %s' % edited_base)
        

        # Check if positions of guide and quality control is valid
        if not (args.guide start >= 0 and args.guide_end < read_length): 
            raise ValueError("The guide start and end position is not valid. Start position should be nonnegative integer ({} provided), and the ending position should be shorter than the read length {} ({} provided)".format(args.guide_start, read_length, args.guide_end))
        
        if not (args.qstart_R1 < read_length and args.qstart_R1 < read_length):
            raise ValueError("The start position of base quality filter is not nonnegative ({} for R1, {} for R2 provided)".format(args.qstart_R1, args.qstart_R2)
        if not (args.qend_R1 < read_length and args.qend_R1 < read_length):
            raise ValueError("The start position of base quality filter is not nonnegative ({} for R1, {} for R2 provided)".format(args.qstart_R1, args.qstart_R2)
        

        #normalize name and remove not allowed characters
        if args.name:   
            clean_name=slugify(args.name)
            if args.name!= clean_name:
                   warn('The specified name %s contained characters not allowed and was changed to: %s' % (args.name,clean_name))
                   args.name=clean_name
                    
        
        get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','')
        
        if not args.name:
                database_id='%s' % get_name_from_fasta(args.fastq)
        
        else:
                database_id=args.name
        
        
        OUTPUT_DIRECTORY='CRISPRessoCount_on_%s' % database_id
        
        if args.output_folder:
                OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
        
        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoCount_RUNNING_LOG.txt')
        
        
        try:
        
                 os.makedirs(OUTPUT_DIRECTORY)
                 logging.getLogger().addHandler(logging.FileHandler(log_filename))
        
                 with open(log_filename,'w+') as outfile:
                     outfile.write('[Command used]:\nCRISPRessoCount %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
                 info('Creating Folder %s' % OUTPUT_DIRECTORY)
                 info('Done!')
        except:
                 warn('Folder %s already exists.' % OUTPUT_DIRECTORY)


        # Check read name is aligned to be the same in R1 and R2 files.
        # Optionally filter the reads based on the base quality
        if args.min_average_read_quality>0 or args.min_single_bp_quality>0:
            info('Filtering reads with average bp quality < %d and single bp quality < %d ...' % (args.min_average_read_quality,args.min_single_bp_quality))

        filtered_R1_filename, filtered_R2_filename = check_read_names_and_filter_quality(args.R1,
            args.R2, 
            output_R1_filename=_jp(os.path.basename(args.R1).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'),
            output_R2_filename=_jp(os.path.basename(args.R2).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz'),
            min_bp_quality=args.min_average_read_quality,
            min_single_bp_quality=args.min_single_bp_quality, 
            qend_R1 = qend_R1, 
            qend_R2 = qend_R2)

        info('Done!')
        
        #count reads 
        N_READS_INPUT=get_n_reads_fastq(args.R1)
        N_READS_INPUT_R2 = get_n_reads_fastq(args.R2)
        if N_READS_INPUT != N_READS_INPUT_r2 : 
            raise InputFASTQError("Number of reads not matched in two input fastq file") 
        N_READS_AFTER_PREPROCESSING = get_n_reads_fastq(filtered_R1_filename)
        N_READS_AFTER_PREPROCESSING_R2 = get_n_reads_fastq(filtered_R2_filename)
        assert N_READS_AFTER_PREPROCESSING == N_READS_AFTER_PREPROCESSING_R2, "Filtered read numbers are different in R1 and R2"

        if N_READS_AFTER_PREPROCESSING == 0:             
            raise NoReadsAfterQualityFiltering('No reads in input or no reads survived the average or single bp quality filtering.')
        else:
            info('Number of reads in input:%d\tNumber of reads after filtering:%d\n' % (N_READS_INPUT, N_READS_AFTER_PREPROCESSING))
 
        if args.sgRNA_file:
            info('Using guides information in %s' % args.sgRNA_file)
            guides_count=dict()
            guides = dict()
            with open(args.sgRNA_file) as infile:
                sgRNA_df = pd.read_csv(infile)
                if not ('gRNA' in data.columns and 'gRNA_barcode' in data.columns):
                    raise InputFileError("Input gRNA file doesn't have the column 'gRNA' or 'gRNA_barcode'.")
                
                duplicate_gRNA = 0
                for i, guide in enumerate(sgRNA_df["gRNA"]):
                    assert len(guide) == 19 or len(guide) == 20, "Guide length not 19-20bp"
                    
                    guide = guide.replace(args.edited_base, "N")
                    
                    if guide in guides.keys(): 
                        guides[guide] += 1
                    else: guides[guide] = 1

                    gRNA_bc = sgRNA_df["gRNA_barcode"][i]
                    if len(gRNA_bc) != args.guide_bc_len:
                        raise InputFileError("Input gRNA barcode length {} in {}th record doesn't match the provided barcode length {}".format(len(gRNA_bc), i, args.guide_bc_len)
                    
                    gRNA_seq = sgRNA_df["gRNA"][i]
                    mask_pos = np.nonzero(gRNA_seq == args.edited_base)
                    masked_seq = gRNA_seq.replace(edited_base, "N")

                    if gRNA_bc in guides_count.keys():
                        if mask_pos in guides_count[gRNA_bc].keys():
                            if masked_seq in guides_count[gRNA_bc][mask_pos].keys():
                                duplicate_gRNA += 1
                            else:
                                guides_count[gRNA_bc][mask_pos][masked_seq] = 0
                        else:
                            guides_count[gRNA_bc][mask_pos] = {masked_seq: 0}                    
                    else: 
                        guides_count[gRNA_bc] = {mask_pos : {masked_seq:0}}

                    # Write masked guide counts
                    f = open(_jp("masked_guide_counts"), "w")
                    for key, val in guides.items():
                        f.write("{}\t{}\n".format(key, val))
                    f.close()
                    info("{} guides with indistinguishable masked sequence + barcode combinations".format(duplicate_gRNA))


        else:
            info('No guide information file specified, counting all the guides')
            guides_count=defaultdict(lambda:0)
            
        
        if filtered_R1_filename.endswith('.gz'):
            infile_R1=gzip.open(filtered_R1_filename)
        else:
            infile_R1=open(filtered_R1_filename)
        if filtered_R2_filename.endswith('.gz'):
            infile_R2 = gzip.open(filtered_R2_filename)
        else:
            infile_R2 = open(filtered_R2_filename)
        
            
        
        
        info('Counting sgRNAs...')
        N_READS=0
        
        iterator_R1 = FastqGeneralIterator(infile_R1)
        iterator_R2 = FastqGeneralIterator(infile_R2)

        for i in range(N_READS_AFTER_PREPROCESSING):
            _, R1_seq, _ = next(iterator_R1)
            _, R2_seq, _ = next(iterator_R2)
            gRNA_barcode = R2_seq[:args.guide_bc_len]
            gRNA_seq = R1_seq[guide_start:guide_end]
            assert len(gRNA_seq) == 20, "gRNA sequence is not of length 20bp"
            
            for mask in guides_count[gRNA_barcode].keys():
                masked_seq = gRNA_seq





        while infile.readline():
            read_seq=infile.readline().strip()
            infile.readline()
            infile.readline()
            
            N_READS+=1
            
            tracrRNA_idx=read_seq.find(args.tracrRNA)
            
            if tracrRNA_idx>=0:
                guide_seq=read_seq[tracrRNA_idx-args.guide_length: tracrRNA_idx]
        
        
                if args.sgRNA_file and not guide_seq in guides_count:
                    pass
                else:
                    guides_count[guide_seq]+=1
                        
        infile.close()
        info('Done!')    

        info('Writing output table...')
        df_guide_counts=pd.Series(guides_count,name='Read_Counts').to_frame()
        df_guide_counts.index.name='Guide_Sequence'
        df_guide_counts['Read_%']=df_guide_counts['Read_Counts']/N_READS*100
        df_guide_counts['RPM']=df_guide_counts['Read_Counts']/N_READS*1000000
        df_guide_counts.head()
        
        df_guide_counts.sort_values(by='Read_Counts',ascending=False).to_csv(_jp('CRISPRessoCount_%s_on_%s.txt' % \
            ('only_ref_guides' if args.sgRNA_file else 'no_ref_guides',args.fastq)),sep='\t')
        
        info('Done!')
        
        if not args.keep_intermediate:
            info('Removing Intermediate files...')
            
            files_to_remove=[]
                             
            if args.min_average_read_quality>0 or args.min_single_bp_quality>0:
                files_to_remove.append(processed_output_filename)
                
            for file_to_remove in files_to_remove:
                         try:
                                 if os.path.islink(file_to_remove):
                                     os.unlink(file_to_remove)
                                 else:                             
                                     os.remove(file_to_remove)
                         except:
                                 warn('Skipping:%s' %file_to_remove)    
    
 
        info('All Done!')
        print r'''
              )                                             )
             (           ________________________          (
            __)__       | __   __            ___ |        __)__
         C\|     \      |/  ` /  \ |  | |\ |  |  |     C\|     \
           \     /      |\__, \__/ \__/ | \|  |  |       \     /
            \___/       |________________________|        \___/
        '''
        sys.exit(0)
    
    except Exception as e:
        error('\n\nERROR: %s' % e)
        sys.exit(-1)

if __name__ == '__main__':
    main()                       
