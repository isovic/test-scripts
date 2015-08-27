#! /usr/bin/python

import os;
import sys;
import subprocess;
import fastqparser;

GOLDEN_BUNDLE_SRC = '../../golden-bundle/src/';

def execute_command(command):
	sys.stderr.write('Executing command: "%s"\n' % command);
	subprocess.call(command, shell=True);

def process_sam(sam_file, reference_path, vcf_mutations_path, coverage_threshold):
	# execute_command('%s/analyzesam.py file 1 %s %s %s' % (GOLDEN_BUNDLE_SRC, sam_file, reference_path, reads_path)); 			# This calculates only the consensus.
	sam_name = os.path.splitext(os.path.basename(sam_file))[0];
	vcf_prefix = 'data/variant/ecoliR7.3/%s' % (sam_name);
	vcf_file = '%s-cov_%d.variant.vcf' % (vcf_prefix, coverage_threshold);
	execute_command('%s/consensus_from_mpileup.py %s %d %s %s' % (GOLDEN_BUNDLE_SRC, reference_path, coverage_threshold, vcf_prefix, sam_file)); 			# This calculates only the consensus.
	process_vcf(vcf_file, vcf_mutations_path);

def process_and_left_align_sam(sam_file, reference_path, vcf_mutations_path, coverage_threshold):
	# Convert to BAM and index it.
	bam_prefix = os.path.splitext(sam_file)[0];
	bam_leftaligned_prefix = bam_prefix + '-leftaligned';
	# execute_command('samtools view -bS %s | samtools sort - %s' % (sam_file, bam_prefix));
	# execute_command('samtools index %s.bam %s.bam.bai' % (bam_prefix, bam_prefix));
	# Leftalign the BAM file.
	execute_command('cat %s.bam | ./bamleftalign -f %s > %s.bam' % (bam_prefix, reference_path, bam_leftaligned_prefix));
	execute_command('samtools index %s.bam %s.bam.bai' % (bam_leftaligned_prefix, bam_leftaligned_prefix));

	process_sam(bam_leftaligned_prefix + '.bam', reference_path, vcf_mutations_path, coverage_threshold);

def leftalign(bam_file):
# ### Left align the BAM file:
# cat $SAM_FILTERED-sorted.bam | ./bamleftalign -f $REFERENCE_PATH > $BAM_FINAL.bam
# samtools index $BAM_FINAL.bam
	pass;

def process_vcf(vcf_file, vcf_mutations_path):
	if (vcf_mutations_path.endswith('.gz') == False and os.path.exists(vcf_mutations_path + '.gz') == False):
		execute_command('bgzip -c %s > %s.gz' % (vcf_mutations_path, vcf_mutations_path));
		execute_command('tabix -p vcf %s.gz' % (vcf_mutations_path));
	if (vcf_mutations_path.endswith('.gz') == False):
		vcf_mutations_path += '.gz';
		# print vcf_mutations_path;

	execute_command('bgzip -c %s > %s.gz' % (vcf_file, vcf_file));
	execute_command('tabix -p vcf %s.gz' % (vcf_file));
	# execute_command('bamsurgeon/etc/evaluator.py -v %s.gz -t %s -m SNV' % (vcf_file, vcf_mutations_path));
	fn_file = '%s-fn.vcf' % (os.path.splitext(vcf_file)[0]);
	fp_file = '%s-fp.vcf' % (os.path.splitext(vcf_file)[0]);
	execute_command('bamsurgeon/etc/evaluator.py -v %s.gz -t %s -m SNV --fn %s --fp %s' % (vcf_file, vcf_mutations_path, fn_file, fp_file));


def vcf_extract_positions(vcf_file):
	fp = open(vcf_file, 'r');
	ret = [];
	for line in fp:
		if (len(line.strip()) == 0 or line[0] == '#'):
			continue;
		split_line = line.split();
		chrom = split_line[0];
		pos = int(split_line[1]);
		ref = split_line[3];
		alt = split_line[4];
		info = split_line[7];
		ret.append([chrom, pos, ref, alt, info]);
	fp.close();
	return ret;

def get_kmers_from_positions(fastq_file, pos_list, k):
	kmers = [];

	[headers, seqs, quals] = fastqparser.read_fastq(fastq_file);
	header_hash = {};
	i = 0;
	while (i < len(headers)):
		header_hash[headers[i]] = i;
		header_hash[headers[i].split()[0]] = i;
		i += 1;

	num_homo = 0;

	i = 0;
	for pos_item in pos_list:
		i += 1;
		chrom = pos_item[0];
		pos = pos_item[1] - 1;
		ref = pos_item[2];
		alt = pos_item[3];
		info = pos_item[4];
		try:
			seq = seqs[header_hash[chrom]];
		except Exception, e:
			sys.stderr.write(str(e) + '\n');
			continue;
		# kstart = (pos - 1) - k/2;
		k_before = k;
		k_after = k;
		# klen = 
		# kend = (pos - 1) + k/2;
		# kmer_before = seq[pos-k_before:pos] if (k_before <= pos) else (' ' * (k_before - pos) + seq[0:pos]);
		# kmer_after = seqs[(pos + 1):(pos+1+k_after)] if ((pos+1+k_after) >= len(seq)) else (seq[(pos+1):len(seq)] + (' ' * (pos+1+k_after - len(seq))));
		kmer_before = seq[pos-k:pos];
		kmer_after = seq[(pos+1):(pos+1+k)];
		kmer = kmer_before + '_' + seq[pos] + '_' + kmer_after;
		kmers.append([kmer, chrom, pos]);

		kmer_ref = kmer_before + '_' + ref + '_' + kmer_after;
		kmer_alt = kmer_before + '_' + alt + '_' + kmer_after;

		num_homo += 1 if (kmer_before[-1] == ref or kmer_after[0] == ref) else 0;

		# if (kmer_before[-1] == ref or kmer_after[0] == ref):
		# 	sys.stdout.write('\th [%d] %s\t%s\t%s\t%d\t%s\n' % (i, kmer, kmer_ref, kmer_alt, (pos + 1), info));
		# else:
		sys.stdout.write('[%d] %s\t%s\t%s\t%d\t%s\n' % (i, kmer, kmer_ref, kmer_alt, (pos + 1), info));

	sys.stdout.write('\n');
	sys.stdout.write('Num homopolimer bases: %d\n' % (num_homo));

	return kmers;


def copy_vcf(path_on_genome):
	execute_command('scp -P 2222 isovic@genome.zesoi.fer.hr:%s data/variant/ecoliR7.3/' % (path_on_genome))
	vcf_file = 'data/variant/ecoliR7.3/%s' % (os.path.basename(path_on_genome));
	return vcf_file;

def main():
	sam_file = 'data/alignment/ecoliR7.3/graphmap-params_final_20150417_mutref.sam';
	sam_file = 'data/alignment/ecoliR7.3/LAST-q1.sam';
	sam_file = 'data/alignment/ecoliR7.3/BWAMEM-real_nanopore.sam';

	reference_path = 'data/reference/mutated_ecoli.fa';
	vcf_mutations_path = 'data/reference/rev_mutants.vcf.gz';
	# reads_path = '/home/isovic/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/reads/ecoliR7.3.fasta';

	# process_sam(sam_file, reference_path, vcf_mutations_path, 20);
	# process_sam(sam_file, reference_path, vcf_mutations_path, 0);
	# process_and_left_align_sam(sam_file, reference_path, vcf_mutations_path, 20);
	# process_and_left_align_sam(sam_file, reference_path, vcf_mutations_path, 0);

	### Filter only unique best alignments and process only those:
	# sam_file_uniquebest = os.path.splitext(sam_file)[0] + '-uniquebest.sam';
	# execute_command('%s/samfilter.py uniquebest %s %s' % (GOLDEN_BUNDLE_SRC, sam_file, sam_file_uniquebest)); 			# This calculates only the consensus.
	# process_sam(sam_file_uniquebest, reference_path, vcf_mutations_path, 0);



	### Process only 2d reads:
	# sam_file_2d = os.path.splitext(sam_file)[0] + '-2d.sam';
	# execute_command('%s/samfilter.py 2d %s %s' % (GOLDEN_BUNDLE_SRC, sam_file, sam_file_2d)); 			# This calculates only the consensus.
	# process_sam(sam_file_2d, reference_path, vcf_mutations_path, 0);


	### Process directly a given VCF file:
	# vcf_file = 'data/variant/ecoliR7.3/consensus-graphmap-params_mex_mutref_5-4-8-2-nowrongcig-cov_20.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/consensus-graphmap-params_mex_mutref_2-nowrongcig-evalue-cov_20.variant.vcf';
	vcf_file = 'data/variant/ecoliR7.3/consensus-graphmap-params_mex_mutref_2-nowrongcig-cov_20.variant.vcf';

	# vcf_file = 'data/variant/ecoliR7.3/LAST-q1-cov_20.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/LAST-q1-uniquebest-cov_20.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/BWAMEM-real_nanopore-cov_20.variant.vcf';
	vcf_file = 'data/variant/ecoliR7.3/BWAMEM-real_nanopore-uniquebest-cov_20.variant.vcf';

	# vcf_file = 'data/variant/ecoliR7.3/LAST-q1-uniquebest-cov_0.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/LAST-q1-cov_20.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/BWAMEM-real_nanopore-leftaligned-cov_20.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/BWAMEM-real_nanopore-uniquebest-cov_20.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/BWAMEM-real_nanopore-uniquebest-cov_0.variant.vcf';
	# vcf_file = 'data/variant/ecoliR7.3/BWAMEM-real_nanopore-cov_20.variant.vcf';

	### Simulated reads on a mutated reference:
	# vcf_file = 'data/variant/consensus-bwamem-params_real_nanopore-cov_20.variant.vcf';
	# vcf_file = 'data/variant/consensus-graphmap.sam-cov_20.variant.vcf';
	# vcf_file = 'data/variant/consensus-graphmap.sam.sam-evalue-cov_20.variant.vcf';	

	# process_vcf(vcf_file, vcf_mutations_path);



	### Test for extracting kmers from FP SNPs:
	# process_vcf('CYP2D6/graphmap-20150419_mex_CYP2D6-CYP2D6-2d-0.99.vcf', 'CYP2D6/correct_snps.vcf');
	# pos_list = vcf_extract_positions('CYP2D6/graphmap-20150419_mex_CYP2D6-CYP2D6-2d-0.99-fp.vcf');
	# get_kmers_from_positions('data/reference/hg19-chr22-GRCh37.p13-NC_000022.10.fa', pos_list, 4);



	##### Simulated reads on a mutated reference #####
	### Retrieve a VCF from Genome, and verify the mutations:
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-myers-sim_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-myers-sim_reads_on_mut_ref-evalue-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-seqan-sim_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-seqan-sim_reads_on_mut_ref-evalue-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-BWAMEM-sim_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-LAST-sim_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-BWAMEM-sim_reads_on_mut_ref-uniquebest-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-LAST-sim_reads_on_mut_ref-uniquebest-cov_20.variant.vcf'), vcf_mutations_path);



	##### Real reads on mutated reference #####
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-myers-real_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-myers-real_reads_on_mut_ref-evalue-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-seqan-real_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-seqan-real_reads_on_mut_ref-evalue-cov_20.variant.vcf'), vcf_mutations_path);

	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-seqan-1111-real_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-graphmap-seqan-1111-real_reads_on_mut_ref-evalue-cov_20.variant.vcf'), vcf_mutations_path);
	

	####################################
	####################################
	####################################
	### Using all E. Coli K-12 reads for testing!
	#####
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-2d-evalue-cov_20.variant.vcf'), vcf_mutations_path);	### Ovo vec daje solidne rezultate!!! tpcount, fpcount, subrecs, trurecs: 4575 124 4699 6401, precision, recall, F1 score: 0.973611406682,0.714732073114,0.824324324324

		# alignments_file: data/alignment//GraphMap-all_ecoli-2d-evalue.sam
		# mpileup_file: data/alignment/GraphMap-all_ecoli-2d-evalue.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4699
		# insertion_count: 3
		# deletion_count: 237
		# num_undercovered_bases: 0
		# num_called_bases: 4639678
		# num_correct_bases: 4634976
		# average_coverage: 42.49
		# 
		# Executing command: "bamsurgeon/etc/evaluator.py -v data/variant/ecoliR7.3/consensus-GraphMap-all_ecoli-2d-evalue-cov_20.variant.vcf.gz -t data/reference/rev_mutants.vcf.gz -m SNV --fn data/variant/ecoliR7.3/consensus-GraphMap-all_ecoli-2d-evalue-cov_20.variant-fn.vcf --fp data/variant/ecoliR7.3/consensus-GraphMap-all_ecoli-2d-evalue-cov_20.variant-fp.vcf"
		# tpcount, fpcount, subrecs, trurecs:
		# 4575 124 4699 6401
		# precision, recall, F1 score: 0.973611406682,0.714732073114,0.824324324324
	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-2d-evalue-1e-100-cov_20.variant.vcf'), vcf_mutations_path);	### Ovo vec daje solidne rezultate!!! tpcount, fpcount, subrecs, trurecs: 4575 124 4699 6401, precision, recall, F1 score: 0.973611406682,0.714732073114,0.824324324324

		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_ecoli-2d-evalue-1e-100.sam
		# mpileup_file: data/alignment/GraphMap-all_ecoli-2d-evalue-1e-100.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4796
		# insertion_count: 3
		# deletion_count: 266
		# num_undercovered_bases: 651
		# num_called_bases: 4639027
		# num_correct_bases: 4634228
		# average_coverage: 41.63

		# tpcount, fpcount, subrecs, trurecs:
		# 4650 146 4796 6401
		# precision, recall, F1 score: 0.969557964971,0.726448992345,0.830579619541
	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-1e-100-cov_20.variant.vcf'), vcf_mutations_path);	### Ovo vec daje solidne rezultate!!! tpcount, fpcount, subrecs, trurecs: 4575 124 4699 6401, precision, recall, F1 score: 0.973611406682,0.714732073114,0.824324324324

		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_ecoli-evalue-1e-100.sam
		# mpileup_file: data/alignment/GraphMap-all_ecoli-evalue-1e-100.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 2789
		# insertion_count: 0
		# deletion_count: 46
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4636886
		# average_coverage: 136.35

		# tpcount, fpcount, subrecs, trurecs:
		# 2768 21 2789 6401
		# precision, recall, F1 score: 0.992470419505,0.432432432432,0.60239390642

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-1e-6-cov_20.variant.vcf'), vcf_mutations_path);

		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_ecoli-evalue-1e-6.sam
		# mpileup_file: data/alignment/GraphMap-all_ecoli-evalue-1e-6.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 2438
		# insertion_count: 0
		# deletion_count: 38
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4637237
		# average_coverage: 145.37
		# tpcount, fpcount, subrecs, trurecs:
		# 2425 13 2438 6401
		# precision, recall, F1 score: 0.994667760459,0.378847055148,0.548704604593

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-cov_20.variant.vcf'), vcf_mutations_path);
		# tpcount, fpcount, subrecs, trurecs:
		# 2373 12 2385 6401
		# precision, recall, F1 score: 0.994968553459,0.370723324481,0.540177555201

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-cov_20.variant.vcf'), vcf_mutations_path);
		# tpcount, fpcount, subrecs, trurecs:
		# 685 0 685 6401
		# precision, recall, F1 score: 1.0,0.10701452898,0.193338978267

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-0-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_ecoli-evalue-0.sam
		# mpileup_file: data/alignment/GraphMap-all_ecoli-evalue-0.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 3191
		# insertion_count: 0
		# deletion_count: 62
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4636484
		# average_coverage: 116.77
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 3160 31 3191 6401
		# precision, recall, F1 score: 0.99028517706,0.493672863615,0.658882402002

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-2d-evalue-1e-300-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: GraphMap-all_ecoli-2d-evalue-1e-300.sam
		# mpileup_file: ./GraphMap-all_ecoli-2d-evalue-1e-300.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4963
		# insertion_count: 3
		# deletion_count: 321
		# num_undercovered_bases: 2680
		# num_called_bases: 4636998
		# num_correct_bases: 4632032
		# average_coverage: 39.09

		# tpcount, fpcount, subrecs, trurecs:
		# 4793 170 4963 6401
		# precision, recall, F1 score: 0.96574652428,0.748789251679,0.843541006688

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5643
		# insertion_count: 0
		# deletion_count: 34
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4634032
		# average_coverage: 42.14

		# tpcount, fpcount, subrecs, trurecs:
		# 5313 330 5643 6401
		# precision, recall, F1 score: 0.941520467836,0.83002655835,0.88226502823

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4986
		# insertion_count: 0
		# deletion_count: 15
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4634689
		# average_coverage: 47.22
		#
		# tpcount, fpcount, subrecs, trurecs:
		# 4802 184 4986 6401
		# precision, recall, F1 score: 0.963096670678,0.750195281987,0.84341793273

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-0-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-0.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-0.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5804
		# insertion_count: 0
		# deletion_count: 37
		# num_undercovered_bases: 2268
		# num_called_bases: 4637407
		# num_correct_bases: 4631603
		# average_coverage: 39.62
		#
		# tpcount, fpcount, subrecs, trurecs:
		# 5397 407 5804 6401
		# precision, recall, F1 score: 0.929875947622,0.843149507889,0.884391642769

	####################################
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-300.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-300.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5794
		# insertion_count: 0
		# deletion_count: 37
		# num_undercovered_bases: 2268
		# num_called_bases: 4637407
		# num_correct_bases: 4631613
		# average_coverage: 39.82
		#
		# tpcount, fpcount, subrecs, trurecs:
		# 5391 403 5794 6401
		# precision, recall, F1 score: 0.930445288229,0.842212154351,0.884132841328

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-LAST-all_2d_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//LAST-all_2d_reads_on_mut_ref.sam
		# mpileup_file: data/alignment/LAST-all_2d_reads_on_mut_ref.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 6974
		# insertion_count: 0
		# deletion_count: 25
		# num_undercovered_bases: 421
		# num_called_bases: 4639254
		# num_correct_bases: 4632280
		# average_coverage: 50.90

		# tpcount, fpcount, subrecs, trurecs:
		# 5279 1695 6974 6401
		# precision, recall, F1 score: 0.756954402065,0.824714888299,0.78938317757

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-LAST-all_2d_reads_on_mut_ref-uniquebest-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//LAST-all_2d_reads_on_mut_ref-uniquebest.sam
		# mpileup_file: data/alignment/LAST-all_2d_reads_on_mut_ref-uniquebest.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5698
		# insertion_count: 0
		# deletion_count: 36
		# num_undercovered_bases: 543
		# num_called_bases: 4639130
		# num_correct_bases: 4633432
		# average_coverage: 41.10

		# Executing command: "bamsurgeon/etc/evaluator.py -v data/variant/ecoliR7.3/consensus-LAST-all_2d_reads_on_mut_ref-uniquebest-cov_20.variant.vcf.gz -t data/reference/rev_mutants.vcf.gz -m SNV --fn data/variant/ecoliR7.3/consensus-LAST-all_2d_reads_on_mut_ref-uniquebest-cov_20.variant-fn.vcf --fp data/variant/ecoliR7.3/consensus-LAST-all_2d_reads_on_mut_ref-uniquebest-cov_20.variant-fp.vcf"
		# tpcount, fpcount, subrecs, trurecs:
		# 5404 294 5698 6401
		# precision, recall, F1 score: 0.948402948403,0.844243087018,0.893296966691

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-BWAMEM-all_2d_reads_on_mut_ref-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//BWAMEM-all_2d_reads_on_mut_ref.sam
		# mpileup_file: data/alignment/BWAMEM-all_2d_reads_on_mut_ref.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 6009
		# insertion_count: 0
		# deletion_count: 1356
		# num_undercovered_bases: 1140
		# num_called_bases: 4638535
		# num_correct_bases: 4632526
		# average_coverage: 40.74
		#
		# tpcount, fpcount, subrecs, trurecs:
		# 5423 586 6009 6401
		# precision, recall, F1 score: 0.902479613912,0.847211373223,0.87397260274

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-BWAMEM-all_2d_reads_on_mut_ref-uniquebest-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//BWAMEM-all_2d_reads_on_mut_ref-uniquebest.sam
		# mpileup_file: data/alignment/BWAMEM-all_2d_reads_on_mut_ref-uniquebest.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 6027
		# insertion_count: 0
		# deletion_count: 1414
		# num_undercovered_bases: 1847
		# num_called_bases: 4637826
		# num_correct_bases: 4631799
		# average_coverage: 40.05
		#
		# tpcount, fpcount, subrecs, trurecs:
		# 5428 599 6027 6401
		# precision, recall, F1 score: 0.900613904098,0.847992501172,0.873511425813

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-1111-evalue-1e-100-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-1111-evalue-1e-100.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-1111-evalue-1e-100.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4712
		# insertion_count: 2
		# deletion_count: 333
		# num_undercovered_bases: 618
		# num_called_bases: 4639059
		# num_correct_bases: 4634345
		# average_coverage: 41.66
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 4565 147 4712 6401
		# precision, recall, F1 score: 0.968803056027,0.713169817216,0.821560334743



	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-BWAMEM-all_ecoli-uniquebest-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//BWAMEM-all_ecoli-uniquebest.sam
		# mpileup_file: data/alignment/BWAMEM-all_ecoli-uniquebest.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4670
		# insertion_count: 0
		# deletion_count: 40
		# num_undercovered_bases: 8
		# num_called_bases: 4639667
		# num_correct_bases: 4634997
		# average_coverage: 123.77
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 4549 121 4670 6401
		# precision, recall, F1 score: 0.97408993576,0.71067020778,0.821786649806	

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-LAST-all_ecoli-uniquebest-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//LAST-all_ecoli-uniquebest.sam
		# mpileup_file: data/alignment/LAST-all_ecoli-uniquebest.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 4360
		# insertion_count: 0
		# deletion_count: 1
		# num_undercovered_bases: 8
		# num_called_bases: 4639667
		# num_correct_bases: 4635307
		# average_coverage: 141.52
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 4311 49 4360 6401
		# precision, recall, F1 score: 0.98876146789,0.673488517419,0.801226651798

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-cov_20.variant.vcf'), vcf_mutations_path);
		# tpcount, fpcount, subrecs, trurecs:
		# 685 0 685 6401
		# precision, recall, F1 score: 1.0,0.10701452898,0.193338978267

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-cov_20.variant.vcf'), vcf_mutations_path);
		# tpcount, fpcount, subrecs, trurecs:
		# 2373 12 2385 6401
		# precision, recall, F1 score: 0.994968553459,0.370723324481,0.540177555201

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-0-cov_20.variant.vcf'), vcf_mutations_path);
		# tpcount, fpcount, subrecs, trurecs:
		# 3160 31 3191 6401
		# precision, recall, F1 score: 0.99028517706,0.493672863615,0.658882402002

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_ecoli-evalue-0-mapq30-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_ecoli-evalue-0-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_ecoli-evalue-0-mapq30.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 3272
		# insertion_count: 0
		# deletion_count: 66
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4636403
		# average_coverage: 109.20
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 3242 30 3272 6401
		# precision, recall, F1 score: 0.990831295844,0.506483361975,0.67031944588		

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-0-mapq30-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-0-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-0-mapq30.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5830
		# insertion_count: 0
		# deletion_count: 41
		# num_undercovered_bases: 3712
		# num_called_bases: 4635963
		# num_correct_bases: 4630133
		# average_coverage: 39.07
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5407 423 5830 6401
		# precision, recall, F1 score: 0.927444253859,0.844711763787,0.884146839997

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-mapq30-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-mapq30.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5318
		# insertion_count: 0
		# deletion_count: 23
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4634357
		# average_coverage: 44.38
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5069 249 5318 6401
		# precision, recall, F1 score: 0.953177886423,0.791907514451,0.865090878061

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1-mapq30-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1-mapq30.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5421
		# insertion_count: 0
		# deletion_count: 25
		# num_undercovered_bases: 117
		# num_called_bases: 4639558
		# num_correct_bases: 4634137
		# average_coverage: 43.62
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5143 278 5421 6401
		# precision, recall, F1 score: 0.948717948718,0.803468208092,0.870072745728

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq30-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq30.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5666
		# insertion_count: 0
		# deletion_count: 36
		# num_undercovered_bases: 1517
		# num_called_bases: 4638158
		# num_correct_bases: 4632492
		# average_coverage: 41.50
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5323 343 5666 6401
		# precision, recall, F1 score: 0.93946346629,0.831588814248,0.882240822077

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq30-cov_0.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq30.bam.mpileup
		# coverage_threshold: 0
		# snp_count: 5674
		# insertion_count: 0
		# deletion_count: 36
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4634001
		# average_coverage: 41.51
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5323 351 5674 6401
		# precision, recall, F1 score: 0.938138879098,0.831588814248,0.8816563147

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-10-mapq30-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-10-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-10-mapq30.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5487
		# insertion_count: 0
		# deletion_count: 25
		# num_undercovered_bases: 146
		# num_called_bases: 4639529
		# num_correct_bases: 4634042
		# average_coverage: 43.10
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5194 293 5487 6401
		# precision, recall, F1 score: 0.946601057044,0.81143571317,0.873822341857

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-10-mapq30-cov_0.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-10-mapq30.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-10-mapq30.bam.mpileup
		# coverage_threshold: 0
		# snp_count: 5488
		# insertion_count: 0
		# deletion_count: 25
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4634187
		# average_coverage: 43.11

		# tpcount, fpcount, subrecs, trurecs:
		# 5194 294 5488 6401
		# precision, recall, F1 score: 0.946428571429,0.81143571317,0.873748843469

	####################################
	# process_vcf(copy_vcf('graphmap/simulated-ref/data/alignment/analysis-intermediate/consensus-GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq3-cov_20.variant.vcf'), vcf_mutations_path);
		# [Consensus statistics]
		# alignments_file: data/alignment//GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq3.sam
		# mpileup_file: data/alignment/GraphMap-all_2d_reads_on_mut_ref-seqan-evalue-1e-100-mapq3.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 5648
		# insertion_count: 0
		# deletion_count: 34
		# num_undercovered_bases: 0
		# num_called_bases: 4639675
		# num_correct_bases: 4634027
		# average_coverage: 42.10
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5315 333 5648 6401
		# precision, recall, F1 score: 0.941041076487,0.83033900953,0.88223089053



	# process_vcf('data/alignment/all_2d/GraphMap-consvar_20150502_anchor-evalue-1e-100-cov_20.variant.vcf', vcf_mutations_path);
	# process_vcf('data/alignment/all_ecoli/GraphMap-consvar_20150502_anchor-evalue-1e-100-cov_20.variant.vcf', vcf_mutations_path);

	# process_vcf('data/alignment/all_2d/GraphMap-consvar_20150502_anchorsq-evalue-1e-100-cov_20.variant.vcf', vcf_mutations_path);
		# [4639600] snps = 6172, insertions = 0, deletions = 52, undercovered = 25245, coverage = 33.18alignments_file: data/alignment/all_2d/GraphMap-consvar_20150502_anchorsq-evalue-1e-100.sam
		# mpileup_file: data/alignment/all_2d/GraphMap-consvar_20150502_anchorsq-evalue-1e-100.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 6172
		# insertion_count: 0
		# deletion_count: 52
		# num_undercovered_bases: 25320
		# num_called_bases: 4614355
		# num_correct_bases: 4608183
		# average_coverage: 33.18
		# 
		# tpcount, fpcount, subrecs, trurecs:
		# 5505 667 6172 6401
		# precision, recall, F1 score: 0.891931302657,0.860021871583,0.875685993796

	# process_vcf('data/alignment/all_2d/GraphMap-consvar_20150502_anchorsq-cov_20.variant.vcf', vcf_mutations_path);

	# process_vcf('data/alignment/all_2d/GraphMap-consvar_20150503_anchorsq-evalue-1e-100-cov_20.variant.vcf', vcf_mutations_path);
		# Processing file "data/alignment/all_2d/GraphMap-consvar_20150503_anchorsq-evalue-1e-100.sam"...
		# Coverage threshold: 20
		# [4639600] snps = 6172, insertions = 0, deletions = 52, undercovered = 25245, coverage = 33.18alignments_file: data/alignment/all_2d/GraphMap-consvar_20150503_anchorsq-evalue-1e-100.sam
		# mpileup_file: data/alignment/all_2d/GraphMap-consvar_20150503_anchorsq-evalue-1e-100.bam.mpileup
		# coverage_threshold: 20
		# snp_count: 6172
		# insertion_count: 0
		# deletion_count: 52
		# num_undercovered_bases: 25320
		# num_called_bases: 4614355
		# num_correct_bases: 4608183
		# 
		# average_coverage: 33.18
		# tpcount, fpcount, subrecs, trurecs:
		# 5505 667 6172 6401
		# precision, recall, F1 score: 0.891931302657,0.860021871583,0.875685993796	

	# process_vcf('data/alignment/ecoliR7.3/GraphMap-consvar_20150503_anchor-evalue-1e-100-cov_20.variant.vcf', vcf_mutations_path);

	# process_vcf('/home/isovic/work/eclipse-workspace/data/island-test/data/alignment/GraphMap-all_2d-pacbioasm-anchor-cov_20.variant.vcf', 'data/reference/rev_mutants-ecoliK12_polished_assembly.vcf.gz');
		# tpcount, fpcount, subrecs, trurecs:
		# 4545 132 4677 6506
		# precision, recall, F1 score: 0.971776779987,0.698585920689,0.812840919252

	# process_vcf('/home/isovic/work/eclipse-workspace/data/island-test/data/alignment/GraphMap-all_2d-pacbioasm-anchor-evalue-0-cov_20.variant.vcf', 'data/reference/rev_mutants-ecoliK12_polished_assembly.vcf.gz');
		# 4781 183 4964 6506
		# precision, recall, F1 score: 0.963134568896,0.734860129112,0.833653007847	

	# process_vcf('/home/isovic/work/eclipse-workspace/data/island-test/data/alignment/GraphMap-all_2d-pacbioasm-anchor-evalue-0-mapq4-cov_0.variant.vcf', 'data/reference/rev_mutants-ecoliK12_polished_assembly.vcf.gz');
		# tpcount, fpcount, subrecs, trurecs:
		# 4810 1575 6385 6506
		# precision, recall, F1 score: 0.753328112764,0.739317553028,0.746257078582
	


if __name__ == "__main__":
	main();









################################
### Results
################################

"""
BWA-MEM sa svim alignmentima, coverage 20 threshold: precision, recall, F1 score: 0.753882117215,0.705358537728,0.728813559322
BWA-MEM ako samo uzmem jedan alignment za jedan read (isfiltriram duplikate), coverage 20 threshold: precision, recall, F1 score: 0.735437683046,0.706139665677,0.720490954013

LAST sa svim alignmentima, coverage 20 threshold: precision, recall, F1 score: 0.798553144129,0.672551163881,0.730156037992
LAST sa samo jednim alignmentom za jedan read, coverage 20 threshold: precision, recall, F1 score: 0.843985325352,0.682862052804,0.754922279793
LAST sa samo jednim alignmentom za jedan read, coverage 0 threshold: precision, recall, F1 score: 0.838722020279,0.684892985471,0.754041967664
"""



"""
##### Real reads on mutated reference results #####

GraphMap with Myers:
	tpcount, fpcount, subrecs, trurecs:
	1000 22 1022 6401
	precision, recall, F1 score: 0.978473581213,0.156225589752,0.269432843864
GraphMap with Myers and E-value filtering:
	tpcount, fpcount, subrecs, trurecs:
	2841 297 3138 6401
	precision, recall, F1 score: 0.905353728489,0.443836900484,0.595659922424

GraphMap with SeqAn and 5-4-8-6 parameters:
	tpcount, fpcount, subrecs, trurecs:
	1820 82 1902 6401
	precision, recall, F1 score: 0.956887486856,0.284330573348,0.438395760568
GraphMap with SeqAn and E-value filtering:
	tpcount, fpcount, subrecs, trurecs:
	3031 305 3336 6401
	precision, recall, F1 score: 0.908573141487,0.473519762537,0.622573687994

BWA-MEM:
	tpcount, fpcount, subrecs, trurecs:
	4515 1474 5989 6401
	precision, recall, F1 score: 0.753882117215,0.705358537728,0.728813559322
BWA-MEM uniquebest:
	tpcount, fpcount, subrecs, trurecs:
	4520 1626 6146 6401
	precision, recall, F1 score: 0.735437683046,0.706139665677,0.720490954013

LAST:
	tpcount, fpcount, subrecs, trurecs:
	4305 1086 5391 6401
	precision, recall, F1 score: 0.798553144129,0.672551163881,0.730156037992
LAST uniquebest:
	tpcount, fpcount, subrecs, trurecs:
	4371 808 5179 6401
	precision, recall, F1 score: 0.843985325352,0.682862052804,0.754922279793

"""


"""
##### Simulated reads on mutated reference results #####

GraphMap with Myers:
	tpcount, fpcount, subrecs, trurecs:
	4640 1 4641 6401
	precision, recall, F1 score: 0.999784529196,0.724886736447,0.840427458794
GraphMap with Myers and E-value filtering:
	tpcount, fpcount, subrecs, trurecs:
	4650 1 4651 6401
	precision, recall, F1 score: 0.999784992475,0.726448992345,0.841476655809

GraphMap with SeqAn and 5-4-8-6 parameters:
	tpcount, fpcount, subrecs, trurecs:
	5925 3 5928 6401
	precision, recall, F1 score: 0.999493927126,0.925636619278,0.961148511639
GraphMap with SeqAn and E-value filtering:
	tpcount, fpcount, subrecs, trurecs:
	5930 2 5932 6401
	precision, recall, F1 score: 0.999662845583,0.926417747227,0.961647612098

BWA-MEM:
	tpcount, fpcount, subrecs, trurecs:
	5700 10 5710 6401
	precision, recall, F1 score: 0.998248686515,0.890485861584,0.941293039386
BWA-MEM only one alignment per read (uniquebest):
	tpcount, fpcount, subrecs, trurecs:
	5686 9 5695 6401
	precision, recall, F1 score: 0.998419666374,0.888298703328,0.940145502646

LAST:
	tpcount, fpcount, subrecs, trurecs:
	5869 113 5982 6401
	precision, recall, F1 score: 0.981109996657,0.916887986252,0.947912460632
LAST uniquebest:
	tpcount, fpcount, subrecs, trurecs:
	5897 2 5899 6401
	precision, recall, F1 score: 0.999660959485,0.921262302765,0.958861788618

"""
