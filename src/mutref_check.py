#! /usr/bin/python

import os;
import sys;
import subprocess;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
sys.path.append(SCRIPT_PATH + '/../tools/samscripts/src/');

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
	execute_command('%s/../tools/bamsurgeon/etc/evaluator.py -v %s.gz -t %s -m SNV --fn %s --fp %s' % (SCRIPT_PATH, vcf_file, vcf_mutations_path, fn_file, fp_file));


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

def main():
	# sam_file = 'data/alignment/ecoliR7.3/graphmap-params_final_20150417_mutref.sam';
	# sam_file = 'data/alignment/ecoliR7.3/LAST-q1.sam';
	# sam_file = 'data/alignment/ecoliR7.3/BWAMEM-real_nanopore.sam';

	reference_path = 'data/mutated-reference/mutated_ecoli.fa';
	vcf_mutations_path = 'data/mutated-reference/rev_mutants.vcf.gz';

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-GraphMap-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-GraphMap-anchor-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-BLASR-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-BWAMEM-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-DALIGNER-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-LAST-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-marginAlign-mutecoli_ecoliR7.3-graphmap_anchor-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-marginAlign-mutecoli_ecoliR7.3-graphmap-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-marginAlign-mutecoli_ecoliR7.3-last-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);

	vcf_file = 'data/out/mutecoli_ecoliR7.3_1/analysis-intermediate/consensus-GraphMap-anchor_new-mutecoli_ecoliR7.3-cov_20.variant.vcf';
	process_vcf(vcf_file, vcf_mutations_path);



if __name__ == "__main__":
	main();
