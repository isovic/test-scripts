#! /usr/bin/python

import os;
import sys;
import subprocess;
import multiprocessing;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
tools_path = '%s/../tools' % (SCRIPT_PATH);
SAMSCRIPTS = '%s/../tools/samscripts/src/';

sys.path.append('%s/../tools/samscripts/src/' % (SCRIPT_PATH));
import fastqparser;



def main():
	setup_tools();
	# RUN_CONSENSUS_TEST();
	# RUN_MUTATED_REFERENCE_TEST();
	# RUN_SV_TEST();
	RUN_AMPLICON_TEST();

	# RUN_MUTATED_REFERENCE_ADDITIONAL_TESTS();

def RUN_CONSENSUS_TEST():
	run_all_mappers_only(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), ('%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % SCRIPT_PATH), 'ecoliR7.3', '%s/../data/out/fig3cd/' % (SCRIPT_PATH), 'nanopore');
	evaluate_alignments(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), ('%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % SCRIPT_PATH), 'ecoliR7.3', '%s/../data/out/fig3cd/' % (SCRIPT_PATH));

def RUN_MUTATED_REFERENCE_TEST():
	run_all_mappers_only(('%s/../data/mutated-reference/mutated_ecoli.fa' % SCRIPT_PATH), ('%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % SCRIPT_PATH), 'mutecoli_ecoliR7.3', '%s/../data/out/mutecoli_ecoliR7.3_1/' % (SCRIPT_PATH), 'nanopore');
	execute_command('cp -R %s/../data/out/mutecoli_ecoliR7.3_1 %s/../data/out/mutecoli_ecoliR7.3_on_original_ref' % (SCRIPT_PATH, SCRIPT_PATH))
	### Evaluate mappings on the original reference.
	evaluate_alignments(('%s/../data/mutated-reference/mutated_ecoli.fa' % SCRIPT_PATH), ('%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % SCRIPT_PATH), 'mutecoli_ecoliR7.3', '%s/../data/out/mutecoli_ecoliR7.3_1/' % (SCRIPT_PATH));
	### Evaluate mappings on the mutated reference on the original non-mutated reference.
	evaluate_alignments(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), ('%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % SCRIPT_PATH), 'mutecoli_ecoliR7.3', '%s/../data/out/mutecoli_ecoliR7.3_on_original_ref' % (SCRIPT_PATH));

def RUN_MUTATED_REFERENCE_ADDITIONAL_TESTS():
	# generate_mutated_reference(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), 0.0001, 'temp/mutated-refs/');
	pass;

def RUN_SV_TEST():
	### First run the mappers on the original reference, to detect the differences that normaly exist and need to be omitted from further comparisons.
	run_all_mappers_only(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), ('%s/../data/reads-all_2d_for_sv/all_2d_for_sv.fastq' % SCRIPT_PATH), 'all_2d_for_sv', '%s/../data/out/all_2d_for_sv-normal_ref/' % (SCRIPT_PATH), 'nanopore');
	### Second, run the mappers on the modified reference.
	run_all_mappers_only(('%s/../data/reference_for_sv/escherichia_coli-indel_events.fa' % SCRIPT_PATH), ('%s/../data/reads-all_2d_for_sv/all_2d_for_sv.fastq' % SCRIPT_PATH), 'all_2d_for_sv', '%s/../data/out/all_2d_for_sv-indel_ref/' % (SCRIPT_PATH), 'nanopore');

def RUN_AMPLICON_TEST():
	reference = '%s/../data/amplicons-f1000/reference/ref_chr6_chr22-hg19_v38.fa' % (SCRIPT_PATH);
	reads = '%s/../data/amplicons-f1000/reads/reads_2d-f1000.fastq' % (SCRIPT_PATH);
	dataset_name = 'f1000amplicons_2d';
	out_path = '%s/../data/out/f1000amplicons_2d/' % (SCRIPT_PATH);
	### Map all the amplicon reads to the chr6 and chr22 references.
	# run_all_mappers_only(reference, reads, 'nanopore', out_path, dataset_name, do_not_recalc=True, is_circular=False);

	REGION_CYP2D6 = ['gi|224589814|ref|NC_000022.10|:42522077-42527144', 'CYP2D6'];
	REGION_HLAA = ['gi|224589818|ref|NC_000006.11|:29909854-29913805', 'HLA-A'];
	REGION_HLAB = ['gi|224589818|ref|NC_000006.11|:31321279-31325303', 'HLA-B'];

	regions = [REGION_CYP2D6, REGION_HLAA, REGION_HLAB];
	sam_out_folder = '%s/inregion/';
	sam_path = '%s/marginAlign-nanopore-nospecialchars-with_AS.sam' % (out_path);
	for region in regions:
		region_name = region[1];
		sys.stderr.write('Running region %s:' % (region[1]));

		region_to_use = region;
		if ('marginalign' in os.path.basename(sam_path).lower()):
			region_to_use = [re.sub('[^0-9a-zA-Z]', '_', region[0]), region[1]];

		### First prepare the alignments for variant calling. This includes filtering the uniquely aligned reads, taking only 2d reads, and taking only reads that fully span the region.
		filter_spanning_reads(True, region_to_use, reads, sam_path, sam_out_folder, reference_path=None, leftalign=False);




def execute_command(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_get_stdout(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	[out, err] = p.communicate()
	sys.stderr.write('\n');
	return [out, err];

def execute_command_w_dryrun(dry_run, command):
	sys.stderr.write('[Executing] "%s"\n' % command);
	if (dry_run == False):
		subprocess.call(command, shell=True);

def measure_command_wrapper(out_filename):
#	if (USE_BASICDEFINES_ == True):
#		return basicdefines.measure_command(out_filename);
#	else:
	return '/usr/bin/time --format "Command line: %%C\\nReal time: %%e s\\nCPU time: -1.0 s\\nUser time: %%U s\\nSystem time: %%S s\\nMaximum RSS: %%M kB\\nExit status: %%x" --quiet -o %s ' % out_filename;

def setup_tools():
	if (not os.path.exists('%s/../tools' % (SCRIPT_PATH))):
		execute_command('mkdir %s/../tools' % (SCRIPT_PATH));
	
	# if (not os.path.exists('%s/../tools/marginAlign' % (SCRIPT_PATH))):
	# 	execute_command('cd %s/../tools; git clone https://github.com/isovic/marginAlign.git; cd marginAlign; git submodule update --init; make' % (SCRIPT_PATH));

	# if (not os.path.exists('%s/../tools/wrapper-dev' % (SCRIPT_PATH))):
	# 	execute_command('cd %s/../tools; git clone https://github.com/isovic/wrapper-dev.git' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/samscripts' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/samscripts.git' % (SCRIPT_PATH));
	
	if (not os.path.exists('%s/../tools/aligneval' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/aligneval.git; cd aligneval; ./setup.py aligners; ./setup.py tools' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/mutatrix/' % (SCRIPT_PATH))):
		execute_command('cd %s/../packs; tar -xvf mutatrix.tar.gz' % (SCRIPT_PATH));
		execute_command('mv packs/mutatrix tools/' % (SCRIPT_PATH));



def evaluate_alignments(reference, reads, dataset_name, out_path):
	out_collect_file = '%s/collected.csv' % (out_path);

	out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file hcalc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

def collect_alignments(reference, reads, dataset_name, out_path):
	out_collect_file = '%s/collected.csv' % (out_path);

	out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file hcollect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

	out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
	execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

### Params:
### 	orig_reference is the reference to which to map to.
###		orig_reads is the file containing the reads to map to.
###		dataset_name is the name that will be added to the suffix of the output SAM file.
###		out_path is the path to the folder
def run_all_mappers_only(orig_reference, orig_reads, dataset_name, out_path, machine_name, do_not_recalc=True, is_circular=True):
	if (not os.path.exists(out_path)):
		sys.stderr.write('Creating output path: "%s".\n' % out_path);
		execute_command('mkdir -p %s' % out_path);
	num_threads = multiprocessing.cpu_count() / 2;
	reads_basename = os.path.basename(os.path.splitext(orig_reads)[0]);
	out_collect_file = '%s/collected.csv' % (out_path);

	# dataset_name = os.path.splitext(os.path.basename(orig_reads))[0];
	sys.stderr.write('Dataset name: "%s".\n' % (dataset_name));

	out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_daligner.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
	if ((not os.path.exists(out_sam))):
		execute_command('%s/aligneval/wrappers/wrapper_graphmap.py run %s %s nanopore%s %s %s' % (tools_path, orig_reads, orig_reference, ('circ' if (is_circular == True) else ''), out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
	if (not (os.path.exists(out_sam))):
		execute_command('%s/aligneval/wrappers/wrapper_graphmap.py run %s %s anchor%s %s anchor-%s' % (tools_path, orig_reads, orig_reference, ('circ' if (is_circular == True) else ''), out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_lastal.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_bwamem.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_blasr.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));



	out_sam = '%s/marginAlign-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_marginalign.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/marginAlignGraphMap-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_marginaligngraphmap.py run %s %s %s %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/marginAlignGraphMap-anchor-%s.sam' % (out_path, dataset_name);
	if (not os.path.exists(out_sam)):
		execute_command('%s/aligneval/wrappers/wrapper_marginaligngraphmap.py run %s %s anchor %s %s' % (tools_path, orig_reads, orig_reference, out_path, 'anchor-' + dataset_name));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));






### Mutates the given reference, and writes the mutations in a vcf file.
def generate_mutated_reference(reference_path, snp_rate, out_path):
	reference_path = os.path.abspath(reference_path);
	out_path = os.path.abspath(out_path);
	if (not os.path.exists(out_path)):
		os.makedirs(out_path);

	out_prefix = '%s/mutated_%s_%f' % (out_path, os.path.splitext(os.path.basename(reference_path))[0], snp_rate);
	out_vcf = os.path.abspath('%s.vcf' % (out_prefix));
	out_rev_vcf = '%s/rev_%s.vcf' % (out_path, os.path.basename(out_prefix));
	ref_ext = os.path.splitext(reference_path)[-1];
	out_ref_file = '%s%s' % (out_prefix, ref_ext);

	sys.stderr.write('Mutating the reference using Mutatrix, output VCF file: "%s".\n' % (out_vcf));
	execute_command('cd %s; %s/mutatrix/mutatrix --snp-rate %f --population-size 1 --microsat-min-len 0 --mnp-ratio 0 --indel-rate 0 --indel-max 0 %s > %s' % (out_path, tools_path, snp_rate, reference_path, out_vcf));

	sys.stderr.write('Reversing the SNP bases in the VCF file, output VCF file: "%s".\n' % (out_rev_vcf));
	execute_command(r"cat %s | awk -F '\t' 'BEGIN {OFS = FS} {if ($0 == /^#.*/) print ; else {a=$4; $4=$5; $5=a; print } }' > %s" % (out_vcf, out_rev_vcf));

	sys.stderr.write('Compressing and indexing the VCF file.\n');
	execute_command('bgzip -c %s > %s.gz' % (out_rev_vcf, out_rev_vcf));
	execute_command('tabix -p vcf %s.gz' % (out_rev_vcf));

	### Mutatrix splits all reference sequences into separate files. This part of code joins them back into one FASTA file.
	[headers, lengths] = fastqparser.get_headers_and_lengths(reference_path);
	print headers;
	all_files = ['%s/1:%s:0%s' % (out_path, header.split(' ')[0], ref_ext) for header in headers];
	if (os.path.exists(out_ref_file)):
		os.rename(out_ref_file, '%s.bak' % (out_ref_file));
	for ref_file in all_files:
		### Take care of the special characters.
		escaped_ref_file = ref_file.replace('|', '\|');
		execute_command('cat %s >> %s' % (escaped_ref_file, out_ref_file));
		if (len(ref_file) > 0 and ('*' in ref_file) == False):
			print 'Removing file: "%s".' % (ref_file);
			os.remove(ref_file);


















#######################
### These are additional functions for the amplicon test (filtering and such stuff).
#######################
def filter_spanning_reads(dry_run, region, reads_path, sam_in_path, sam_out_folder, reference_path=None, leftalign=False):
	if not os.path.exists(sam_out_folder):
		sys.stderr.write('Creating folder "%s".\n' % (sam_out_folder));
		os.makedirs(sam_out_folder);

	# if (('last' in os.path.basename(sam_in_path).lower()) or ('marginalign' in os.path.basename(sam_in_path).lower())):
	# 	sys.stderr.write('[0] Adding quality values to alignments...\n');
	# 	execute_command_w_dryrun(dry_run, '%s/samfilter.py addqv %s %s %s-addedqv.sam' % (SAMSCRIPTS, sam_in_path, reads_path, os.path.splitext(sam_in_path)[0]));
	# 	sys.stderr.write('\n');
	# 	sam_in_path = '%s-addedqv.sam' % (os.path.splitext(sam_in_path)[0]);

#	# temp_sam_region = '%s/%s-%s.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]), region[1]);
#	# temp_region_uniquebest = '%s/%s-uniquebest.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_sam_region)[0]));
#	# temp_region_uniquebest_evalue = '%s/%s-evalue-1.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_uniquebest_region)[0]));
#	# temp_before_2d = temp_uniquebest_region;

	temp_sam_uniquebest = '%s/%s-uniquebest.sam' % (sam_out_folder, os.path.basename(os.path.splitext(sam_in_path)[0]));
	temp_uniquebest_region = '%s/%s-%s.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_sam_uniquebest)[0]), region[1]);
	temp_region_uniquebest_evalue = '%s/%s-evalue-1.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_uniquebest_region)[0]));
	temp_before_2d = temp_uniquebest_region;

	### Extract only unique reads.
	sys.stderr.write('[1] Extracting only the unique reads with the best score "%s"...\n' % (temp_sam_uniquebest));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py uniquebest %s %s' % (SAMSCRIPTS, sam_in_path, temp_sam_uniquebest));
	sys.stderr.write('\n');

	### Extract only the reads which fully cover the region:
	sys.stderr.write('[2] Extracting only reads that fully span region "%s" to file "%s"...\n' % (region[0], temp_uniquebest_region));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py regionfull "%s" %s %s' % (SAMSCRIPTS, region[0], temp_sam_uniquebest, temp_uniquebest_region));
	# execute_command_w_dryrun('%s/samfilter.py regionpartial "%s" %s %s' % (SAMSCRIPTS, region[0], temp_sam_uniquebest, temp_uniquebest_region));
	sys.stderr.write('\n');

	if ('graphmap' in os.path.basename(sam_in_path).lower()):
		sys.stderr.write('[2.1] Filtering by E-value (GraphMap specific) to file "%s"...\n' % (temp_region_uniquebest_evalue));
		execute_command_w_dryrun(dry_run, '%s/samfilter.py evalue 1 %s %s' % (SAMSCRIPTS, temp_uniquebest_region, temp_region_uniquebest_evalue));
		sys.stderr.write('\n');
		temp_before_2d = temp_region_uniquebest_evalue;

	temp_1d = '%s/%s-1d.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_before_2d)[0]));
	sys.stderr.write('[3] Extracting only the 1d reads to file "%s"...\n' % (temp_1d));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py 1d %s %s' % (SAMSCRIPTS, temp_before_2d, temp_1d));
	sys.stderr.write('\n');

	temp_2d = '%s/%s-2d.sam' % (sam_out_folder, os.path.basename(os.path.splitext(temp_before_2d)[0]));
	sys.stderr.write('[4] Extracting only the 2d reads to file "%s"...\n' % (temp_2d));
	execute_command_w_dryrun(dry_run, '%s/samfilter.py 2d %s %s' % (SAMSCRIPTS, temp_before_2d, temp_2d));
	sys.stderr.write('\n');

	### Convert the SAM file to BAM format:
	sys.stderr.write('[5.1] Converting all aligned reads in the region to BAM format from file "%s"...' % (temp_before_2d));
	execute_command_w_dryrun(dry_run, '%s/convert_to_bam.sh %s' % (SAMSCRIPTS, os.path.splitext(temp_before_2d)[0]));
	bam_all = '%s-sorted.bam' % os.path.splitext(temp_before_2d)[0];
	sys.stderr.write('\n');
	sys.stderr.write('[5.2] Converting the 1d aligned reads to BAM format from file "%s"...' % (temp_1d));
	execute_command_w_dryrun(dry_run, '%s/convert_to_bam.sh %s' % (SAMSCRIPTS, os.path.splitext(temp_1d)[0]));
	bam_1d = '%s-sorted.bam' % os.path.splitext(temp_1d)[0];
	sys.stderr.write('\n');
	sys.stderr.write('[5.3] Converting the 2d aligned reads to BAM format from file "%s"...' % (temp_2d));
	execute_command_w_dryrun(dry_run, '%s/convert_to_bam.sh %s' % (SAMSCRIPTS, os.path.splitext(temp_2d)[0]));
	bam_2d = '%s-sorted.bam' % os.path.splitext(temp_2d)[0];
	sys.stderr.write('\n');

	bam_all_reads_in_region = bam_all;
	bam_1d_reads_in_region = bam_1d
	bam_2d_reads_in_region = bam_2d

	if (leftalign == True):
		if (reference_path == None):
			sys.stderr.write('ERROR: Reference not specified, cannot leftalign the BAM file.\n');
		else:
			sys.stderr.write('[6.1] Left aligning all aligned reads in the region from file "%s"...' % (bam_all));
			bam_all_leftaligned = '%s-leftaligned.bam' % (os.path.splitext(bam_all)[0]);
			execute_command_w_dryrun(dry_run, 'cat %s | %s/bamleftalign -f %s > %s' % (bam_all, TOOLS_PATH, reference_path, bam_all_leftaligned));
			execute_command_w_dryrun(dry_run, 'samtools index %s' % (bam_all_leftaligned));
			sys.stderr.write('\n');

			sys.stderr.write('[6.2] Left aligning all aligned reads in the region from file "%s"...' % (bam_1d));
			bam_1d_leftaligned = '%s-leftaligned.bam' % (os.path.splitext(bam_1d)[0]);
			execute_command_w_dryrun(dry_run, 'cat %s | %s/bamleftalign -f %s > %s' % (bam_1d, TOOLS_PATH, reference_path, bam_1d_leftaligned));
			execute_command_w_dryrun(dry_run, 'samtools index %s' % (bam_1d_leftaligned));
			sys.stderr.write('\n');

			sys.stderr.write('[6.3] Left aligning all aligned reads in the region from file "%s"...' % (bam_2d));
			bam_2d_leftaligned = '%s-leftaligned.bam' % (os.path.splitext(bam_2d)[0]);
			execute_command_w_dryrun(dry_run, 'cat %s | %s/bamleftalign -f %s > %s' % (bam_2d, TOOLS_PATH, reference_path, bam_2d_leftaligned));
			execute_command_w_dryrun(dry_run, 'samtools index %s' % (bam_2d_leftaligned));
			sys.stderr.write('\n');

			bam_all_reads_in_region = bam_all_leftaligned;
			bam_1d_reads_in_region = bam_1d_leftaligned
			bam_2d_reads_in_region = bam_2d_leftaligned

	return [bam_all_reads_in_region, bam_1d_reads_in_region, bam_2d_reads_in_region];



# cat ../data-in/mutated-reference/mutants.vcf | awk -F '\t' 'BEGIN {OFS = FS} {if ($0 == /^#.*/) print ; else {a=$4; $4=$5; $5=a; print } }' > ../data-in/mutated-reference/rev_mutants.vcf
# bgzip -c ../data-in/mutated-reference/rev_mutants.vcf > ../data-in/mutated-reference/rev_mutants.vcf.gz
# tabix -p vcf ../data-in/mutated-reference/rev_mutants.vcf.gz



# 	### Run marginAlign with GraphMap.
# 	out_sam = '%s/marginAlign-%s-graphmap.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		memtime_file = '%s.memtime' % (os.path.splitext(out_sam)[0]);
# #		output_model_file = '%s/../data/hmm/hmm-%s-graphmap.txt' % (SCRIPT_PATH, dataset_name);
# 		# output_model_file = '%s/../data/hmm/hmm-ecoliR7.3-graphmap.txt' % (SCRIPT_PATH);
# 		output_model_file = '%s/../data/hmm/hmm-%s-graphmap.txt' % (SCRIPT_PATH, reads_basename);
# 		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
# 		if (os.path.exists(jobtree)):
# 			execute_command('rm -r %s' % (jobtree));
# 		if (not os.path.exists(output_model_file)):
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmap --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		else:
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmap --jobTree %s --inputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

# 	### Run marginAlign default.
# 	out_sam = '%s/marginAlign-%s-last.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		# memtime_file = '%s/marginAlign-%s-last.memtime' % (out_path, dataset_name);
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		memtime_file = '%s.memtime' % (os.path.splitext(out_sam)[0]);
# #		output_model_file = '%s/../data/hmm/hmm-%s-last.txt' % (SCRIPT_PATH, dataset_name);
# 		# output_model_file = '%s/../data/hmm/hmm-ecoliR7.3-last.txt' % (SCRIPT_PATH);
# 		output_model_file = '%s/../data/hmm/hmm-%s-last.txt' % (SCRIPT_PATH, reads_basename);
# 		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
# 		if (os.path.exists(jobtree)):
# 			execute_command('rm -r %s' % (jobtree));
# 		if (not os.path.exists(output_model_file)):
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		else:
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --jobTree %s --inputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

# 	### Run marginAlign with GraphMap.
# 	out_sam = '%s/marginAlign-%s-graphmap_anchor.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		memtime_file = '%s.memtime' % (os.path.splitext(out_sam)[0]);
# #		output_model_file = '%s/../data/hmm/hmm-%s-graphmap_anchor.txt' % (SCRIPT_PATH, dataset_name);
# #		output_model_file = '%s/../data/hmm/hmm-ecoliR7.3-graphmap_anchor.txt' % (SCRIPT_PATH);
# 		output_model_file = '%s/../data/hmm/hmm-%s-graphmap_anchor.txt' % (SCRIPT_PATH, reads_basename);
# 		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
# 		if (os.path.exists(jobtree)):
# 			execute_command('rm -r %s' % (jobtree));
# 		if (not os.path.exists(output_model_file)):
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmapanchor --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		else:
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmapanchor --jobTree %s --inputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
		


# ### Params:
# ### 	orig_reference is the reference to which to map to.
# ###		orig_reads is the file containing the reads to map to.
# ###		dataset_name is the name that will be added to the suffix of the output SAM file.
# ###		out_path is the path to the folder
# def run_all_mappers(orig_reference, orig_reads, dataset_name, out_path, do_not_recalc=True):
# 	if (not os.path.exists(out_path)):
# 		# execute_command('mkdir -p %s/../data/out/fig3cd/' % (SCRIPT_PATH));
# 		execute_command('mkdir -p %s' % out_path);

# # /usr/bin/time marginAlign/marginAlign $READS $REFERENCE data/marginAlign-graphmap-ecoliR7.3-all.sam --graphmap --jobTree ./jobTree-2d-wgraphmap --em --maxThreads=12 --logInfo --defaultMemory=100000000000 --defaultCpu=12
# 	num_threads = multiprocessing.cpu_count() / 2;

# 	# orig_reads = '%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % (SCRIPT_PATH);
# 	# orig_reference = '%s/../data/reference/escherichia_coli.fa' % (SCRIPT_PATH);
# 	# reads = '%s/../data/reads-ecoliR7.3/ecoliR7.3-marginalign.fastq' % (SCRIPT_PATH);
# 	# reference = '%s/../data/reference/escherichia_coli-nospecchar.fa' % (SCRIPT_PATH);
# 	reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 	reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);

# 	reads_basename = os.path.basename(os.path.splitext(orig_reads)[0]);

# 	out_collect_file = '%s/collected.csv' % (out_path);

# 	# dataset_name = os.path.splitext(os.path.basename(orig_reads))[0];
# 	sys.stderr.write('Dataset name: "%s".\n' % (dataset_name));
	
# 	if (not os.path.exists(reads)):
# 		execute_command('%s/../tools/samscripts/src/fastqfilter.py marginalign %s %s' % (SCRIPT_PATH, orig_reads, reads));
# 	if (not os.path.exists(reference)):
# 		execute_command('%s/../tools/samscripts/src/fastqfilter.py specialchars %s %s' % (SCRIPT_PATH, orig_reference, reference));



# 	out_sam = '%s/DALIGNER-%s.sam' % (out_path, dataset_name);
# 	if (not os.path.exists(out_sam)):
# 		reads = orig_reads[0:-1] + 'a';
# 		if (orig_reads[-1] == 'q'):
# 			execute_command('%s/../../../golden-bundle/src/fastq2fasta.py %s %s' % (SCRIPT_PATH, orig_reads, reads));
# 		execute_command('%s/aligneval/wrappers/wrapper_daligner.py run %s %s nanopore %s %s' % (tools_path, reads, orig_reference, os.path.dirname(out_sam), dataset_name));
# 		execute_command('%s/samscripts/src/alignmentstats.py file hcalc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file hcollect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file hcalc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

# 	out_sam = '%s/GraphMap-%s.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		execute_command('%s/aligneval/wrappers/wrapper_graphmap.py run %s %s nanoporecirc %s %s' % (tools_path, orig_reads, orig_reference, os.path.dirname(out_sam), dataset_name));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));

# 	out_sam = '%s/GraphMap-anchor-%s.sam' % (out_path, dataset_name);
# 	if (not (os.path.exists(out_sam))):
# 		execute_command('%s/aligneval/wrappers/wrapper_graphmap.py run %s %s anchorcirc %s anchor-%s' % (tools_path, orig_reads, orig_reference, os.path.dirname(out_sam), dataset_name));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));

# #	out_sam = '%s/GraphMap-anchor_new-%s.sam' % (out_path, dataset_name);
# #	if (not (os.path.exists(out_sam))):
# #		reads = orig_reads;
# #		execute_command('/home/isovic/graphmap/git/graphmap/bin/graphmap-not_release -a anchor -B 0 -b 3 -v 5 -r %s -d %s -o %s' % (orig_reference, orig_reads, out_sam));
# #		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, reads, out_collect_file));
# #	else:
# #		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# ##		execute_command('%s/samscripts/src/alignmentstats.py file hcollect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, reads, out_collect_file));
# #		execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, reads, out_collect_file));
# #		# execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, reads, out_collect_file));

# 	out_sam = '%s/LAST-%s.sam' % (out_path, dataset_name);
# 	if (not os.path.exists(out_sam)):
# 		execute_command('%s/aligneval/wrappers/wrapper_lastal.py run %s %s nanopore %s %s' % (tools_path, orig_reads, orig_reference, os.path.dirname(out_sam), dataset_name));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));

# 	out_sam = '%s/BWAMEM-%s.sam' % (out_path, dataset_name);
# 	if (not os.path.exists(out_sam)):
# 		execute_command('%s/aligneval/wrappers/wrapper_bwamem.py run %s %s nanopore %s %s' % (tools_path, orig_reads, orig_reference, os.path.dirname(out_sam), dataset_name));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));

# 	out_sam = '%s/BLASR-%s.sam' % (out_path, dataset_name);
# 	if (not os.path.exists(out_sam)):
# 		execute_command('%s/aligneval/wrappers/wrapper_blasr.py run %s %s nanopore %s %s' % (tools_path, orig_reads, orig_reference, os.path.dirname(out_sam), dataset_name));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, orig_reference, orig_reads, out_collect_file));



# 	### Run marginAlign with GraphMap.
# 	out_sam = '%s/marginAlign-%s-graphmap.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		memtime_file = '%s.memtime' % (os.path.splitext(out_sam)[0]);
# #		output_model_file = '%s/../data/hmm/hmm-%s-graphmap.txt' % (SCRIPT_PATH, dataset_name);
# 		# output_model_file = '%s/../data/hmm/hmm-ecoliR7.3-graphmap.txt' % (SCRIPT_PATH);
# 		output_model_file = '%s/../data/hmm/hmm-%s-graphmap.txt' % (SCRIPT_PATH, reads_basename);
# 		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
# 		if (os.path.exists(jobtree)):
# 			execute_command('rm -r %s' % (jobtree));
# 		if (not os.path.exists(output_model_file)):
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmap --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		else:
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmap --jobTree %s --inputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

# 	### Run marginAlign default.
# 	out_sam = '%s/marginAlign-%s-last.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		# memtime_file = '%s/marginAlign-%s-last.memtime' % (out_path, dataset_name);
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		memtime_file = '%s.memtime' % (os.path.splitext(out_sam)[0]);
# #		output_model_file = '%s/../data/hmm/hmm-%s-last.txt' % (SCRIPT_PATH, dataset_name);
# 		# output_model_file = '%s/../data/hmm/hmm-ecoliR7.3-last.txt' % (SCRIPT_PATH);
# 		output_model_file = '%s/../data/hmm/hmm-%s-last.txt' % (SCRIPT_PATH, reads_basename);
# 		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
# 		if (os.path.exists(jobtree)):
# 			execute_command('rm -r %s' % (jobtree));
# 		if (not os.path.exists(output_model_file)):
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		else:
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --jobTree %s --inputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));

# 	### Run marginAlign with GraphMap.
# 	out_sam = '%s/marginAlign-%s-graphmap_anchor.sam' % (out_path, dataset_name);
# 	if ((not os.path.exists(out_sam))):
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		memtime_file = '%s.memtime' % (os.path.splitext(out_sam)[0]);
# #		output_model_file = '%s/../data/hmm/hmm-%s-graphmap_anchor.txt' % (SCRIPT_PATH, dataset_name);
# #		output_model_file = '%s/../data/hmm/hmm-ecoliR7.3-graphmap_anchor.txt' % (SCRIPT_PATH);
# 		output_model_file = '%s/../data/hmm/hmm-%s-graphmap_anchor.txt' % (SCRIPT_PATH, reads_basename);
# 		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
# 		if (os.path.exists(jobtree)):
# 			execute_command('rm -r %s' % (jobtree));
# 		if (not os.path.exists(output_model_file)):
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmapanchor --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		else:
# 			execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmapanchor --jobTree %s --inputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
# 		execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 	else:
# 		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));
# 		reads = '%s-marginalign.fastq' % (os.path.splitext(orig_reads)[0]);
# 		reference = '%s-nospecchar.fa' % (os.path.splitext(orig_reference)[0]);
# 		if (do_not_recalc == True):
# 			execute_command('%s/samscripts/src/alignmentstats.py file collect %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
# 		else:
# 			execute_command('%s/samscripts/src/alignmentstats.py file calc %s %s %s 20 >> %s' % (tools_path, out_sam, reference, reads, out_collect_file));
		





if __name__ == "__main__":
	main()
