#! /usr/bin/python

import os;
import sys;
import subprocess;
import multiprocessing;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
tools_path = '%s/../tools' % (SCRIPT_PATH);

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

def RUN_SV_TEST():
	run_all_mappers_only(('%s/../data/reference/escherichia_coli.fa' % SCRIPT_PATH), ('%s/../data/reads-all_2d_for_sv/all_2d_for_sv.fastq' % SCRIPT_PATH), 'all_2d_for_sv', '%s/../data/out/all_2d_for_sv-normal_ref/' % (SCRIPT_PATH), 'nanopore');
	run_all_mappers_only(('%s/../data/reference_for_sv/escherichia_coli-indel_events.fa' % SCRIPT_PATH), ('%s/../data/reads-all_2d_for_sv/all_2d_for_sv.fastq' % SCRIPT_PATH), 'all_2d_for_sv', '%s/../data/out/all_2d_for_sv-indel_ref/' % (SCRIPT_PATH), 'nanopore');

def RUN_AMPLICON_TEST():
	reference = '%s/../data/amplicons-f1000/reference/ref_chr6_chr22-hg19_v38.fa' % (SCRIPT_PATH);
	reads = '%s/../data/amplicons-f1000/reads/reads_2d-f1000.fastq' % (SCRIPT_PATH);
	dataset_name = 'f1000amplicons';
	out_path = '%s/../data/out/f1000amplicons/' % (SCRIPT_PATH);
	run_all_mappers_only(reference, reads, 'nanopore', out_path, dataset_name, do_not_recalc=True, is_circular=False);





def execute_command(command):
	sys.stderr.write('[Executing] %s\n' % (command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_get_stdout(command):
	sys.stderr.write('[Executing] %s\n' % (command));
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	[out, err] = p.communicate()
	sys.stderr.write('\n');

	return [out, err];

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
		execute_command('mv packs/mutatrix tools/mutatrix' % (SCRIPT_PATH));



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
		execute_command('%s/aligneval/wrappers/wrapper_marginaligngraphmap.py run %s %s anchor %s %s' % (tools_path, orig_reads, orig_reference, machine_name, out_path, 'anchor-' + dataset_name));
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
