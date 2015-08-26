#! /usr/bin/python

import os;
import sys;
import subprocess;
import multiprocessing;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

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
	
	if (not os.path.exists('%s/../tools/marginAlign' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/marginAlign.git; cd marginAlign; git submodule update --init; make' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/wrapper-dev' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/wrapper-dev.git' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../tools/samscripts' % (SCRIPT_PATH))):
		execute_command('cd %s/../tools; git clone https://github.com/isovic/samscripts.git' % (SCRIPT_PATH));

	if (not os.path.exists('%s/../data/out/fig3cd/' % (SCRIPT_PATH))):
		execute_command('mkdir -p %s/../data/out/fig3cd/' % (SCRIPT_PATH));

def main():
	setup_tools();
	
	tools_path = '%s/../tools' % (SCRIPT_PATH);

# /usr/bin/time marginAlign/marginAlign $READS $REFERENCE data/marginAlign-graphmap-ecoliR7.3-all.sam --graphmap --jobTree ./jobTree-2d-wgraphmap --em --maxThreads=12 --logInfo --defaultMemory=100000000000 --defaultCpu=12
	num_threads = multiprocessing.cpu_count() / 2;

	orig_reads = '%s/../data/reads-ecoliR7.3/ecoliR7.3.fastq' % (SCRIPT_PATH);
	orig_reference = '%s/../data/reference/escherichia_coli.fa' % (SCRIPT_PATH);
	reads = '%s/../data/reads-ecoliR7.3/ecoliR7.3-nospecchar.fastq' % (SCRIPT_PATH);
	reference = '%s/../data/reference/escherichia_coli-nospecchar.fa' % (SCRIPT_PATH);
	
	if (not os.path.exists(reads)):
		execute_command('%s/../tools/samscripts/src/fastqfilter.py specialchars %s %s' % (SCRIPT_PATH, orig_reads, reads));
	if (not os.path.exists(reference)):
		execute_command('%s/../tools/samscripts/src/fastqfilter.py specialchars %s %s' % (SCRIPT_PATH, orig_reference, reference));
	
	out_sam = '%s/../data/out/fig3cd/marginAlign-ecoliR7.3_all-last.sam' % (SCRIPT_PATH);
	if (not os.path.exists(out_sam)):
		memtime_file = '%s/../data/out/fig3cd/marginAlign-ecoliR7.3_all-last.memtime' % (SCRIPT_PATH);
		output_model_file = '%s/../data/hmm-ecoliR7.3-last.txt' % (SCRIPT_PATH);
		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
		if (os.path.exists(jobtree)):
			execute_command('rm -r %s' % (jobtree));
		execute_command('%s %s/marginAlign/marginAlign %s %s %s --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));

	out_sam = '%s/../data/out/fig3cd/marginAlign-ecoliR7.3_all-graphmap.sam' % (SCRIPT_PATH);
	if (not os.path.exists(out_sam)):
		memtime_file = '%s/../data/out/fig3cd/marginAlign-ecoliR7.3_all-graphmap.memtime' % (SCRIPT_PATH);
		output_model_file = '%s/../data/hmm-ecoliR7.3-graphmap.txt' % (SCRIPT_PATH);
		jobtree = '%s/../jobTree' % (SCRIPT_PATH);
		if (os.path.exists(jobtree)):
			execute_command('rm -r %s' % (jobtree));
		execute_command('%s %s/marginAlign/marginAlign %s %s %s --graphmap --jobTree %s --em --outputModel=%s --maxThreads=%d --logInfo --defaultMemory=100000000000 --defaultCpu=%d' % (measure_command_wrapper(memtime_file), tools_path, reads, reference, out_sam, jobtree, output_model_file, num_threads, num_threads));
	else:
		sys.stderr.write('Warning: File "%s" already exists. Please use another name. Skipping.\n' % (out_sam));



if __name__ == "__main__":
	main()
