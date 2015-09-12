#! /bin/sh

### Find all the variants.
# ../src/consensus.py mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa 20 out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/GraphMap-mutated_ref_draftlike_ecolinmeth out/mutated_ref_draftlike_ecolinmeth/GraphMap-mutated_ref_draftlike_ecolinmeth.sam
# ../src/consensus.py mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa 20 out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/GraphMap-anchor-mutated_ref_draftlike_ecolinmeth out/mutated_ref_draftlike_ecolinmeth/GraphMap-anchor-mutated_ref_draftlike_ecolinmeth.sam
# ../src/consensus.py mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa 20 out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/LAST-mutated_ref_draftlike_ecolinmeth out/mutated_ref_draftlike_ecolinmeth/LAST-mutated_ref_draftlike_ecolinmeth.sam
# ../src/samfilter.py uniquebest out/mutated_ref_draftlike_ecolinmeth/LAST-mutated_ref_draftlike_ecolinmeth.sam out/mutated_ref_draftlike_ecolinmeth/LAST-mutated_ref_draftlike_ecolinmeth-uniquebest.sam
# ../src/consensus.py mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa 20 out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/LAST-mutated_ref_draftlike_ecolinmeth-uniquebest out/mutated_ref_draftlike_ecolinmeth/LAST-mutated_ref_draftlike_ecolinmeth-uniquebest.sam

### GATK needs a dictionary of the reference. To generate this, Picard tools needs to be downloaded, e.g. from: https://github.com/broadinstitute/picard/releases/download/1.138/picard-tools-1.138.zip
java -jar picard-tools-1.138/picard.jar CreateSequenceDictionary R= mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa O= mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.dict
samtools faidx mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa

### Run the creation of the alternate reference:
# java -jar GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa -o output-GraphMap.fasta -V out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/GraphMap-mutated_ref_draftlike_ecolinmeth-cov_20.variant.vcf
java -jar GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa -o output-GraphMap-anchor.fasta -V out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/GraphMap-anchor-mutated_ref_draftlike_ecolinmeth-cov_20.variant.vcf
# java -jar GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa -o output-LAST.fasta -V out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/LAST-mutated_ref_draftlike_ecolinmeth-cov_20.variant.vcf
java -jar GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa -o output-LAST-uniquebest.fasta -V out/mutated_ref_draftlike_ecolinmeth/analysis-intermediate/LAST-mutated_ref_draftlike_ecolinmeth-uniquebest-cov_20.variant.vcf

# dnadiff mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa output-GraphMap.fasta -p comparison-GraphMap
# dnadiff mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.fa output-LAST.fasta -p comparison-LAST
dnadiff reference/escherichia_coli.fa output-GraphMap.fasta -p comparison-GraphMap
dnadiff reference/escherichia_coli.fa output-GraphMap-anchor.fasta -p comparison-GraphMap-anchor
dnadiff reference/escherichia_coli.fa output-LAST.fasta -p comparison-LAST
dnadiff reference/escherichia_coli.fa output-LAST-uniquebest.fasta -p comparison-LAST-uniquebest

echo "Truth SNPs:"
grep ";TYPE=snp" mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.vcf | wc -l
echo "Truth insertions:"
grep ";TYPE=ins" mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.vcf | wc -l
echo "Truth deletions:"
grep ";TYPE=del" mutated-refs/draftlike/mutated_escherichia_coli_snp0.000600_indel0.006700.vcf | wc -l



# Bug za popraviti:
# gi|48994873|gb|U00096.2|	531803	.	C	{}	1000	PASS	DP=28;TYPE=snp
# SNP	pos = 531803	ref = gi|48994873|gb|U00096.2|	coverage = 28	non_indel_cov_curr = 0	most_common_base_count = 0	ref_base = C	cons_base = {}	base_counts = []	insertion_counts = {}	deletion_counts = {}	gi|48994873|gb|U00096.2|	531803	C	28	****************************	*/2+,+0---)/-3-,.//0,,..1+..
