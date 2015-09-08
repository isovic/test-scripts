#! /bin/sh

### Description of the VCF output: the REF field is the original reference, while the ALT field
### is the induced variant in the mutated reference. Be careful to switch the two before comparing
### the variant calling results to the mutated reference.
REFERENCE=data/references/escherichia_coli.fa
mutatrix/mutatrix -n 1 -m 0 -M 0 -i 0 -X 0 $REFERENCE > data/references/mutants.vcf

### The mutated reference is stored in the path from where the script was run, and the filename is something like:
### 1:fasta_header_from_original_reference

#### This reverses only SNPs. Does not work for indels.
cat ../data-in/mutated-reference/mutants.vcf | awk -F '\t' 'BEGIN {OFS = FS} {if ($0 == /^#.*/) print ; else {a=$4; $4=$5; $5=a; print } }' > ../data-in/mutated-reference/rev_mutants.vcf
bgzip -c ../data-in/mutated-reference/rev_mutants.vcf > ../data-in/mutated-reference/rev_mutants.vcf.gz
tabix -p vcf ../data-in/mutated-reference/rev_mutants.vcf.gz



# Slightly more complicatet reversing of mutations for the case where there are also indel mutations.
#
# python /home/wilma/local/src/src.git/vcf_reverse.py mutatrix.vcf.gz |\
#    sed -e 's,del,INS,' -e 's,ins,del,' -e 's,INS,ins,'
# Afterwards you will still need to add the vcf header (#CHROM...).
# Otherwise the evaluator (see below) will choke on it.
# 
# For evaluation I used the script from bamsurgeon:
# /mnt/pnsg10_home/wilma/local/src/bamsurgeon.git/etc/evaluator.py










### Description of the VCF output: the REF field is the original reference, while the ALT field
### is the induced variant in the mutated reference. Be careful to switch the two before comparing
### the variant calling results to the mutated reference.
# mutatrix/mutatrix -n 1 -m 0 -M 0 -i 0 -X 0 data/reference/ecoliK12_polished_assembly.fasta > data/reference/mutants-ecoliK12_polished_assembly.vcf

### The mutated reference is stored in the path from where the script was run, and the filename is something like:
### 1:fasta_header_from_original_reference

#### This reverses only SNPs. Does not work for indels.
# cat data/reference/mutants-ecoliK12_polished_assembly.vcf | awk -F '\t' 'BEGIN {OFS = FS} {if ($0 == /^#.*/) print ; else {a=$4; $4=$5; $5=a; print } }' > data/reference/rev_mutants-ecoliK12_polished_assembly.vcf
# bgzip -c data/reference/rev_mutants-ecoliK12_polished_assembly.vcf > data/reference/rev_mutants-ecoliK12_polished_assembly.vcf.gz
# tabix -p vcf data/reference/rev_mutants-ecoliK12_polished_assembly.vcf.gz
