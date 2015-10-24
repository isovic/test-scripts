#! /bin/sh

# sam_path=done-all_2d_for_sv-normal_ref
sam_path=all_2d_for_sv-indel_ref
# sam_path=all_2d_for_sv-normal_ref
reference=data/reference/escherichia_coli.fa



samname=GraphMap-anchor-all_2d_for_sv
cp data/out/done-$sam_path/$samname.sam data/out/$sam_path/$samname.sam
cp data/out/done-$sam_path/$samname.memtime data/out/$sam_path/$samname.memtime
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
samname2=${samname}-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=BLASR-all_2d_for_sv
cp data/out/done-$sam_path/$samname.sam data/out/$sam_path/$samname.sam
cp data/out/done-$sam_path/$samname.memtime data/out/$sam_path/$samname.memtime
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
samname2=${samname}-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=LAST-all_2d_for_sv
cp data/out/done-$sam_path/$samname.sam data/out/$sam_path/$samname.sam
cp data/out/done-$sam_path/$samname.memtime data/out/$sam_path/$samname.memtime
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
samname2=${samname}-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=DALIGNER-all_2d_for_sv
cp data/out/done-$sam_path/$samname.sam data/out/$sam_path/$samname.sam
cp data/out/done-$sam_path/$samname.memtime data/out/$sam_path/$samname.memtime
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
samname2=${samname}-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=GraphMap-all_2d_for_sv
cp data/out/done-$sam_path/$samname.sam data/out/$sam_path/$samname.sam
cp data/out/done-$sam_path/$samname.memtime data/out/$sam_path/$samname.memtime
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
samname2=${samname}-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=BWAMEM-all_2d_for_sv
cp data/out/done-$sam_path/$samname.sam data/out/$sam_path/$samname.sam
cp data/out/done-$sam_path/$samname.memtime data/out/$sam_path/$samname.memtime
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/$samname.sam data/out/$sam_path/${samname}-uniquebest.sam
samname2=${samname}-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

exit

samname=marginAlign-all_2d_for_sv
../samscripts/src/samfilter.py generateAS $reference data/out/$sam_path/${samname}.sam data/out/$sam_path/${samname}-with_AS.sam
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/${samname}-with_AS.sam data/out/$sam_path/${samname}-with_AS-uniquebest.sam
samname2=${samname}-with_AS-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=marginAlignGraphMap-all_2d_for_sv
../samscripts/src/samfilter.py generateAS $reference data/out/$sam_path/${samname}.sam data/out/$sam_path/${samname}-with_AS.sam
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/${samname}-with_AS.sam data/out/$sam_path/${samname}-with_AS-uniquebest.sam
samname2=${samname}-with_AS-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

samname=marginAlignGraphMap-anchor-all_2d_for_sv
../samscripts/src/samfilter.py generateAS $reference data/out/$sam_path/${samname}.sam data/out/$sam_path/${samname}-with_AS.sam
../samscripts/src/samfilter.py uniquebest data/out/$sam_path/${samname}-with_AS.sam data/out/$sam_path/${samname}-with_AS-uniquebest.sam
samname2=${samname}-with_AS-uniquebest
../samscripts/src/svconsensus.py data/reference/escherichia_coli.fa data/out/$sam_path/$samname2.sam data/out/$sam_path/analysis-final/$samname2.structvars.csv

exit

