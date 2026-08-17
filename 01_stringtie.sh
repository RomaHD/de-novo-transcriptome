##run de novo assembly for each bam file
cd /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/bam_files
for file in `ls *.bam`
do /omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/stringtie/stringtie  ${file} --rf -o /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/stringtie_output/${file}.gtf -p 5 \
-m 200 -f 0.05 -c 1.5 -p 20
echo $file
done


##merge de novo assembled trancripts (both cell lines)
cd /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/all_gtfs/
/omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/stringtie/stringtie *.gtf --merge -G ../../de_novo_rna2/reference_data/gencode.v29lift37.annotation.gtf \
-p 20 -o mergedTranscripts.gtf 

# move mergedTranscripts.gtf to /stringtie_output

##gffcompare 
cd /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/
mkdir gffcompare_output
../de_novo_rna2/gffcompare/gffcompare -R -r ../de_novo_rna2/reference_data/gencode.v29lift37.annotation.gtf \
-o gffcompare_output ./stringtie_output/mergedTranscripts.gtf   

cd gffcompare_output
module load Perl/5.38.2
perl ../../de_novo_rna2/gffcompare_output/gff3sort.pl gffcompare_output.annotated.gtf  > gffCompare.annotated.sorted.gtf

## next step: manually change chromosome names (then count abundances)

## estimate transcript abundances and generate read coverage tables

module load StringTie/2.2.3-GCC-14.1.0
# COLO800
cd /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/bam_files
#mkdir  ../stringtie_output2
for file in *.bam; do
	filename="${file%.*}"
	mkdir  ../stringtie_output2/$filename
	StringTie -p 20 -G ../gffcompare_output/gffCompare.annotated.sorted_ver1.gtf -e -o ../work/tmp.gtf -A ../stringtie_output2/${file}.tsv -b ../stringtie_output2/${filename} ${file}
	echo $file
done

#A375
cd /omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/bam_files
for file in *.bam; do
	filename="${file%.*}"
	mkdir  ../stringtie_output2/$filename
	../stringtie/stringtie --rf -p 20 -G ../../de_novo_colo800/gffcompare_output/gffCompare.annotated.sorted_ver1.gtf -e -o ../work/tmp.gtf -A ../stringtie_output2/${file}.tsv -b ../stringtie_output2/${filename} ${file}
	echo $file
done

