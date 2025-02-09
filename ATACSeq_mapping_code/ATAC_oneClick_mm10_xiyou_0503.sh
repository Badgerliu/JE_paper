#!/bin/bash

# Author Huan Liu
# Author Organization Wuhan University School of Stomatology
# Usage: Working with ATAC-seq raw fq.gz files and resulted in bam, bw and peak files for whole peaks and NFRs, EXCLUDING MNR, DNR etc.
# Date: 2022-10-02

now=`date +'%Y-%m-%d %H:%M:%S'`
start_time=$(date --date="$now" +%s);
 




# Part 0 Activate Conda Environment
source deactivate 
conda deactivate
echo "list all the conda env avail"
conda env list
conda activate ATAC
echo "list the currently activated conda environment"
conda env list

# Part 1 Trimming with Trimmomatic (v0.39).
echo "First, we will work with trimming using Trimmomatic0.39"
mkdir ./trimmed
mkdir ./paired
mkdir ./unpaired

for i in $(ls *.fq.gz | rev |cut -c 10- |rev |uniq)
	do
		java -jar /home/liuhuan/apps/Trimmomatic/trimmomatic-0.39.jar PE -phred33 ${i}_R1.fq.gz ${i}_R2.fq.gz ${i}_R1_paired.fq.gz ${i}_R1_unpaired.fq.gz ${i}_R2_paired.fq.gz ${i}_R2_unpaired.fq.gz ILLUMINACLIP:/home/liuhuan/apps/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:7
		mv ${i}_R*_paired* ./paired
		mv ${i}_R*_unpaired* ./unpaired
		echo done with ${i}
	done

mv *.fq.gz ./trimmed
cd ./paired
mv *_paired.fq.gz ..
cd ..

# Part 2 Mapping with Bowtie2, using mm10 ref genome.
echo "Now, we will start mapping with bowtie2."
mkdir ./mapped
for i in $(ls *_paired.fq.gz|rev | cut -c 17- |rev| uniq)
	do
		bowtie2 -x /home/liuhuan/refergenome/mm10_scATAC/fasta/genome -1 ${i}_R1_paired.fq.gz -2 ${i}_R2_paired.fq.gz -X 1500 -p 16 -S ${i}.sam
		sed '/chrM/d;/random/d;/chrUn/d' ${i}.sam > ${i}_m.sam
		samtools view -bT /home/liuhuan/refergenome/mm10_scATAC/fasta/genome.fa ${i}_m.sam > ${i}_m.bam
		rm *.sam
		samtools sort -o ${i}.bam ${i}_m.bam
		picard MarkDuplicates I=${i}.bam O=${i}_m_d.bam M=dups.txt REMOVE_DUPLICATES=true
		samtools index ${i}_m_d.bam && rm ${i}.bam && rm ${i}_m.bam
		mv ${i}_m_d.bam* ./mapped
		echo OK done with ${i}
	done


mv *_paired.fq.gz ./paired

cd ./mapped # so here you will be in ./mapped folder

# Part 3 QC After mapping
echo "Now move to deeptools virtual environment"

chmod -R 755 /home/liuhuan/apps/ATAC_scripts/ #permission to all ATAC-scripts

echo "QC for insert distribution"
/home/liuhuan/apps/ATAC_scripts/insert_xiyou.sh
mkdir ./insert_size_distribution
mv *_insert_size_* ./insert_size_distribution # Pack all size distribution results where you could analyze your library quality.
echo "Done with insert size distribution. Make sure you will see the right nucleosome distribution pattern. If not, repeat your experiment."

# Part 4 BigWig files generation and peak calling for whole peaks of ATAC-seq
echo "Do bw generation for ATAC-seq peaks"
/home/liuhuan/apps/ATAC_scripts/generate_bw_mm10_xiyou.sh
echo "OK done with bw files for ATAC-seq peaks"
echo "do peak calling for ATAC-seq peaks"
/home/liuhuan/apps/ATAC_scripts/shift_and_call_peak_mm10_xiyou.sh
echo "OK done with peak calling for ATAC-seq peaks"


# Part 5 Working with NFR, including isolation of NFR reads, generation of bw files and peak calling.
echo "Start with NFR identification"
for j in $(ls *_m_d.bam| rev |cut -c 9- |rev |uniq)
	do
		python /home/liuhuan/apps/ATAC_scripts/size_selection_xiyou.py ${j}_m_d.bam
		echo done with size_selection for ${j}
	done

mkdir ./NFR
rm -rf  monoNucReads.* # We don't want MNR.

echo "Let's move all the nucFreeReads.* files to NFR folder for SB conversion"
mv nucFreeReads.* ./NFR # Now, you will be in ./NFR folder

cd ./NFR
/home/liuhuan/apps/ATAC_scripts/SB_conversion_mm10_xiyou.sh
rm -rf *.sam


echo "Do bw generation for NFRs"
/home/liuhuan/apps/ATAC_scripts/generate_bw_mm10_xiyou.sh
echo "OK done with bw files for NFRs"
echo "do peak calling for NFRs"
/home/liuhuan/apps/ATAC_scripts/shift_and_call_peak_mm10_xiyou.sh
echo "OK done with peak calling for NFRs"

 
now=`date +'%Y-%m-%d %H:%M:%S'`
end_time=$(date --date="$now" +%s);
echo "used time:"$((end_time-start_time))"s"