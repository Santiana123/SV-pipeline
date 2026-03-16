#!/bin/bash
#PBS -j oe 
#PBS -q queue_name
#PBS -V 
#PBS -l nodes=1:ppn=40

cd $PBS_O_WORKDIR
source ~/.bashrc
conda activate SV

threads=40
name=
ref=$PBS_O_WORKDIR/
query=$PBS_O_WORKDIR/
reads=$PBS_O_WORKDIR/
read_group="@RG\tID:${name}\tSM:${name}\tLB:${name}-lib\tPL:PACBIO"

###pbmm2-cutesv2 pipline
mkdir 1.pbmm2-cutesv2
cd 1.pbmm2-cutesv2
pbmm2 index ${ref} ${name}.mmi
pbmm2 align --sort -j $threads ${name}.mmi ${reads} ${name}.bam
samtools view -@ $threads -b -q 20 -F 1284 ${name}.bam -o ${name}.f.bam && \
samtools index -@ $threads ${name}.f.bam
rm ${name}.bam
mkdir work-place
cuteSV ${name}.f.bam ${ref} ${name}.vcf ./work-place -t $threads -l 50 -L 100000 --genotype \
--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
--max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
rm -r work-place
rm ${name}.f.bam
cd ../
###minimap2-sniffles2 pipline
mkdir 2.minimap2-sniffles2
cd 2.minimap2-sniffles2
minimap2 -t $threads -ax map-pb ${ref} ${reads} > ${name}.sam
samtools sort -@ $threads -o ${name}.bam ${name}.sam
samtools view -@ $threads -b -q 20 -F 1284 ${name}.bam -o ${name}.f.bam && \
samtools index -@ $threads ${name}.f.bam
rm ${name}.sam ${name}.bam
sniffles --input ${name}.f.bam --vcf ${name}.vcf --reference ${ref} -t 1
rm ${name}.f.bam
cd ../
###ngmlr-svim pipline
mkdir 3.ngmlr-svim
cd 3.ngmlr-svim
ngmlr -r ${ref} -q ${reads} -o ${name}.sam -t ${threads} -x pacbio --bam-fix -i 0.85 -R 0.5
sambamba view -S -f bam -t $threads -o ${name}.bam ${name}.sam
sambamba sort -t $threads -o ${name}.s.bam ${name}.bam
samtools view -@ $threads -b -q 20 -F 1284 ${name}.s.bam -o ${name}.f.bam
samtools index -@ $threads ${name}.f.bam
rm ${name}.sam ${name}.bam
mkdir work-place
svim alignment --min_mapq 20 --min_sv_size 50 ./work-place ${name}.f.bam ${ref}
mv work-place/variants.vcf ./
rm ${name}.f.bam
cd ../
###winnowmap-pbsv pipline
mkdir 4.winnowmap-pbsv
cd 4.winnowmap-pbsv
winnowmap -x map-pb -a -Y -R $read_group -t $threads $ref $reads -o ${name}.sam
samtools sort -@ $threads -o ${name}.bam ${name}.sam
samtools view -@ $threads -b -q 20 -F 1284 ${name}.bam -o ${name}.f.bam
samtools index -@ $threads ${name}.f.bam
rm ${name}.sam ${name}.bam
pbsv discover ${name}.f.bam ${name}.svsig.gz
pbsv call -m 50 $ref ${name}.svsig.gz ${name}.vcf
rm ${name}.f.bam
cd ../
###miniamp2-SVIM-asm pipline
mkdir 5.minimap2-SVIM-asm
cd 5.minimap2-SVIM-asm
minimap2 -t $threads -ax map-pb ${ref} ${query} > ${name}.sam
samtools sort -@ $threads -o ${name}.bam ${name}.sam
samtools view -@ $threads -b -q 20 -F 1284 ${name}.bam -o ${name}.f.bam && \
samtools index -@ $threads ${name}.f.bam
rm ${name}.sam ${name}.bam
mkdir work-place
svim-asm haploid ./work-place ${name}.f.bam ${ref} --min_sv_size 50
mv work-place/variants.vcf ./
rm ${name}.f.bam
cd ../
###nucmer-SyRI pipline
mkdir 6.nucmer-SyRI
cd 6.nucmer-SyRI
nucmer -c 500 -l 100 --maxmatch ${ref} ${query} -p ${name}-SV
delta-filter -m -i 90 -l 100 ./${name}-SV.delta > out.filtered.delta
show-coords -THrd ./out.filtered.delta >./out.filtered.coords
syri -c ./out.filtered.coords -d ./out.filtered.delta -r $ref -q $query
awk -F'\t' '{if (length($4) > 50 || length($5) > 50) print $0;}' syri.vcf > ${name}.vcf
cd ../
###SURVIVOR merge
mkdir 7.merge
cd 7.merge
cp ../1.pbmm2-cutesv2/${name}.vcf ./cutesv.vcf
cp ../2.minimap2-sniffles2/${name}.vcf ./sniffles2.vcf
cp ../3.ngmlr-svim/variants.vcf ./svim.vcf
cp ../4.winnowmap-pbsv/${name}.vcf ./pbsv.vcf
cp ../5.minimap2-SVIM-asm/variants.vcf ./svim-asm.vcf
cp ../6.nucmer-SyRI/${name}.vcf ./SyRI.vcf
ls *.vcf > sample
SURVIVOR merge sample 1000 3 1 1 0 50 merge.vcf
