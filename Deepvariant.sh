#!/bin/bash
#PBS -j oe 
#PBS -q FAFU3
#PBS -V 
#PBS -l nodes=1:ppn=32

cd $PBS_O_WORKDIR
threads=32
ref=ZB-XYh.N+Yh.fa
fq=reads.fq.gz
name=ZB_vs_SVP11

source /public/home/yuejingjing/.bashrc
source activate /public/home/yuejingjing/mambaforge/envs/test

cp $ref ./reference.fasta
samtools faidx reference.fasta

minimap2 -ax map-hifi -t ${threads} -R "@RG\tID:${name}\tSM:${name}\tPL:pacbio" reference.fasta ${fq} | samtools view -@ ${threads} -b - > ${name}.bam
samtools view -@ ${threads} -b -q 30 ${name}.bam | \
samtools sort -@ ${threads} -o ${name}.f.s.bam && \
samtools index -@ ${threads} ${name}.f.s.bam
rm ${name}.bam

mkdir -p split
mkdir -p vcf
mkdir -p log_$name

cut -f1 reference.fasta.fai | grep -v '^$' | while read -r contig; do
mkdir -p log_$name/"${contig}"
done
perl /public/home/yuejingjing/tyh/script/splitChrByFa.pl reference.fasta
for fa in fabyChr/*.fasta; do samtools faidx "$fa"; done

cut -f1 reference.fasta.fai | parallel -j4 "samtools view -@8 -b ${name}.f.s.bam {} > split/${name}.{}.bam"
cut -f1 reference.fasta.fai | parallel -j4 "samtools index -@8 split/${name}.{}.bam"

cut -f1 reference.fasta.fai | parallel -j4 "
      contig={};
      singularity run deepvariant_1.6.1.sif \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref=fabyChr/{}.fasta \
        --reads=split/${name}.{}.bam \
        --output_vcf=vcf/${name}.{}.vcf.gz \
        --intermediate_results_dir=log_${name}/{} \
        --num_shards=8 \
        --sample_name=${name} \
         1>split/${name}.{}.log \
        2>split/${name}.{}.err
    "

cut -f1 reference.fasta.fai | parallel "zcat vcf/${name}.{}.vcf.gz | grep -E '^#|1/1|0/1|1/0' > vcf/${name}.{}.filtered.vcf"
cut -f1 reference.fasta.fai | parallel "bgzip -c vcf/${name}.{}.filtered.vcf > vcf/${name}.{}.filtered.vcf.gz && tabix -p vcf vcf/${name}.{}.filtered.vcf.gz"

ls vcf/${name}.*.filtered.vcf.gz | sort > vcf.list
bcftools concat -Oz -o ${name}.ccs.concat.vcf.gz -f vcf.list
#rm -r split log_$name ${name}.f.s.bam ${name}.f.s.bam.bai
