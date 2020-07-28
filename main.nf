#!/usr/bin/env nextflow

/*
================================================================================
                           Germline Variant Calling
================================================================================
Galway Genomics
July 2019
--------------------------------------------------------------------------------
@Homepage
https://github.com/BarryDigby/Germline_VC
--------------------------------------------------------------------------------
*/

params.genome = "/data/bdigby/WES/assets/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.dict = "/data/bdigby/WES/assets/GRCh38_full_analysis_set_plus_decoy_hla.dict"
params.fai = "/data/bdigby/WES/assets/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
params.index = "/data/bdigby/WES/assets/index/GRCH38*"
params.targets = "/data/bdigby/WES/assets/20130108.exome.targets.bed"
params.dbsnp = "/data/bdigby/WES/assets/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz"
params.dbsnp_idx = "/data/bdigby/WES/assets/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz.tbi"
params.indels = "/data/bdigby/WES/assets/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz"
params.indels_idx = "/data/bdigby/WES/assets/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz.tbi"
params.reads = "/data/bdigby/WES/reads/*trim_R{1,2}.fastq.gz"
Channel
        .fromFilePairs( params.reads )
        .set{ reads_ch }


/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

process MapReads{
        tag "${base}"

        input:
        tuple val(base), file(reads) from reads_ch
        path(idx) from params.index
        path(genome) from params.genome
        path(fai) from params.fai

        output:
        tuple val(base), file("${base}.bam") into bamMapped

        script:
        readGroup = "@RG\\tID:HT52VDMXX\\tPU:HT52VDMXX:1\\tSM:METIN\\tLB:METIN\\tPL:illumina"
        """
        bwa mem -K 100000000 -R \"${readGroup}\" -t 4 -M $genome $reads \
        samtools sort --threads 4 - > ${base}.bam
        """
}
