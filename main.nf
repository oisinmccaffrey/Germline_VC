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
@Usage
nextflow -bg -q run BarryDigby/Germline_VC -profile standard, singularity \
--refDir /path/to/genome/files
--------------------------------------------------------------------------------
*/

params.fasta = Channel.fromPath("$params.refDir/*fa").getVal()
params.fai = Channel.fromPath("$params.refDir/*fa.fai").getVal()
params.dict = Channel.fromPath("$params.refDir/*dict").getVal()

params.amb = Channel.fromPath("$params.refDir/*fa.amb").getVal()
params.ann = Channel.fromPath("$params.refDir/*fa.ann").getVal()
params.bwt = Channel.fromPath("$params.refDir/*fa.bwt").getVal()
params.pac = Channel.fromPath("$params.refDir/*fa.pac").getVal()
params.sa = Channel.fromPath("$params.refDir/*fa.sa").getVal()

params.targets = Channel.fromPath("$params.refDir/*exome.targets.bed").getVal()
params.dbsnp = Channel.fromPath("$params.refDir/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz").getVal()
params.dbsnp_idx = Channel.fromPath("$params.refDir/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz.tbi").getVal()
params.indels = Channel.fromPath("$params.refDir/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz").getVal()
params.indels_idx = Channel.fromPath("$params.refDir/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz.tbi").getVal()

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
        tuple file(fasta), file(fai) from Channel.value([params.fasta, params.fai])
        tuple file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])
        
        output:
        tuple val(base), file("${base}.bam") into bamMapped

        script:
        readGroup = "@RG\\tID:HT52VDMXX\\tPU:HT52VDMXX:1\\tSM:METIN\\tLB:METIN\\tPL:illumina"
        """
        bwa mem -K 100000000 -R \"${readGroup}\" -t 4 -M $fasta $reads | \
        samtools sort --threads 4 - > ${base}.bam
        """
}
