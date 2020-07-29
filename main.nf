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
params.dbsnp = Channel.fromPath("$params.refDir/*dbSNP142_human_GRCh38.snps.vcf.gz").getVal()
params.dbsnp_idx = Channel.fromPath("$params.refDir/*dbSNP142_human_GRCh38.snps.vcf.gz.tbi").getVal()
params.indels = Channel.fromPath("$params.refDir/Mills_and_1000G_gold_standard.*.vcf.gz").getVal()
params.indels_idx = Channel.fromPath("$params.refDir/Mills_and_1000G_gold_standard.*.vcf.gz.tbi").getVal()

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
        
        input:
        tuple val(base), file(reads) from reads_ch
        tuple file(fasta), file(fai) from Channel.value([params.fasta, params.fai])
        tuple file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])
        
        output:
        tuple val(base), file("${base}.bam") into bamMapped

        script:
        readGroup = "@RG\\tID:HT52VDMXX\\tPU:HT52VDMXX:1\\tSM:METIN\\tLB:METIN\\tPL:illumina"
        """
        bwa mem -K 100000000 -R \"${readGroup}\" -t 8 -M $fasta $reads | \
        samtools sort --threads 8 - > ${base}.bam
        """
}


process MarkDuplicates{

	input:
	tuple val(base), file(bam) from bamMapped

	output:
	tuple val(base), file("${base}.md.bam"), file("${base}.md.bam.bai") into bam_duplicates_marked

	script:
	"""
	gatk --java-options -Xmx8g \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 50000 \
        --INPUT $bam \
        --METRICS_FILE ${base}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${base}.md.bam
    
        mv ${base}.md.bai ${base}.md.bam.bai
	"""
}


process BQSR{

	input:
	tuple val(base), file(bam), file(bai) from bam_duplicates_marked
	tuple file(fasta), file(fai), file(dict), file(intervals) from Channel.value([params.fasta, params.fai, params.dict, params.targets])
	tuple file(dbsnp), file(dbsnp_idx), file(indels), file(indels_idx) from Channel.value([params.dbsnp, params.dbsnp_idx, params.indels, params.indels_idx])

	output:
	tuple val(base), file("${base}.recal.bam"), file("${base}.recal.bam.bai") into BQSR_bams

	script:
	"""
	gatk --java-options -Xmx8g \
	BaseRecalibrator \
	-I $bam \
	-O ${base}.recal.table \
	--tmp-dir . \
	-R $fasta \
	-L $intervals \
	--known-sites $dbsnp \
	--known-sites $indels

	gatk --java-options -Xmx8g \
	ApplyBQSR \
	-I $bam
	-O ${base}.recal.bam \
	-R $fasta \
	-L $intervals \
	--bqsr-recal-file ${base}.recal.table

	samtools index ${base}.recal.bam ${base}.recal.bam.bai
	"""
}

