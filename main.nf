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
--refDir /path/to/genome/files --outDir /path/to/publish
--------------------------------------------------------------------------------
*/

/*
 Reference Files
*/

params.fasta = Channel.fromPath("$params.refDir/*fasta").getVal()
params.fai = Channel.fromPath("$params.refDir/*fasta.fai").getVal()
params.dict = Channel.fromPath("$params.refDir/*dict").getVal()

params.amb = Channel.fromPath("$params.refDir/*fasta.amb").getVal()
params.ann = Channel.fromPath("$params.refDir/*fasta.ann").getVal()
params.bwt = Channel.fromPath("$params.refDir/*fasta.bwt").getVal()
params.pac = Channel.fromPath("$params.refDir/*fasta.pac").getVal()
params.sa = Channel.fromPath("$params.refDir/*fasta.sa").getVal()

/*
 Exome Intervals, files
*/

params.intlist = Channel.fromPath("$params.refDir/exome/*.bed.interval_list").getVal()
params.bed = Channel.fromPath("$params.refDir/exome/*.bed").getVal()
params.bedgz = Channel.fromPath("$params.refDir/exome/*.bed.gz").getVal()
params.bedgztbi = Channel.fromPath("$params.refDir/exome/*.bed.gz.tbi").getVal()

/*
 dbSNP, known Indels
*/

params.dbsnp = Channel.fromPath("$params.refDir/dbsnp*.gz").getVal()
params.dbsnptbi = Channel.fromPath("$params.refDir/dbsnp*.tbi").getVal()
params.mills = Channel.fromPath("$params.refDir/Mills_KG*.gz").getVal()
params.millstbi = Channel.fromPath("$params.refDir/Mills_KG*.gz.tbi").getVal()

/*
 Annotation cache, database versions
*/

params.vepcache = Channel.fromPath("/data/VEP/GRCh37/homo_sapiens/99_GRCh37/*").getVal()
params.vepversion = "99"
params.snpeffcache = Channel.fromPath("/data/snpEff/GRCh37/data/GRCh37.87/*").getVal()
params.snpeffversion = "GRCh37.87"

// Not sure where to use these files, omit for now 
//params.omni = Channel.fromPath("$params.refDir/KG_omni*.gz").getVal()
//params.otbi = Channel.fromPath("$params.refDir/KG_omni*.gz.tbi").getVal()
//params.kgp1 = Channel.fromPath("$params.refDir/KG_phase1*.gz").getVal()
//params.ktbi = Channel.fromPath("$params.refDir/KG_phase1*.gz.tbi").getVal()
//params.hpmp = Channel.fromPath("$params.refDir/hapmap*.gz").getVal()
//params.htbi = Channel.fromPath("$params.refDir/hapmap*.gz.tbi").getVal()

//params.gps = Channel.fromPath("$params.refDir/exome/af-only-gnomad.*.vcf.gz").getVal()
//params.gpstbi = Channel.fromPath("$params.refDir/exome/af-only-gnomad.*.vcf.gz.tbi").getVal()

/*
 FASTQ reads
*/

params.reads = "/data/bdigby/WES/reads/*trim_R{1,2}.fastq.gz"
Channel
        .fromFilePairs( params.reads )
        .set{ reads_ch }

/*
 Initialise outDir
*/

params.outDir = ""

/*
================================================================================
                                  PREPROCESSING
================================================================================
*/


process MapReads{
        
	publishDir path: "$params.outDir/analysis/bwa", mode: "copy"
	
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

	publishDir path: "$params.outDir/analysis/mark_dups", mode: "copy"

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

	publishDir path: "$params.outDir/analysis/bqsr", mode: "copy"

	input:
	tuple val(base), file(bam), file(bai) from bam_duplicates_marked
	tuple file(fasta), file(fai), file(dict), file(intlist) from Channel.value([params.fasta, params.fai, params.dict, params.intlist])
	tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
	tuple file(mills), file(millstbi) from Channel.value([params.mills, params.millstbi])

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
	-L $intlist \
	--known-sites $dbsnp \
	--known-sites $mills 

	gatk --java-options -Xmx8g \
	ApplyBQSR \
	-I $bam \
	-O ${base}.recal.bam \
	-R $fasta \
	-L $intlist \
	--bqsr-recal-file ${base}.recal.table

	samtools index ${base}.recal.bam ${base}.recal.bam.bai
	"""
}


/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/


process HaplotypeCaller {

	publishDir path: "$params.outDir/analysis/haplotypecaller", mode: "copy"
	
	input:
	tuple val(base), file(bam), file(bai) from BQSR_bams
	tuple file(fasta), file(fai), file(dict), file(intlist) from Channel.value([params.fasta, params.fai, params.dict, params.intlist])
	tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
	
	output:
	tuple val(base), file("${base}.g.vcf") into gvcfHaplotypeCaller
	
	script:
	"""
	gatk --java-options -Xmx8g \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -L $intlist \
        --D $dbsnp \
        -O ${base}.g.vcf \
        -ERC GVCF
	"""
}


process GenotypeGVCFs {

	publishDir path: "$params.outDir/analysis/genotypeGVCF", mode: "copy"
	
	input:
	tuple val(base), file(gvcf) from gvcfHaplotypeCaller
	tuple file(fasta), file(fai), file(dict), file(intlist) from Channel.value([params.fasta, params.fai, params.dict, params.intlist])
	tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
	
	output:
	tuple val(base), file("${base}.vcf") into vcfGenotypeGVCFs
	
	script:
	"""
	gatk --java-options -Xmx8g \
	IndexFeatureFile \
        -I ${gvcf}
	
	gatk --java-options -Xmx8g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L $intlist \
        --D $dbsnp \
        -V ${gvcf} \
        -O ${base}.vcf
	"""
}


/*
================================================================================
                                 ANNOTATION
================================================================================
*/


(vcfSnpEff, vcfVEP) = vcfGenotypeGVCFs.into(2)


process snpEff {

	publishDir path: "$params.outDir/analysis/snpEff", mode: "copy"
	
	input:
	tuple val(base), file(vcf) from vcfSnpEff
	file(dataDir) from params.snpeffcache
	val(snpeffDB) from params.snpeffversion
	
	output:
	tuple val(base), file("${base}_snpEff.genes.txt"), file("${base}_snpEff.html"), file("${base}_snpEff.csv") into snpeffReport
        tuple val(base), file("${base}_snpEff.ann.vcf") into snpeffVCF
	
	script:
	cache = "-dataDir ${dataDir}"
	"""
	snpEff -Xmx8g \
        ${snpeffDB} \
        -csvStats ${base}_snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${base}_snpEff.ann.vcf
    
    	mv snpEff_summary.html ${base}_snpEff.html
	"""
}


process CompressVCFsnpEff {

    	publishDir path: "$params.outDir/analysis/snpEff", mode: "copy"

    	input:
        tuple val(base), file(vcf) from snpeffVCF

    	output:
        tuple val(base), file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFsnpEffOut

    	script:
    	"""
    	bgzip < ${vcf} > ${vcf}.gz
    	tabix ${vcf}.gz
    	"""
}

process VEP {

    	publishDir path: "$params.outDir/analysis/VEP", mode: "copy"

    	input:
        tuple val(base), file(vcf), file(idx) from vcfVEP
        file(dataDir) from params.vepcache
        val(vepversion) from params.vepversion

    	output:
        tuple val(base), file("${base}_VEP.ann.vcf") into vepVCF
        file("${base}_VEP.summary.html") into vepReport


    	script:
    	"""
    	vep \
    	-i ${vcf} \
    	-o ${base}_VEP.ann.vcf \
    	--assembly GRCh37 \
    	--species homo_sapiens \
    	--cache \
    	--cache_version ${vepversion} \
    	--dir_cache ${dataDir} \
    	--everything \
    	--filter_common \
    	--fork 8 \
    	--format vcf \
    	--per_gene \
    	--stats_file ${base}_VEP.summary.html \
    	--total_length \
    	--vcf
    	"""
}


process CompressVCFvep {

    	publishDir path: "$params.outDir/analysis/VEP", mode: "copy"

    	input:
        tuple val(base), file(vcf) from vepVCF

    	output:
        tuple val(base), file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFOutVEP

    	script:
    	"""
    	bgzip < ${vcf} > ${vcf}.gz
    	tabix ${vcf}.gz
    	"""
}
