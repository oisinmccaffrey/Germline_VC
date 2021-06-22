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
================================================================================
                           Pipeline Parameters
================================================================================
*/

/*
 Reference Genome Files
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
 Exome Intervals, bed files
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
 Annotation tools cache, database versions
*/

params.vep_cache = "/data/VEP/GRCh37"
params.vep_version = "99"
params.snpeff_cache = "/data/snpEff"
params.snpeff_db = "GRCh37.75"


/*
 VEP Plugin files (CADD, LoFTool, ExAC)
*/

params.cadd_wg_snvs = Channel.fromPath("/data/VEP/GRCh37/Plugin_files/whole_genome_SNVs.tsv.gz").getVal()
params.cadd_wg_snvs_tbi = Channel.fromPath("/data/VEP/GRCh37/Plugin_files/whole_genome_SNVs.tsv.gz.tbi").getVal()
params.cadd_indels = Channel.fromPath("/data/VEP/GRCh37/Plugin_files/InDels.tsv.gz").getVal()
params.cadd_indels_tbi = Channel.fromPath("/data/VEP/GRCh37/Plugin_files/InDels.tsv.gz.tbi").getVal()
params.lof = Channel.fromPath("/data/VEP/VEP_plugins/LoFtool_scores.txt").getVal()
params.exac = Channel.fromPath("/data/VEP/VEP_plugins/ExAC.r0.3.1.sites.vep.vcf.gz").getVal()
params.exac_tbi = Channel.fromPath("/data/VEP/VEP_plugins/ExAC.r0.3.1.sites.vep.vcf.gz.tbi").getVal()


/*
 FASTQ reads
*/

params.reads = "/data/omccaffrey/SRR/SRR099957_2.fastq.gz"
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
	tuple val(base), file("${base}.bam") into bamMappedBamQC

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
	file("${base}.bam.metrics") into duplicates_marked_report

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


duplicates_marked_report = duplicates_marked_report.dump(tag:'MarkDuplicates')


process BQSR{

	publishDir path: "$params.outDir/analysis/bqsr", mode: "copy"

	input:
	tuple val(base), file(bam), file(bai) from bam_duplicates_marked
	tuple file(fasta), file(fai), file(dict), file(intlist) from Channel.value([params.fasta, params.fai, params.dict, params.intlist])
	tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
	tuple file(mills), file(millstbi) from Channel.value([params.mills, params.millstbi])

	output:
	tuple val(base), file("${base}.recal.bam"), file("${base}.recal.bam.bai") into BQSR_bams
	tuple val(base), file("${base}.recal.bam") into bam_recalibrated_qc
	file("${base}.recal.stats.out") into samtoolsStatsReport
	file("${base}.recal.table") into baseRecalibratorReport

	script:
	"""
	gatk --java-options -Xmx8g \
	BaseRecalibrator \
	-I $bam \
	-O ${base}.recal.table \
	-L $intlist \
	--tmp-dir . \
	-R $fasta \
	--known-sites $dbsnp \
	--known-sites $mills

	gatk --java-options -Xmx8g \
	ApplyBQSR \
	-I $bam \
	-O ${base}.recal.bam \
	-L $intlist \
	-R $fasta \
	--bqsr-recal-file ${base}.recal.table

	samtools index ${base}.recal.bam ${base}.recal.bam.bai
	samtools stats ${base}.recal.bam > ${base}.recal.stats.out
	"""
}

samtoolsStatsReport = samtoolsStatsReport.dump(tag:'SAMToolsStats')


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
        -D $dbsnp \
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
        -D $dbsnp \
        -V ${gvcf} \
        -O ${base}.vcf
	"""
}


process Split_SNPs_Indels{

	publishDir path: "$params.outDir/analysis/splits", mode: "copy"

	input:
	tuple val(base), file(vcf) from vcfGenotypeGVCFs
	tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

	output:
	tuple val(base), file('*.snps.vcf.gz') into snps_vcf
	tuple val(base), file('*.indels.vcf.gz') into indels_vcf

	script:
	"""
	gatk SelectVariants \
	-R $fasta \
    	-V $vcf \
	-O ${base}.snps.vcf.gz \
    	-select-type SNP

	gatk SelectVariants \
	-R $fasta \
    	-V $vcf \
    	-O ${base}.indels.vcf.gz \
    	-select-type INDEL
	"""
}


process Filter_SNPs{

	publishDir path: "$params.outDir/analysis/splits", mode: "copy"

	input:
	tuple val(base), file(vcf) from snps_vcf
	tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

	output:
	tuple val(base), file("${base}_filtsnps.vcf") into snps_filtered

	script:
	"""
	gatk --java-options -Xmx8g \
	IndexFeatureFile \
        -I ${vcf}

	gatk VariantFiltration \
	-R $fasta \
	-V $vcf \
	-O ${base}_filtsnps.vcf \
	--filter-expression "QD < 2.0" \
	--filter-name "filterQD_lt2.0" \
	--filter-expression "MQ < 25.0" \
	--filter-name "filterMQ_lt25.0" \
	--filter-expression "SOR > 3.0" \
	--filter-name "filterSOR_gt3.0" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "filterMQRankSum_lt-12.5" \
	--filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "filterReadPosRankSum_lt-8.0"
	"""
}


process Filter_Indels{

	publishDir path: "$params.outDir/analysis/splits", mode: "copy"

	input:
	tuple val(base), file(vcf) from indels_vcf
	tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

	output:
	tuple val(base), file("${base}_filtindels.vcf") into indels_filtered

	script:
	"""
	gatk --java-options -Xmx8g \
	IndexFeatureFile \
        -I ${vcf}

	gatk VariantFiltration \
	-R $fasta \
	-V $vcf \
	-O ${base}_filtindels.vcf \
	--filter-expression "QD < 2.0" \
	--filter-name "filterQD" \
	--filter-expression "SOR > 10.0" \
	--filter-name "filterSOR_gt10.0" \
	--filter-expression "ReadPosRankSum < -20.0" \
	--filter-name "filterReadPosRankSum"
	"""
}



process Merge_VCFs {

	publishDir path: "$params.outDir/analysis/splits", mode: "copy"

	input:
	tuple val(base), file(snps) from snps_filtered
	tuple val(base), file(indels) from indels_filtered

	output:
	tuple val(base), file("${base}.vcf.gz") into filtered_vcf

	script:
	"""
	gatk MergeVcfs \
        -I= $snps \
        -I= $indels \
        -O= ${base}.vcf.gz
	"""
}



(vcfsnpEff, bcfstats, vcfstats) = filtered_vcf.into(3)


/*
================================================================================
                                 ANNOTATION
================================================================================
*/


process snpEff{

    publishDir path: "$params.outDir/analysis/snpEff", mode: "copy"

    input:
    tuple val(base), file(vcf) from vcfsnpEff
    val(cache) from params.snpeff_cache
	  val(database) from params.snpeff_db

    output:
    tuple val(base), file("${base}.snpeff.vcf") into snpeff_out

    script:
    cache = "-dataDir ${cache}"
    """
    snpEff -Xmx8G \\
        ${database} \\
        ${cache} \\
        -nostats \\
        -noLog \\
        -lof \\
        -canon \\
        -ud 0 \\
        $vcf > ${base}.snpeff.vcf
    """

}


process EXaC {

    publishDir path: "$params.outDir/analysis/EXaC", mode:'copy'

    input:
    tuple val(base), file(vcf) from snpeff_out

    output:
    tuple val(base), file("${base}.snpeff.exac.vcf") into exac_out

    script:
    """
    java -Xmx4g -jar ${projectDir}/bin/CmdLineAnnotator-1.21.1.jar \\
    -a exac \\
    -s /data/VEP/GRCh37/Plugin_files/ExAC.r0.3.1.sites.vep.vcf.gz \\
    -i $vcf \\
    -o ${base}.snpeff.exac.vcf
    """

}


process CADD {

    publishDir path: "$params.outDir/analysis/CADD", mode:'copy'

    input:
    tuple val(base), file(vcf) from exac_out

    output:
    tuple val(base), file("${base}.snpeff.exac.cadd.vcf") into cadd_out

    script:
    """
    java -Xmx4g -jar ${projectDir}/bin/CmdLineAnnotator-1.21.1.jar \\
    -a cadd \\
    -s /data/VEP/GRCh37/Plugin_files/whole_genome_SNVs.tsv.gz \\
    -i $vcf \\
    -o ${base}.snpeff.exac.cadd.vcf
    """
}

process GAVIN_toCADD {

    publishDir path: "$params.outDir/analysis/upload_to_CADD", mode:'copy'

    input:
    tuple val(base), file(vcf) from cadd_out

    output:
    tuple val(base), file("${base}.toCadd.tsv") into tocadd_out

    script:
    """
    java -Xmx8G -jar ${projectDir}/bin/GAVIN-Plus-1.3beta.jar \\
    -i $vcf \\
    -o ${base}.gavin_firstpass.vcf \\
    -m CREATEFILEFORCADD \\
    -a ${base}.toCadd.tsv \\
    -c ${projectDir}/assets/clinvar.patho.fix.11oct2016.vcf.gz \\
    -d ${projectDir}/assets/CGD_11oct2016.txt.gz \\
    -f ${projectDir}/assets/FDR_allGenes_r1.0.tsv \\
    -g ${projectDir}/assets/GAVIN_calibrations_r0.3.tsv
    """

}
