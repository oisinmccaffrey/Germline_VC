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

params.reads = "/data/bdigby/WES/reads/METIN_YUSA_EZELSOY_R{1,2}_001.fastq.gz"
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


/*
  !!! omitting bed intervals until supplied with Illumina capture kit 
  !!! -L option from analysis until discussed with Pilib 
  !!! (BQSR, ApplyBQSR, haplotypecaller, genotypecaller had -L flags)
*/


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
	file(fasta) from params.fasta
	
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
	file(fasta) from params.fasta
	
	output:
	tuple val(base), file("${base}_filtsnps.vcf") into snps_filtered
	
	script:
	"""
	gatk VariantFiltration \
	-R $fasta \
	-V $vcf \
	-O ${base}_filtsnps.vcf \
	--filterExpression "QD < 2.0" \
	--filterName "filterQD_lt2.0" \
	--filterExpression "MQ < 25.0" \
	--filterName "filterMQ_lt25.0" \
	--filterExpression "SOR > 3.0" \
	--filterName "filterSOR_gt3.0" \
	--filterExpression "MQRankSum < -12.5" \
	--filterName "filterMQRankSum_lt-12.5" \
	--filterExpression "ReadPosRankSum < -8.0" \
	--filterName "filterReadPosRankSum_lt-8.0"
	"""
}


process Filter_Indels{

	publishDir path: "$params.outDir/analysis/splits", mode: "copy"
	
	input:
	tuple val(base), file(vcf) from indels_vcf
	file(fasta) from params.fasta
	
	output:
	tuple val(base), file("${base}_filtindels.vcf") into indels_filtered
	
	script:
	"""
	gatk VariantFiltration \
	-R $fasta \
	-V $vcf \
	-O ${base}_filtindels.vcf \
	--filterExpression "QD < 2.0" \
	--filterName "filterQD" \
	--filterExpression "SOR > 10.0" \
	--filterName "filterSOR_gt10.0" \
	--filterExpression "ReadPosRankSum < -20.0" \
	--filterName "filterReadPosRankSum"
	"""
}



process Merge_VCFs {

	publishDir path: "$params.outDir/analysis/splits", mode: "copy"
	
	input:
	tuple val(base), file(snps) from snps_filtered
	tuple val(base), file(indels) from indels_filtered
	
	output:
	tuple val(base), file(vcf) into filtered_vcf
	
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

/*
 Annotation Strategy: 
 	Run VCF through snpEff, then run through VEP with Plugins for rich 
  	annotation and GAVIN compatability (requires snpEff, CADD, ExAC ann)
*/


process snpEff {

	publishDir path: "$params.outDir/analysis/snpEff", mode: "copy"
	
	input:
	tuple val(base), file(vcf) from vcfsnpEff
	val(cache) from params.snpeff_cache
	val(database) from params.snpeff_db
	
	output:
	set file("${base}_snpEff.genes.txt"), file("${base}_snpEff.html"), file("${base}_snpEff.csv") into snpeffReport
        tuple val(base), file("${base}_snpEff.ann.vcf") into snpeffVCF

	script:
	cache = "-dataDir ${cache}"
	"""
	snpEff -Xmx8g \
        ${database} \
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


snpeffReport = snpeffReport.dump(tag:'snpEff')


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


process VEPsnpEff {

    	publishDir path: "$params.outDir/analysis/snpEff", mode: "copy"

    	input:
        tuple val(base), file(vcf), file(vcf_tbi) from compressVCFsnpEffOut
        val(dataDir) from params.vep_cache
        val(vepversion) from params.vep_version
	file(fasta) from params.fasta
	tuple file(cadd_snv), file(cadd_snv_tbi) from Channel.value([params.cadd_wg_snvs, params.cadd_wg_snvs_tbi])
	tuple file(cadd_indels), file(cadd_indels_tbi) from Channel.value([params.cadd_indels, params.cadd_indels_tbi])
	tuple file(exac), file(exac_tbi) from Channel.value([params.exac, params.exac_tbi])
	file(lof) from params.lof

	
    	output:
        tuple val(base), file("${base}_VEP.ann.vcf") into vepVCF
        file("${base}_VEP.summary.html") into vepReport


    	script:
	ExAC = "--plugin ExAC,ExAC.r0.3.1.sites.vep.vcf.gz"
	CADD = "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz"
	LOF = "--plugin LoFtool,LoFtool_scores.txt"
	genesplicer = "--plugin GeneSplicer,/opt/conda/envs/Germline_VC/bin/genesplicer,/opt/conda/envs/Germline_VC/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${base}"
    	"""
    	vep \
    	-i ${vcf} \
    	-o ${base}_VEP.ann.vcf \
    	--assembly GRCh37 \
    	--species homo_sapiens \
	${ExAC} \
	${CADD} \
	${LOF} \
	${genesplicer} \
	--offline \
    	--cache \
	--fasta $fasta \
    	--cache_version ${vepversion} \
    	--dir_cache ${dataDir} \
    	--everything \
    	--filter_common \
    	--fork 4 \
    	--format vcf \
    	--per_gene \
    	--stats_file ${base}_VEP.summary.html \
    	--total_length \
    	--vcf
	
	rm -rf ${base}
    	"""
}


vepReport = vepReport.dump(tag:'VEP')


process CompressVCFvep {

    	publishDir path: "$params.outDir/analysis/combined_annot", mode: "copy"

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



/*
================================================================================
				Quality Control
================================================================================
*/

process BamQC {

    	publishDir path: "$params.outDir/analysis/bamQC", mode: "copy"

    	input:
    	tuple val(base), file(bam) from bam_recalibrated_qc
    	file(targetBED) from params.bed

    	output:
    	file("${bam.baseName}") into bamQCReport

    	script:
    	use_bed = "-gff ${targetBED}"
    	"""
    	qualimap --java-mem-size=8G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        $use_bed \
        -nt 8 \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${bam.baseName} \
        -outformat HTML
    	"""
}


bamQCReport = bamQCReport.dump(tag:'BamQC')


process BcftoolsStats {

    	publishDir path: "$params.outDir/analysis/quality", mode: "copy"

    	input:
    	tuple val(base), file(vcf) from bcfstats

    	output:
    	file ("*.bcf.tools.stats.out") into bcftoolsReport

    	script:
    	"""
    	bcftools stats ${vcf} > ${base}.bcf.tools.stats.out
    	"""
}


bcftoolsReport = bcftoolsReport.dump(tag:'BCFTools')


process Vcftools {

    	publishDir path: "$params.outDir/analysis/quality", mode: "copy"

    	input:
    	tuple val(base), file(vcf) from vcfstats

    	output:
    	file ("${base}.*") into vcftoolsReport

    	script:
    	"""
    	vcftools \
    	--gzvcf ${vcf} \
    	--TsTv-by-count \
    	--out ${base}.vcf
    
    	vcftools \
    	--gzvcf ${vcf} \
    	--TsTv-by-qual \
    	--out ${base}.vcf
    
    	vcftools \
    	--gzvcf ${vcf} \
    	--FILTER-summary \
    	--out ${base}.vcf
    	"""
}


vcftoolsReport = vcftoolsReport.dump(tag:'VCFTools')



process MultiQC {

    	publishDir path: "$params.outDir/analysis/MultiQC", mode: "copy"

    	input:
        file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
        file ('BCFTools/*') from bcftoolsReport.collect().ifEmpty([])
        file ('MarkDuplicates/*') from duplicates_marked_report.collect().ifEmpty([])
        file ('DuplicatesMarked/*.recal.table') from baseRecalibratorReport.collect().ifEmpty([])
        file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
        file ('snpEff/*') from snpeffReport.collect().ifEmpty([])
        file ('VCFTools/*') from vcftoolsReport.collect().ifEmpty([])

    	output:
    	file ("*multiqc_report.html") into ch_multiqc_report
    	file ("*_data")

    	script:
    	rtitle = "--title Galway_Genomics"
	rfilename = "--filename Galway_Genomics_multiqc_report"
    	"""
    	multiqc -f ${rtitle} ${rfilename} .
    	"""
}
