#!/usr/bin/env nextflow

/*
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quality control
 - 2:   BBDuk for adapter trimming
 - 3: BWA alignment against host reference genome
 - 3.1: Post-alignment processing and format conversion
----------------------------------------------------------------------------------------
*/

def helpMessage() {
	log.info"""
	==============================================
	Heinz Lab ChIP-seq Pipeline
	==============================================
	Usage:

	nextflow run chip-seq-pipeline/main.nf --reads '*.fastq.gz' --genome hg38

	Mandatory arguments:
	 --reads                       Path to input data (must be surrounded with quotes).
	 --genome                      Name of iGenomes reference.

	Options:
	 --singleEnd                   Specifies that the input is single end reads.
	 --allow_multi_align           Secondary alignments and unmapped reads are also reported in addition to primary alignments.
	 --extendReadsLen [int]        Length to extend reads. (Default: 100)

	Trimming options:
	 --length [int]                Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18

	Other options:
	 --outdir                      The output directory where the results will be saved. Default: results
	 --skip_qc                     Skip all QC steps aside from MultiQC.
	 --skip_fastqc                 Skip FastQC.
	 --list_genomes                Similar to --help. Used to list all available genome keywords.

	References:
	 --saveReference               Save the generated reference files the the Results directory.
	 --saveAlignedIntermediates    Save the intermediate BAM files from the Alignment step prior to sorting - not done by default.
	 --fasta                       Path to fasta reference.
	""".stripIndent()
}

def genome_listMessage() {
	log.info"""
	==============================================
	              Available Genomes
	==============================================

	Homo sapiens:
	 hg38                          UCSC
	 GRCh37                        Ensembl

	Mus musculus:
	 mm10                          UCSC
	 GRCm38                        Ensembl

	Rattus norvegicus:
	 rn6                           UCSC
	 Rnor_6.0                      Ensembl
	""".stripIndent()
}

// Show help message
params.help = false
if (params.help){
	helpMessage()
	exit 0
}

// Show list of genomes
params.list_genomes = false
if (params.list_genomes){
	genome_listMessage()
	exit 0
}

// Reference path configuration
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.adapters = params.genome ? params.genomes[ params.genome ].adapters ?: false: false
params.chrom_sizes = params.genome ? params.genomes[ params.genome ].chrom_sizes ?: false : false

multiqc_config = file(params.multiqc_config)
chrom_sizes = file(params.chrom_sizes)

// Validate inputs
if( params.adapters ){
	adapters = file(params.adapters)
	if ( !adapters.exists() ) exit 1, "adapters file not found: ${params.adapters}"
}
if( params.gtf ){
	gtf = file(params.gtf)
	if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}
if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
} else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
} else {
    exit 1, "No reference genome specified!"
}

custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Create a channel for input reads
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_bbduk }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_bbduk }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_bbduk }
}

// Header log info
log.info """=======================================================
ChIP-seq pipeline - svenner lab
======================================================="""
def summary = [:]
summary['Run Name']            = custom_runName ?: workflow.runName
summary['Reads']               = params.reads
summary['Genome']              = params.genome
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Trim Min Length']     = params.length
if(params.bwa_index)           summary['BWA Index'] = params.bwa_index
else if(params.fasta)          summary['Fasta Ref'] = params.fasta
if(params.gtf)                 summary['GTF Annotation'] = params.gtf
summary['Adapters']            = params.adapters
summary['Chrom Sizes']         = params.chrom_sizes
summary['Multiple alignments'] = params.allow_multi_align
summary['Extend Reads']        = "$params.extendReadsLen bp"
summary['Output dir']          = params.outdir
summary['Working dir']         = workflow.workDir
summary['Container']           = params.container
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Script dir']          = workflow.projectDir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================="

/*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index && fasta){
    process makeBWAindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BWAIndex" into bwa_index

        script:
        """
        bwa index -a bwtsw $fasta
        mkdir BWAIndex && mv ${fasta}* BWAIndex
        """
    }
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - BBDuk - Not optimized for pairEnded
 */
process bbduk {
    tag "$name"
    publishDir "${params.outdir}/bbduk", mode: 'copy'

    input:
    set val(name), file(reads) from raw_reads_bbduk
    file adapters from adapters

    output:
    file '*.gz' into trimmed_reads

    script:
    tg_length = "--length ${params.length}"
    c_r1 = params.clip_R1 > 0 ? "--clip_R1 ${params.clip_R1}" : ''
    tpc_r1 = params.three_prime_clip_R1 > 0 ? "--three_prime_clip_R1 ${params.three_prime_clip_R1}" : ''
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    if (params.singleEnd) {
        """
        bbduk.sh in=$reads out=${prefix}_trimmed.fastq.gz ref=$adapters -ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20
        """
    } else {
        """
        bbduk.sh in=$reads out=${prefix}_trimmed.fastq.gz ref=$adapters -ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 $reads[0] $reads[1]
        """
    }
}

/*
 * STEP 3.1 - align with bwa
 */
process bwa {
    tag "$prefix"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bwa" : params.outdir }, mode: 'copy',
               saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file reads from trimmed_reads
    file index from bwa_index.first()

    output:
    file '*.bam' into bwa_bam

    script:
    prefix = reads[0].toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "| samtools view -b -q 1 -F 4 -F 256"
    """
    bwa mem -t 8 -M ${index}/genome.fa $reads | samtools view -bT $index - $filtering > ${prefix}.bam
    """
}

/*
 * STEP 3.2 - post-alignment processing
 */

process samtools {
    tag "${bam.baseName}"
    publishDir path: "${params.outdir}/bwa", mode: 'copy',
               saveAs: { filename ->
                   if (filename.indexOf(".stats.txt") > 0) "stats/$filename"
                   else params.saveAlignedIntermediates ? filename : null
               }

    input:
    file bam from bwa_bam

    output:
    file '*.sorted.bam' into bam_picard, bam_for_mapped
    file '*.sorted.bam.bai' into bwa_bai, bai_for_mapped
    file '*.sorted.bed' into bed_total
    file '*.stats.txt' into samtools_stats

    script:
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    bedtools bamtobed -i ${bam.baseName}.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${bam.baseName}.sorted.bed
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}.stats.txt
    """
}


/*
 * STEP 3.3 - Statistics about mapped and unmapped reads against ref genome
 */

process bwa_mapped {
    tag "${input_files[0].baseName}"
    publishDir "${params.outdir}/bwa/mapped", mode: 'copy'

    input:
    file input_files from bam_for_mapped.collect()
    file bai from bai_for_mapped.collect()

    output:
    file 'mapped_refgenome.txt' into bwa_mapped

    script:
    """
    for i in $input_files
    do
      samtools idxstats \${i} | awk -v filename="\${i}" '{mapped+=\$3; unmapped+=\$4} END {print filename,"\t",mapped,"\t",unmapped}'
    done > mapped_refgenome.txt
    """
}

/*
 * STEP 4 Picard
 */

process picard {
    tag "$prefix"
    publishDir "${params.outdir}/picard", mode: 'copy'

    input:
    file bam from bam_picard

    output:
    file '*.dedup.sorted.bam' into bam_dedup_ssp, bam_dedup_ngsplot, bam_dedup_deepTools, bam_dedup_macs, bam_dedup_saturation
    file '*.dedup.sorted.bam.bai' into bai_dedup_deepTools, bai_dedup_ngsplot, bai_dedup_macs, bai_dedup_ssp
    file '*.dedup.sorted.bed' into bed_dedup
    file '*.picardDupMetrics.txt' into picard_reports

    script:
    prefix = bam[0].toString() - ~/(\.sorted)?(\.bam)?$/
    if( task.memory == null ){
        log.warn "[Picard MarkDuplicates] Available memory not known - defaulting to 6GB ($prefix)"
        avail_mem = 6000
    } else {
        avail_mem = task.memory.toMega()
        if( avail_mem <= 0){
            avail_mem = 6000
            log.warn "[Picard MarkDuplicates] Available memory 0 - defaulting to 6GB ($prefix)"
        } else if( avail_mem < 250){
            avail_mem = 250
            log.warn "[Picard MarkDuplicates] Available memory under 250MB - defaulting to 250MB ($prefix)"
        }
    }
    """
    picard MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${prefix}.dedup.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=${prefix}.picardDupMetrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        PROGRAM_RECORD_ID='null'
    samtools sort ${prefix}.dedup.bam -o ${prefix}.dedup.sorted.bam
    samtools index ${prefix}.dedup.sorted.bam
    bedtools bamtobed -i ${prefix}.dedup.sorted.bam | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.dedup.sorted.bed
    """
}

/*
 * STEP 5 Read_count_statistics
 */

process countstat {
    tag "${input[0].baseName}"
    publishDir "${params.outdir}/countstat", mode: 'copy'

    input:
    file input from bed_total.mix(bed_dedup).toSortedList()

    output:
    file 'read_count_statistics.txt' into countstat_results

    script:
    """
    countstat.pl $input
    """
}

/*
 * STEP 6 ssp quality metrics
 */

process ssp {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/ssp", mode: "copy"

    when:
    !params.skip_qc

    input:
    file bam from bam_dedup_ssp
    file bai from bai_dedup_ssp
    file chrom_sizes

    output:
    file 'sspout/*.{txt,pdf}' into ssp_results

    script:
    prefix = bam[0].toString() - ~/(\.sorted)?(\.bam)?$/
    if (!params.singleEnd) {
        """
        ssp -i $bam -o ${prefix} --gt $chrom_sizes -p 2 --pair
        """
    } else {
        """
        ssp -i $bam -o ${prefix} --gt $chrom_sizes -p 2
        """
    }
}

/*
 * STEP 7 deeptools quality metrics
 */

process deepTools {
    tag "${bam[0].baseName}"
    publishDir "${params.outdir}/deepTools", mode: 'copy'

    input:
    file bam from bam_dedup_deepTools.toSortedList()
    file bai from bai_dedup_deepTools.toSortedList()

    output:
    file '*.{txt,pdf,png,npz,bw}' into deepTools_results
    file '*.txt' into deepTools_multiqc

    script:
    if (!params.singleEnd) {
        """
        bamPEFragmentSize \\
            --binSize 1000 \\
            --bamfiles $bam \\
            --histogram fragment_length_distribution_histogram.png \\
            --plotTitle "Fragment Length Distribution"
        """
    }
    if(bam instanceof Path){
        log.warn("Only 1 BAM file - skipping multiBam deepTool steps")
        """
        plotFingerprint \\
            -b $bam \\
            --plotFile ${bam.baseName}_fingerprints.pdf \\
            --outRawCounts ${bam.baseName}_fingerprint.txt \\
            --extendReads ${params.extendReadsLen} \\
            --skipZeros \\
            --ignoreDuplicates \\
            --numberOfSamples 50000 \\
            --binSize 500 \\
            --plotFileFormat pdf \\
            --plotTitle "${bam.baseName} Fingerprints" \\
            --outQualityMetrics ${bam.baseName}_fingerprints_metrics.txt \\
            --smartLabels
        bamCoverage \\
           -b $bam \\
           --extendReads ${params.extendReadsLen} \\
           --normalizeUsing RPKM \\
           -o ${bam}.bw
        """
    } else {
        """
        plotFingerprint \\
            -b $bam \\
            --plotFile fingerprints.pdf \\
            --outRawCounts fingerprint.txt \\
            --extendReads ${params.extendReadsLen} \\
            --skipZeros \\
            --ignoreDuplicates \\
            --numberOfSamples 50000 \\
            --binSize 500 \\
            --plotFileFormat pdf \\
            --plotTitle "Fingerprints" \\
            --outQualityMetrics all_fingerprints_metrics.txt \\
            --smartLabels
        for bamfile in ${bam}
        do
            bamCoverage \\
              -b \$bamfile \\
              --extendReads ${params.extendReadsLen} \\
              --normalizeUsing RPKM \\
              -o \${bamfile}.bw
        done
        multiBamSummary \\
            bins \\
            --binSize 10000 \\
            --bamfiles $bam \\
            -out multiBamSummary.npz \\
            --extendReads ${params.extendReadsLen} \\
            --ignoreDuplicates \\
            --centerReads
        plotCorrelation \\
            -in multiBamSummary.npz \\
            -o scatterplot_PearsonCorr_multiBamSummary.png \\
            --outFileCorMatrix scatterplot_PearsonCorr_multiBamSummary.txt \\
            --corMethod pearson \\
            --skipZeros \\
            --removeOutliers \\
            --plotTitle "Pearson Correlation of Read Counts" \\
            --whatToPlot scatterplot
        plotCorrelation \\
            -in multiBamSummary.npz \\
            -o heatmap_SpearmanCorr_multiBamSummary.png \\
            --outFileCorMatrix heatmap_SpearmanCorr_multiBamSummary.txt \\
            --corMethod spearman \\
            --skipZeros \\
            --plotTitle "Spearman Correlation of Read Counts" \\
            --whatToPlot heatmap \\
            --colorMap RdYlBu \\
            --plotNumbers
        plotPCA \\
            -in multiBamSummary.npz \\
            -o pcaplot_multiBamSummary.png \\
            --plotTitle "Principal Component Analysis Plot" \\
            --outFileNameData pcaplot_multiBamSummary.txt
        """
    }
}

/*
 * STEP 7 MultiQC
 */

process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('samtools/*') from samtools_stats.collect()
    file ('picard/*') from picard_reports.collect()
    file ('deeptools/*') from deepTools_multiqc.collect()

    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config . 2>&1
    """
}
