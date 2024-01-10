#!/usr/bin/env nextflow

// input 
params.reads = "$projectDir/inputData/*_{1,2}.fq"
params.transcripts = "$projectDir/inputData/transcriptome.fa"
params.outdir  = "$projectDir/results"

// location of programs
params.kallisto="$projectDir/kallisto.yml"
params.fastqc="$projectDir/fastqc.yml"
params.multiqc="$projectDir/multiqc.yml"


// index creation
process INDEX {
    conda "${params.kallisto}"
    
    input:
    path transcriptome

    output:
    path "${transcriptome.baseName}"
   
    script:
    """
    kallisto index -i kallisto_index $transcriptome
    """
}

// Read quantification
process QUANTIFICATION {
    conda "${params.kallisto}"
    
    input:
    path kallisto_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    kallisto quant --threads $task.cpus -i $kallisto_index -o $sample_id ${reads[0]} ${reads[1]}
    """
}

// fastqc report
process FASTQC {
    conda "${params.fastqc}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

// mutiqc report - combine fastqc reports
process MULTIQC {
    conda "${params.multiqc}"
    publishDir "${params.outdir}/multiqc", mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}


// creates channels and input channel to processes
workflow{
    Channel.fromFilePairs(params.reads, checkIfExists: true).set{ readpairs_ch }
    Channel.fromPath(params.transcripts, checkIfExists: true).set{ transcripts_ch }
    index_ch=INDEX(transcripts_ch)
    quant_ch=QUANTIFICATION(index_ch,readpairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}
