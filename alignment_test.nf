#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.fastq_folder = "/path/to/fastqs/**.fastq.gz"
params.output_dir = "path/to/output/directory"
params.ref = "path/to/reference"

/*
 * Process to merge multiple fastqs into one single fastq. 
 */
process merge_fastq {
    input:
    path x

    output:
    path merged.fastq.gz

    ```
    copy /b *.fastq.gz $params.output_dir/merged.fastq.gz
    ```
}

/*
 * Process to align reads to reference using bwa mem and sort and index the bam file
 */
process align {

    conda '/path/to/an/existing/env/directory'
    cpus 8
    executor 'sge'
    memory '8 GB'
    queue 'all.q'

    input:
    path merged.fastq.gz

    output:
    path aligned.sorted.bam

    ```
    bwa index -p prefix $params.ref
    bwa mem -M -t threads reference merged.fastq.gz > aligned.sam
    samtools fixmate -O bam aligned.sam aligned.bam
    samtools sort -T temp.prefix -O bam -@ threads -o aligned.bam aligned.sorted.bam
    ```   
}

process mark_duplicates {

    conda '/path/to/an/existing/env/directory'
    cpus 8
    executor 'sge'
    memory '8 GB'
    queue 'all.q'

    input:
    path aligned.sorted.bam

    output:
    path merged.sorted.bam

    ```
    java -Xmx16g -jar picard.jar MarkDuplicates ASSUME_SORTED=true INPUT=path aligned.sorted.bam OUTPUT=merged.sorted.bam METRICS_FILE=metrics.out
    samtools index -b merged.sorted.bam
    ```
}

workflow {
    reads = Channel.fromPath( $params.fastq_folder )
    
}