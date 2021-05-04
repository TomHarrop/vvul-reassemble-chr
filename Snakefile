#!/usr/bin/env python3

max_threads = workflow.cores

assembly = 'data/GCA_014466185.1_ASM1446618v1_genomic.fna'
reads = 'data/vvul.fq.gz'

pigz = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
bwa = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'


rule markdup:
    input:
        'output/020_alignment/aln.sort.bam'
    output:
        'output/020_alignment/aln.bam'
    log:
        'output/logs/markdup.log'
    container:
        samtools
    shell:
        'samtools markdup '
        f'-@ {max_threads // 3} '
        '-s '
        '{input} '
        '{output} '
        '2> {log}'

rule sort:
    input:
        'output/020_alignment/aln.fixmate.bam'
    output:
        pipe('output/020_alignment/aln.sort.bam')
    log:
        'output/logs/sort.log'
    container:
        samtools
    shell:
        'samtools sort '
        f'-@ {max_threads // 3} '
        '{input} '
        '>> {output} '
        '2> {log}'

rule fixmate:
    input:
        'output/020_alignment/aln.sam'
    output:
        pipe('output/020_alignment/aln.fixmate.bam')
    log:
        'output/logs/fixmate.log'
    container:
        samtools
    shell:
        'samtools fixmate '
        '-m '
        f'-@ {max_threads // 3} '
        '{input} '
        '- '
        '>> {output} '
        '2> {log}'

rule map_reads:
    input:
        expand('output/010_index/index.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        reads = 'output/000_tmp/pe_reads.fq'
    output:
        temp('output/020_alignment/aln.sam')
    params:
        prefix = 'output/010_index/index'
    log:
        'output/logs/map_readsmem.log'
    threads:
        workflow.cores
    container:
        bwa
    shell:
        'bwa mem '
        '-t {threads} '
        '-p '
        '{params.prefix} '
        '{input.reads} '
        '> {output} '
        '2> {log}'

rule index_assembly:
    input:
        fasta = assembly
    output:
        expand('output/010_index/index.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/010_index/index'
    log:
        'output/logs/index_assembly.log'
    threads:
        1
    container:
        bwa
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input.fasta} '
        '2> {log} '


rule tmp_unzip_short:
    input:
        reads
    output:
        temp('output/000_tmp/pe_reads.fq')
    threads:
        3
    container:
        pigz
    shell:
        'pigz -dc '
        '{input} '
        '>{output}'