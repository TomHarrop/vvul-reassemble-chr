#!/usr/bin/env python3

import pandas


def rr_input(wildcards):
    fai = checkpoints.fa_index.get(**wildcards).output['fai']
    r = sorted(set(pandas.read_csv(fai, sep='\t', header=None)[0]))
    return expand('output/030_bam-chunks/chunk_{region}.txt',
                  region=r)


def rr_output(wildcards):
    fai = checkpoints.fa_index.get(**wildcards).output['fai']
    r = sorted(set(pandas.read_csv(fai, sep='\t', header=None)[0]))
    return expand('output/040_read-chunks/{region}.fq',
                  region=r)


max_threads = workflow.cores

assembly = 'data/GCA_014466185.1_ASM1446618v1_genomic.fna'
reads = 'data/vvul.fq.gz'

pigz = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
bwa = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
biopython = 'shub://TomHarrop/py-containers:biopython_1.78'


rule target:
    input:
        'output/040_read-chunks/a_file'

rule retrieve_reads:
    input:
        read_ids = rr_input,
        fastq = 'output/000_tmp/pe_reads.fq'
    output:
        'output/040_read-chunks/a_file'
    params:
        outdir = 'output/040_read-chunks',
    log:
        'output/logs/retrieve_reads.log'
    container:
        biopython
    script:
        'src/retrieve_reads.py'

# get the list of read ids for each chunk
rule extract_read_ids:
    input:
        'output/030_bam-chunks/chunk_{region}.sam'
    output:
        'output/030_bam-chunks/chunk_{region}.txt'
    log:
        'output/logs/extract_read_ids.{region}.log'
    container:
        samtools
    shell:
        'samtools view  {input} '
        '| cut -f1 '
        '| sort '
        '| uniq '
        '> {output} '
        '2> {log}'


# subset the BAM by the chunk list
rule chunk_bam:
    input:
        bam = 'output/020_alignment/aln.bam',
        bai = 'output/020_alignment/aln.bam.bai'
    output:
        temp('output/030_bam-chunks/chunk_{region}.sam')
    log:
        'output/logs/chunk_bam.{region}.log'
    container:
        samtools
    shell:
        'samtools view '
        '-h '
        '-F 256 '       # exclude secondary alignments
        '-O SAM '
        '{input.bam} '
        '{wildcards.region} '
        '> {output} '
        '2> {log}'


rule index_bamfile:
    input:
        'output/020_alignment/aln.bam'
    output:
        'output/020_alignment/aln.bam.bai'
    log:
        'output/logs/index_bamfile.log'
    threads:
        2
    container:
        samtools
    shell:
        'samtools index -@ {threads} {input} 2> {log}'

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
        fasta = 'output/010_index/index.fasta'
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


checkpoint fa_index:
    input:
        assembly
    output:
        fa = 'output/010_index/index.fasta',
        fai = 'output/010_index/index.fasta.fai'
    container:
        samtools
    shell:
        'cp {input} {output.fa} '
        '; '
        'samtools faidx {output.fa}'

rule tmp_unzip_short:
    input:
        reads
    output:
        'output/000_tmp/pe_reads.fq'
    threads:
        3
    container:
        pigz
    shell:
        'pigz -dc '
        '{input} '
        '>{output}'

