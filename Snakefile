#!/usr/bin/env python3

import pandas


#############
# FUNCTIONS #
#############


def convert_dmin(from_k, to_k, dmin_from_k, rl):
    '''
    Convert dmin_from_k coverage at from_k to coverage at to_k based on a read
    length of rl
    '''
    dmin_to_k = (dmin_from_k * (to_k - rl - 1)) / (from_k - rl - 1)
    return int(round(dmin_to_k))


def find_completed_assemblies():
    assembly_path = ('output/050_meraculous/'
                     '{region}_k{k}_diplo{diplo}/'
                     'meraculous_final_results/final.scaffolds.fa')
    my_wildcards = snakemake.io.glob_wildcards(assembly_path)
    my_wc_dict = my_wildcards._asdict()
    return(my_wc_dict)


def generate_stats_targets(wildcards):
    completed_assembly_dict = find_completed_assemblies()
    return expand('output/060_stats/{region}.{diplo}.{k}/stats.tsv',
                  zip,
                  **completed_assembly_dict)


###########
# GLOBALS #
###########

max_threads = workflow.cores

assembly = 'data/GCA_014466185.1_ASM1446618v1_genomic.fna'
index = 'data/GCA_014466185.1_ASM1446618v1_genomic.fna.fai'

reads = 'data/vvul.fq.gz'

meraculous_config_file = 'src/meraculous_config.txt'
dmin_file = 'data/dmin.csv'

bbmap = 'https://github.com/deardenlab/container-bbmap/releases/download/0.0.1/container-bbmap.bbmap_38.90.sif'
biopython = 'shub://TomHarrop/py-containers:biopython_1.78'
bwa = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
meraculous = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
pigz = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
r = 'shub://TomHarrop/r-containers:r_4.0.0'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'

########
# MAIN #
########

# get a list of regions
all_regions = sorted(set(pandas.read_csv(index, sep='\t', header=None)[0]))

# csv of individual and dmin at k==31
indiv_dmin = pandas.read_csv(dmin_file,
                             index_col='chr')


rule target:
    input:
        # first set
        # expand(('output/050_meraculous/'
        #         '{region}_k{k}_diplo{diplo}/'
        #         'meraculous_final_results/final.scaffolds.fa'),
        #        region=all_regions,
        #        k=31,
        #        diplo=[0, 1])
        # second set
        expand(('output/050_meraculous/'
                '{region}_k{k}_diplo{diplo}/'
                'meraculous_final_results/final.scaffolds.fa'),
               region=['CM025407.1'],
               k=[31, 51, 71],
               diplo=[0, 1])


#####################
# ASSEMBLY QC RULES #
#####################

# if you ask for assembly_qc, generate_qc_targets() will call
# find_completed_assemblies() and generate assembly_stats and busco targets for
# each completed assembly it finds
rule assembly_qc:
    input:
        # 'output/060_stat-plots/busco_results.pdf',
        'output/070_stat-plots/assembly_stats.pdf'

rule plot_assembly_stats:
    input:
        stats_files = generate_stats_targets
    output:
        plot = 'output/070_stat-plots/assembly_stats.pdf'
    log:
        'output/logs/plot_assembly_stats.log'
    singularity:
        r
    script:
        'src/plot_assembly_stats.R'

rule assembly_stats:
    input:
        contigs = ('output/050_meraculous/'
                   '{region}_k{k}_diplo{diplo}/'
                   'meraculous_final_results/final.scaffolds.fa')
    output:
        stats = 'output/060_stats/{region}.{diplo}.{k}/stats.tsv'
    log:
        'output/logs/assembly_stats.{region}.{diplo}.{k}.log'
    threads:
        1
    container:
        bbmap
    shell:
        'stats.sh '
        'in={input} '
        'format=3 '
        'threads={threads} '
        '> {output} '
        '2> {log}'

# try per-region assembly
rule meraculous:
    input:
        fq = 'output/040_read-chunks/chunk_{region}.fq',
        config = ('output/050_meraculous/'
                  '{region}_k{k}_diplo{diplo}/config.txt'),
    output:
        ('output/050_meraculous/'
         '{region}_k{k}_diplo{diplo}/'
         'meraculous_final_results/final.scaffolds.fa')
    params:
        outdir = 'output/050_meraculous/{region}_k{k}_diplo{diplo}/',
        dmin = '0'
    threads:
        workflow.cores
    log:
        'output/logs/meraculous.{region}.{k}.{diplo}.log'
    container:
        meraculous
    shell:
        'run_meraculous.sh '
        '-dir {params.outdir} '
        '-config {input.config} '
        '-cleanup_level 2 '
        '&> {log}'


rule write_meraculous_config:
    input:
        fq = 'output/040_read-chunks/chunk_{region}.fq',
        config_file = meraculous_config_file
    output:
        config = ('output/050_meraculous/'
                  '{region}_k{k}_diplo{diplo}/config.txt'),
    threads:
        1
    params:
        dmin = lambda wildcards:
            convert_dmin(
                from_k=31,
                to_k=int(wildcards.k),
                dmin_from_k=int(indiv_dmin.loc[wildcards.region, 'dmin_31']),
                rl=int(indiv_dmin.loc[wildcards.region, 'rl'])),
        threads = max_threads - 4
    container:
        biopython
    script:
        'src/write_meraculous_config.py'

rule retrieve_reads:
    input:
        read_ids = expand('output/030_bam-chunks/chunk_{region}.txt',
                          region=all_regions),
        fastq = 'output/000_tmp/pe_reads.fq'
    output:
        expand('output/040_read-chunks/chunk_{region}.fq',
               region=all_regions)
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


rule fa_index:
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

