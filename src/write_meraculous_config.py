#!/usr/bin/env python3

from pathlib import Path


def write_config_file(
        fastq,
        k,
        diplo_mode,
        dmin,
        threads,
        config_string,
        config_file):
    '''
    Accept fastq file, threads config string and output location and write
    config
    '''
    my_fastq = Path(fastq).resolve().as_posix()
    my_conf = config_string.format(my_fastq, k, diplo_mode, dmin, threads)
    with open(config_file, 'wt') as f:
        f.write(my_conf)
    return True


def main():
    with open(meraculous_config_file, 'rt') as f:
        meraculous_config_string = ''.join(f.readlines())

    write_config_file(
        snakemake.input['fq'],
        snakemake.wildcards['k'],
        snakemake.wildcards['diplo'],
        snakemake.params['dmin'],
        snakemake.params['threads'],
        meraculous_config_string,
        snakemake.output['config'])


meraculous_config_file = snakemake.input['config_file']

if __name__ == '__main__':
    main()

