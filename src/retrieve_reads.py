#!/usr/bin/env python3

import logging
import os
from Bio import SeqIO

# set up log
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename=snakemake.log[0],
    level=logging.INFO)

# catch files from snakemake
read_id_list = snakemake.input['read_ids']
read_file = snakemake.input['fastq']
outdir = snakemake.params['outdir']
# read_no = snakemake.params['read_no']

# dev
# read_id_list = ['output/040_read-chunks/chunk_87.txt',
#                 'output/040_read-chunks/chunk_999.txt']
# read_file = 'test/r2.fq'
# outdir = 'test'
# read_no = '2'

# dict of chunk to outfile
chunk_to_outfile = {os.path.basename(x).rstrip('.txt'):
                    os.path.join(outdir,
                                 os.path.basename(x).rstrip('.txt') + '.fq')
                    for x in read_id_list}

# initialise a dict of read_id to chunk
read_to_chunk = dict()

# loop over read_id_list 
for read_id_file in read_id_list:
    my_chunk_id = os.path.basename(read_id_file).rstrip('.txt')
    logging.info(f'Processing {my_chunk_id}')
    # get the list of read ids
    with open(read_id_file, 'rt') as f:
        ids = [x.rstrip('\n') for x in f.readlines()]
    # check if id is already indexed
    for id in ids:
        if id in read_to_chunk:
            read_to_chunk[id].append(my_chunk_id)
        else:
            read_to_chunk[id] = [my_chunk_id]


# open a handle for each chunk
logging.info(f'Opening file handles')
chunk_to_handle = {x: open(chunk_to_outfile[x], 'wt')
                   for x in chunk_to_outfile}

# read through the fastq and write each file to the handles
logging.info(f'Started reading {read_file}')
i = 1
for seq_rec in SeqIO.parse(read_file, 'fastq'):
    if i % 1e6 == 0:
        logging.info(f'Processed {int(i / 1e6)} million reads')
    i += 1
    try:
        write_chunks = read_to_chunk[seq_rec.id]
        for chunk in write_chunks:
            SeqIO.write(seq_rec, chunk_to_handle[chunk], 'fastq')
    except KeyError as e:
        logging.debug(f'Read {seq_rec.id} not in dict(read_to_chunk)')

# close all the handles
logging.info(f'Finished {read_file}')
logging.info(f'Closing file handles')
for chunk in chunk_to_handle:
    chunk_to_handle[chunk].close()
