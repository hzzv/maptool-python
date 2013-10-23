# -*- coding: utf-8 -*-
# Copyright 2013 by Petra Kubincova

import seq_structures as seq
import sys
import bgzf_tool as bgzftool


# Returns int represented by size/8 bytes on input
def bytes_to_number(stream, offset, size):
    number = 0
    size /= 8
    data = bytearray(size)
    stream.readinto(data)
    for i in range(size):
        number <<= 8
        number += data[i]
    return (number, offset+size)

# Reads header file from stream 'f'
# Returns map from genome name to genome id, list of maps from chromosome
# name to chromosome id, list of indices of tuples with information about
# reference sequences (position on chromosome, bases count, pointer to .bgzf
# file and strand)
def read(f):
    genome_map = {}
    chr_maps = []
    index = []
    i = 0
    # Read number of genomes
    genome_count, i = bytes_to_number(f, i, seq.INF_NUM_SIZE)
    # Read into 'genome_map': name and id
    for j in range(genome_count):
        name_len, i = bytes_to_number(f, i, seq.NAME_SIZE_SIZE)
        name = f.read(name_len)
        i += name_len
        genome_id, i = bytes_to_number(f, i , seq.INF_ID_SIZE)
        genome_map[name] = genome_id

    # Read into 'chr_maps': for each genome and chromosome read
    # chromosome name and chromosome id
    for j in range(len(genome_map)):
        chr_count, i = bytes_to_number(f, i, seq.CHR_ID_SIZE)
        chr_maps.append({})
        for k in range(chr_count):
            name_len, i = bytes_to_number(f, i, seq.NAME_SIZE_SIZE)
            name = f.read(name_len)
            i += name_len
            chr_id, i = bytes_to_number(f, i, seq.CHR_ID_SIZE)
            chr_len, i = bytes_to_number(f, i, seq.LEN_SIZE)
            chr_maps[-1][name] = (chr_id, chr_len)

    # Read into 'index': for each chromosome in reference and each
    # reference block read position on chromosome, count of bases
    # and pointer to bgzf file
    for j in range(len(chr_maps[0])):
        chr_id, i = bytes_to_number(f, i, seq.CHR_ID_SIZE)
        ref_count, i = bytes_to_number(f, i, seq.REF_COUNT_SIZE)
        index.append([])
        for k in range(ref_count):
            strand_number, i = bytes_to_number(f, i, seq.STRAND_SIZE)
            plus_strand = True
            if not strand_number:
                plus_strand = False
            chr_pos, i = bytes_to_number(f, i, seq.POS_SIZE)
            bases_count, i = bytes_to_number(f, i, seq.POS_SIZE)
            pointer, i = bytes_to_number(f, i, seq.FILE_OFFSET_SIZE)
            index[j].append((chr_pos, bases_count, pointer, plus_strand))

    f.close()
    return (genome_map, chr_maps, index)
