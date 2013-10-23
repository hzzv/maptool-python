# -*- coding: utf-8 -*-
# Copyright 2013 by Petra Kubincova

import array
import sys
import seq_structures as seq
import bgzf_tool as bgzftool


# Correctly adds genome name 'name' to 'genome_map' if it wasn't already added
def add_genome(genome_map, chr_maps, name):
    if not name[0] in genome_map:
        genome_map[name[0]] = len(genome_map)
        while genome_map[name[0]] >= len(chr_maps):
            chr_maps.append({})

# Parses line 'line' from file of .maf format
# Returns Reference structure or id of informant and Informant structure
def parse_line(line, genome_map, chr_maps):
    # Split the line into non-whitespace parts
    parts = [x for x in line.split(' ') if x != ""]
    name = parts[1].split('.')
    # Complete genome_map and chr_maps
    if len(name) == 1:
        name.append('?')
    elif len(name) > 2:
        name[1] = ".".join(name[1:])
    add_genome(genome_map, chr_maps, name)
    if not name[1] in chr_maps[genome_map[name[0]]]:
        chr_maps[genome_map[name[0]]][name[1]] =\
                    (len(chr_maps[genome_map[name[0]]]), int(parts[5]))
    
    # Extract position, parse sequence
    position = int(parts[2])
    plus_strand = True
    if parts[4] == '-':
        plus_strand = False
        # Conversion of position for '-' strand: not needed
        #position = int(parts[5]) - position - int(parts[3])
    sequence = array.array('c')
    counter = 0
    for letter in parts[6]:
        if letter == '-':
            sequence.append('0')
        else:
            sequence.append('1')
            counter += 1
    
    # Return reference/informant data
    if genome_map[name[0]] == 0:
        return seq.Reference(sequence.tostring(),
            chr_maps[genome_map[name[0]]][name[1]][0], position, plus_strand,
            counter)
    else:
        return genome_map[name[0]], seq.Informant(sequence.tostring(),
            chr_maps[genome_map[name[0]]][name[1]][0], position, plus_strand,
            counter, 0)

# Parses block of alignment
# Returns Reference structure from this block
def parse_block(block_lines, genome_map, chr_maps):
    # If there is something to parse, parse it
    if len(block_lines) <= 0:
        raise BaseException(("no block lines"))
    reference = parse_line(block_lines[0], genome_map, chr_maps)
    for i in range(1,len(block_lines)):
        inf_id, informant = parse_line(block_lines[i], genome_map, chr_maps)
        reference.add_informant(inf_id, informant)
    return reference

def write_bin_block(writer, block, ref_chrs, current_chr):
    block.finalize()
    index_items = []
    for reference in block.ref_blocks:
        index_items.append(seq.IndexItem(reference))
        if reference.chr_id != current_chr:
            current_chr = reference.chr_id
            ref_chrs[current_chr] = 1
        else:
            ref_chrs[current_chr] += 1
    bgzf_starts = bgzftool.write_bgzf_block(writer, block)
    return current_chr, index_items, bgzf_starts

# Parses .maf file, reads it into BinBlock structures
# Returns list of BinBlock structures
def read(filename, genome_map, chr_maps, ref_chrs, writer):
    block_lines = []
    last_bin = current_bin = seq.BinBlock(0)
    next_block = False
    current_chr = -1
    to_index = []
    
    for line in open(filename, 'r'):
        # In case of Windows-like ends of lines:
        line = line.strip() + "\n"
        # Encountered irrelevant line, skip it
        if not (line == "\n" or line[0] == 'a' or line[0] == 's'):
            continue
        # Encountered next block of alignment
        if line[0] == 'a':
            next_block = True
        # Read next line of alignment block, add it to block_lines
        elif line[0] == 's' and next_block:
            block_lines.append(line.strip())
        # Encountered end of alignment block, add it to bin
        elif line == '\n':
            next_block = False
            if len(block_lines) > 0:
                current_bin = current_bin.add(parse_block(block_lines,
                                                          genome_map, chr_maps))
                # If a new bin following the old bin was created
                if last_bin.idn < current_bin.idn:
                    current_chr, index_items, bgzf_starts =\
                        write_bin_block(writer, last_bin, ref_chrs, current_chr)
                    to_index.append((index_items, bgzf_starts))
                    last_bin = current_bin
                # Clear list of lines
                del block_lines[:]
    # Add to bin last alignment block in file
    if len(block_lines) > 0:
        current_bin = current_bin.add(parse_block(block_lines,
                                                  genome_map, chr_maps))
    current_chr, index_items, bgzf_starts =\
        write_bin_block(writer, last_bin, ref_chrs, current_chr)
    to_index.append((index_items, bgzf_starts))
    if last_bin.idn != current_bin.idn:
        current_chr, index_items, bgzf_starts =\
            write_bin_block(writer, current_bin, ref_chrs, current_chr)
        to_index.append((index_items, bgzf_starts))
    bgzftool.close_writer(writer)
    return to_index

# Writes header information from 'genome_map' and 'chr_maps' to 'filestream'
def write_header(filestream, genome_map, chr_maps, ref_chrs, to_index):
    data = array.array("B")
    # Number of genomes
    seq.write_bytes_to(data, len(genome_map), seq.INF_NUM_SIZE)
    # Genome map
    for name in genome_map:
        # Size of name
        seq.write_bytes_to(data, len(name), seq.NAME_SIZE_SIZE)
        # Name itself
        for letter in name:
            data.append(ord(letter))
        # Genome id
        seq.write_bytes_to(data, genome_map[name], seq.INF_ID_SIZE)
    # Chromosome maps:
    for chr_map in chr_maps:
        # Number of chromosomes in this maps
        seq.write_bytes_to(data, len(chr_map), seq.CHR_ID_SIZE)
        for name in chr_map:
            # Size of name
            seq.write_bytes_to(data, len(name), seq.NAME_SIZE_SIZE)
            # Name itself
            for letter in name:
                data.append(ord(letter))
            # Chromosome id
            seq.write_bytes_to(data, chr_map[name][0], seq.CHR_ID_SIZE)
            # Chromosome size
            seq.write_bytes_to(data, chr_map[name][1], seq.LEN_SIZE)
    # Index
    current_chr = -1
    for j in range(len(to_index)):
        for i in range(len(to_index[j][1])):
            if to_index[j][0][i].chr_id != current_chr:
                current_chr = to_index[j][0][i].chr_id
                # Chromosome id and number of references within this chromosome
                seq.write_bytes_to(data, current_chr, seq.CHR_ID_SIZE)
                seq.write_bytes_to(data, ref_chrs[current_chr],
                    seq.REF_COUNT_SIZE)
            # Strand
            seq.write_bytes_to(data, to_index[j][0][i].plus_strand,
                seq.STRAND_SIZE)
            # Position of this reference in chromosome
            seq.write_bytes_to(data, abs(to_index[j][0][i].chr_pos),
                seq.POS_SIZE)
            # Bases count
            seq.write_bytes_to(data, to_index[j][0][i].bases_count,
                seq.POS_SIZE)
            # Pointer
            seq.write_bytes_to(data, to_index[j][1][i],
                seq.FILE_OFFSET_SIZE)
    data.tofile(filestream)
    if filestream != sys.stdout:
        filestream.close()

# Executes all the reading and parsing and writing using 'arguments'
# (as from command line)
def do_all_the_magic(arguments):
    src_file = arguments[1]
    bgzf_file = arguments[2]
    header_file = arguments[3]

    genome_map = {}
    chr_maps = []
    ref_chrs = {}
    writer = bgzftool.get_writer(bgzf_file)
    to_index = read(src_file, genome_map, chr_maps, ref_chrs, writer)

    print "Successfully read", src_file

    print "Compressed data were written into", bgzf_file

    write_header(open(header_file, "wb"), genome_map, chr_maps, ref_chrs,
        to_index)
    print "Binary representation of header was written into", header_file

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: python {0} <from.maf> <to.bgzf> <header.bin>".\
            format(sys.argv[0])
        sys.exit(0)
    
    do_all_the_magic(sys.argv)
