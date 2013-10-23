# -*- coding: utf-8 -*-
# Copyright 2013 by Petra Kubincova

try:
    import Bio.bgzf as bgzf
except ImportError as e:
    print "Missing Biopython module"
    sys.exit(0)
import seq_structures as seq

def write_bgzf(writer, data):
    if len(data) > seq.MAX_LEN:
        raise BaseException()
    else:
        writer.write(data)

def get_writer(filename):
    return bgzf.BgzfWriter(filename)

def close_writer(writer):
    writer.close()

def write_bgzf_block(writer, block):
    data = block.to_bytes()
    bgzf_starts = []
    for reference in data:
        bgzf_starts.append(writer.tell())
        write_bgzf(writer, reference.tostring())
        writer.flush()
    return bgzf_starts

def write_bgzf_blocks(filename, blocks):
    writer = bgzf.BgzfWriter(filename)
    for block in blocks:
        data = block.to_bytes()
        for reference in data:
            block.bgzf_starts.append(writer.tell())
            write_bgzf(writer, reference.tostring())
        writer.flush()
    writer.close()

class ReadBgzf:
    def __init__(self, filename):
        self._filename = filename
        self._reader = bgzf.BgzfReader(filename)
        self._references = {}
    
    def _read_to_number(self, size):
        data = self._reader.read(seq.round_up8(size)/8)
        number = 0
        for value in data:
            if len(value) == 1:
                value = ord(value)
            else:
                value = int(repr(value)[3:-1])
            number <<= 8
            number += int(value)
        return number
    
    def _read_sequence(self, size):
        numbers = []
        size_in_bytes = seq.round_up8(size)/8
        for i in range(size_in_bytes):
            number = self._read_to_number(1)
            binary = bin(number)[2:].zfill(8)
            if size - i*8 >= 8:
                numbers.append(binary)
            else:
                numbers.append(binary[8 - (size - i*8):])
        return ''.join(numbers)
    
    def read_reference(self, virtual_offset):
        if virtual_offset in self._references:
            return self._references[virtual_offset]
        
        self._reader.seek(virtual_offset)
        length = self._read_to_number(seq.LEN_SIZE)
        sequence = self._read_sequence(length)
        inf_num = self._read_to_number(seq.INF_NUM_SIZE)
        reference = seq.Reference(sequence, -2, -2, True)
        # Read informants
        infs = []
        for j in range(inf_num):
            inf_id = self._read_to_number(seq.INF_ID_SIZE)
            inf_block_num = self._read_to_number(seq.BLOCKS_NUM_SIZE)
            infs.append((inf_id, inf_block_num))
        for j in range(len(infs)):
            for k in range(infs[j][1]):
                chr_id = self._read_to_number(seq.CHR_ID_SIZE)
                strand = self._read_to_number(seq.STRAND_SIZE)
                plus_strand = True
                if not strand:
                    plus_strand = False
                chr_pos = self._read_to_number(seq.POS_SIZE)
                seq_pos = self._read_to_number(seq.POS_SIZE)
                seq_pos -= seq.READ_SP_MINUS
                seq_len = self._read_to_number(seq.LEN_SIZE)
                bases_count = self._read_to_number(seq.POS_SIZE)
                sequence = self._read_sequence(seq_len)
                informant = seq.Informant(sequence, chr_id, chr_pos,
                    plus_strand, bases_count, seq_pos)
                reference.add_informant(infs[j][0], informant)
        self._references[virtual_offset] = reference
        return reference
