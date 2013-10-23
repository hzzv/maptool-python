# -*- coding: utf-8 -*-
# Copyright 2013 by Petra Kubincova

import array


MAX_LEN = 64000*8
LEN_SIZE = 32
INF_NUM_SIZE = 8
INF_ID_SIZE = INF_NUM_SIZE
CHR_ID_SIZE = 16
STRAND_SIZE = 8
POS_SIZE = 32 # size of every position and also size of 'bases_count'
BLOCKS_NUM_SIZE = 8
READ_SP_MINUS = 1

NAME_SIZE_SIZE = 8
REF_COUNT_SIZE = 32
FILE_OFFSET_SIZE = 64

BED_NO = 0
BED_PLUS = 1
BED_MINUS = 2

# Generator of bytes in string 'binary', works left-to-right
def byte_generator(binary, size):
    binary = binary.zfill(size)
    for i in range(0, size, 8):
        ret = binary[max(len(binary) - size, 0)+i:i+8]
        yield int(ret, 2)

# Write 'number' in binary form (padded to be 'size' long) to 'data'
def write_bytes_to(data, number, size):
    generator = byte_generator(bin(number)[2:], size)
    for byte in generator:
        data.append(byte)

# "Round up" integer 'number' to nearest multiple of 8 which is bigger
# or equal to 'number'
def round_up8(number):
    return ((number+7)/8)*8


# Block of references written at once in bgzf file
class BinBlock:
    current_size = 0
    blocks_size = 0
    # Needed when reference is too big for 1 block
    ends_in = 0
    
    def __init__(self, idn):
        self.idn = idn
        self.ref_blocks = []
        self.ref_starts = [0]
        self.bgzf_starts = []
        # Whether this block can be manipulated
        self.final = False
        self.current_reference = Reference('', -1, -1, True)
    
    def simple_add(self, reference):
        size = LEN_SIZE + round_up8(len(reference)) + INF_NUM_SIZE +\
            len(reference.informants)*(INF_ID_SIZE + BLOCKS_NUM_SIZE)
        for inf_id in reference.informants:
            size += len(reference.informants[inf_id])*(CHR_ID_SIZE +\
                2*POS_SIZE + LEN_SIZE)
            for informant in reference.informants[inf_id]:
                size += round_up8(len(informant))
        if self.blocks_size + size <= MAX_LEN:
            self.ref_blocks.append(reference)
            self.ref_starts.append(self.ref_starts[-1] + size)
            self.blocks_size += size
            return self
        else:
            if size > MAX_LEN and not self.blocks_size:
                self.ref_blocks.append(reference)
                self.blocks_size = size_increase % MAX_LEN
                self.ends_in = size_increase / MAX_LEN
                return self
            else:
                new_bin_block = BinBlock(self.idn + 1)
                new_bin_block.simple_add(reference)
                finalize(self)
                return new_bin_block
    
    def _undo_reference(self, reference):
        if not self.current_reference:
            raise BaseException()
        if len(self.current_reference) == len(reference):
            print self.current_reference.sequence
            print reference.sequence
            self.current_reference = Reference('', -1, -1, True)
        else:
            self.current_reference.shorten(len(reference))
            for inf_id in reference.informants:
                del self.current_reference.informants[inf_id][-1]
                if len(self.current_reference.informants[inf_id]) == 0:
                    del self.current_reference.informants[inf_id]
    
    def append_current_reference(self):
        self.ref_blocks.append(self.current_reference)
        self.ref_starts.append(self.current_size + self.ref_starts[-1])
        self.blocks_size += self.current_size
        self.current_size = 0
    
    def add(self, reference):
        if self.final:
            raise BaseException()
        
        # First add: initialize 'self.chr_id' and 'self.position'
        if self.current_reference.chr_id == -1:
            self.current_reference = Reference('', reference.chr_id,
                reference.chr_pos, reference.plus_strand, 0)
            self.current_size += LEN_SIZE + INF_NUM_SIZE
        # If 'reference' does not belong into current reference block, add
        # current reference block to 'self.block' and start new reference block
        elif (self.current_reference.chr_id != reference.chr_id or
            self.current_reference.chr_pos + self.current_reference.bases_count\
            != reference.chr_pos or
            self.current_reference.plus_strand != reference.plus_strand):
            self.append_current_reference()
            self.current_reference = Reference('', reference.chr_id,
                reference.chr_pos, reference.plus_strand, 0)
            self.current_size += LEN_SIZE + INF_NUM_SIZE
        
        size_increase = 0
        # For each informant of 'reference': add it to current reference block
        for inf_id in reference.informants:
            for informant in reference.informants[inf_id]:
                if self.current_reference.add_informant(inf_id, informant,
                    len(self.current_reference)):
                    size_increase += INF_ID_SIZE + BLOCKS_NUM_SIZE
                size_increase += CHR_ID_SIZE + 2*POS_SIZE + LEN_SIZE +\
                    round_up8(len(informant))
        
        # To current_reference add sequence of 'reference'
        self.current_reference.add_sequence(Sequence.__str__(reference, False),
            reference.bases_count)
        size_increase += round_up8(len(reference))
        
        # If with 'reference' is whole BinBlock too big, undo adding 'reference'
        if (self.blocks_size + self.current_size + size_increase > MAX_LEN and
                (self.current_size or self.blocks_size)):
            self._undo_reference(reference)
            new_bin_block = BinBlock(self.idn+1)
            new_bin_block.add(reference)
            self.finalize()
            return new_bin_block
        # If with 'reference' is whole BinBlock too big, but it is the only
        # reference in BinBlock, split it into more BinBlock-s
        elif size_increase > MAX_LEN:
            self.ref_blocks.append(current_reference)
            self.blocks_size = size_increase % MAX_LEN
            self.ends_in = size_increase / MAX_LEN
            return self
        else:
            self.current_size += size_increase
            return self
        
    def finalize(self):
        if not self.final:
            if len(self.current_reference):
                self.append_current_reference()
            self.final = True
        return self
    
    def __str__(self):
        ret = "\n----- Block no. " + str(self.idn) + " with " +\
            str(len(self.ref_blocks)) + " references and overall size " +\
            str(self.blocks_size) + " ending in " + str(self.ends_in) + ":"
        # For every reference in block
        for reference in self.ref_blocks:
            ret += "\n" + str(reference) + "\n"
            # And for every informant in that reference
            for inf_id in reference.informants:
                ret += "Informant id: " + str(inf_id) + "\n"
                # And for every piece of alignment in that informant
                for piece in reference.informants[inf_id]:
                    ret += str(piece) + "\n"
        ret += "----- End of block\n"
        return ret
    
    def to_bytes(self):
        data = []
        for reference in self.ref_blocks:
            data.append(array.array("B"))
            data[-1].extend(reference.to_bytes())
            for inf_id in reference.informants:
                id_generator = byte_generator(bin(inf_id)[2:], INF_ID_SIZE)
                num_generator =\
                    byte_generator(bin(len(reference.informants[inf_id]))[2:],
                    BLOCKS_NUM_SIZE)
                for byte in id_generator:
                    data[-1].append(byte)
                for byte in num_generator:
                    data[-1].append(byte)
            for inf_id in reference.informants:
                for informant in reference.informants[inf_id]:
                    data[-1].extend(informant.to_bytes())
            while len(data[-1]) > MAX_LEN:
                data.append(data[-1][MAX_LEN:])
                data[-2] = data[-2][:MAX_LEN]
        return data


# Parent class for Reference and Informant, stores sequence and some
# additional information
class Sequence:
    def __init__(self, sequence, chr_id, chr_pos, plus_strand, bases_count):
        self.more_sequences = False
        self.sequence = sequence
        self.chr_id = chr_id
        self.chr_pos = chr_pos
        self.plus_strand = plus_strand
        self.bases_count = bases_count
    
    def add_sequence(self, sequence):
        if not self.more_sequences:
            self.sequence = [self.sequence]
            self.more_sequences = True
        self.sequence.append(sequence)
    
    def shorten(self, by):
        to_size = len(self) - by
        if not self.more_sequences:
            self.sequence = self.sequence[:to_size]
        else:
            for i in range(len(self.sequence)):
                if len(self.sequence[i]) <= to_size:
                    to_size -= len(self.sequence[i])
                elif to_size > 0:
                    self.sequence[i] = self.sequence[i][:to_size]
                    to_size = 0
                else:
                    self.sequence = self.sequence[:i]
    
    def to_bytes(self):
        # Get reprezentation of sequence in this Sequence in byte array
        binary = Sequence.__str__(self)
        generator = byte_generator(binary, len(binary))
        data = array.array("B")
        for byte in generator:
            data.append(byte)
        return data
    
    def get_sequence(self):
        if not self.more_sequences:
            return self.sequence
        else:
            return "".join(self.sequence)
    
    def __str__(self, pad=False):
        if not self.more_sequences:
            ret = self.sequence
        else:
            ret = ""
            for sequence in self.sequence:
                ret += sequence
            self.sequence = ret
            self.more_sequences = False
        if pad:
            ret += '0' * (8 - (len(ret) % 8))
        return ret
    
    def __len__(self):
        if self.more_sequences:
            ret = 0
            for sequence in self.sequence:
                ret += len(sequence)
            return ret
        else:
            return len(self.sequence)


class Reference(Sequence):
    def __init__(self, sequence, chr_id, chr_pos, plus_strand, bases_count=-1):
        Sequence.__init__(self, sequence, chr_id, chr_pos, plus_strand,
            bases_count)
        self.informants = {}
    
    def add_informant(self, inf_id, informant, seq_pos=-1):
        if not inf_id in self.informants:
            self.informants[inf_id] = []
            ret = True
        else:
            ret = False
        if seq_pos != -1:
            informant.seq_pos = seq_pos
        self.informants[inf_id].append(informant)
        return ret
    
    def add_sequence(self, sequence, bases_count):
        self.bases_count += bases_count
        Sequence.add_sequence(self, sequence)
    
    def to_bytes(self):
        # Create reprezentation of this Reference in byte array
        data = array.array("B")
        length_generator = byte_generator(bin(len(self))[2:], LEN_SIZE)
        inf_num_generator = byte_generator(bin(len(self.informants))[2:],
                                           INF_NUM_SIZE)
        for byte in length_generator:
            data.append(byte)
        data.extend(Sequence.to_bytes(self))
        for byte in inf_num_generator:
            data.append(byte)
        return data
    
    def __str__(self):
        if self.plus_strand:
            strand = '+'
        else:
            strand = '-'
        return "Reference chromosome: " + str(self.chr_id) + ", position: " +\
            str(self.chr_pos) + " on " + strand + " strand, seq: " +\
            Sequence.__str__(self, False) +\
            " of length " + str(Sequence.__len__(self)) + " with " +\
            str(self.bases_count) + " bases and with " +\
            str(len(self.informants)) + " informants"


class Informant(Sequence):
    def __init__(self, sequence, chr_id, chr_pos, plus_strand, bases_count,
        seq_pos=-1):
        Sequence.__init__(self, sequence, chr_id, chr_pos, plus_strand,
            bases_count)
        self.seq_pos = seq_pos
    
    def to_bytes(self):
        # Create reprezentation of this Informant in byte array
        data = array.array("B")
        generator = []
        generator.append(byte_generator(bin(self.chr_id)[2:], CHR_ID_SIZE))
        generator.append([self.plus_strand])
        generator.append(byte_generator(bin(abs(self.chr_pos))[2:], POS_SIZE))
        generator.append(byte_generator(bin(self.seq_pos + READ_SP_MINUS)[2:],
            POS_SIZE))
        generator.append(byte_generator(bin(len(self.sequence))[2:], LEN_SIZE))
        generator.append(byte_generator(bin(self.bases_count)[2:], POS_SIZE))
        for g in generator:
            for byte in g:
                data.append(byte)
        data.extend(Sequence.to_bytes(self))
        return data
    
    def __str__(self):
        strand = '+'
        if not self.plus_strand:
            strand = '-'
        return "Chromosome id: " + str(self.chr_id) + ", sequence: " +\
            Sequence.__str__(self, False) + " of length " +\
            str(len(self.sequence)) + " on position " + str(self.chr_pos) +\
            " on " + strand + " strand, " +\
            " in alignment on position " + str(self.seq_pos) + " with " +\
            str(self.bases_count) + " bases"


# Class containing information about mapping query
class Mapping:
    def __init__(self, args, informant=False, options=False, bed=False):
        if bed:
            self.chromosome = args[0]
            self.bed = (BED_PLUS, BED_MINUS)[args[1] == "-"]
            self.strand = True
            self.left, self.right = int(args[2]), int(args[3])
            self.inf_name = informant
            self.inf_maxgap = self.ref_maxgap = int(options[0])
            self.inner = (True, False)[options[1] == "-outer"]
            self.all_errors = (False, True)[options[2] == "-allerrors"]
            self.always_map = (False, True)[options[2] == "-alwaysmap" or\
                options[2] == "-allerrors"]

class BedQuery:
    def __init__(self, args, line_number=-1):
        self.length = len(args)
        # Defaults: name, score, strand, thick start and end, RGB, exon count,
        # exons sizes and starts
        defaults = [str(line_number), "-1", "+", "-1", "-1", "0", "0", "", ""]
        args += defaults[(len(args) - 3):]
        self.chromosome, self.left, self.right, self.name, self.score,\
            self.strand, self.lthick, self.rthick, self.rgb, self.exon_count,\
            self.sizes, self.starts = args[0:2] + [str(int(args[2]) - 1)] +\
            args[3:7] + [str(int(args[7]) - 1), args[8], int(args[9]), \
            [int(x) for x in str(args[10]).split(",") if x != ""], \
            [int(args[1]) + int(x)\
                for x in str(args[11]).split(",") if x != ""]]
        if len(self.sizes) != len(self.starts) or\
            len(self.sizes) != self.exon_count:
            raise Exception("Invalid input: indicated number of blocks " +\
                "in line '" + " ".join(args) + "' does not match the actual " +\
                "number of blocks.")

    def get_map(self):
        return [self.chromosome, self.strand, self.left, self.right]
    
    def get_thick_map(self):
        return [self.chromosome, self.strand, self.lthick, self.rthick]
        
    def get_left_thick(self):
        return [self.chromosome, self.strand, self.lthick, self.lthick]
    
    def get_right_thick(self):
        return [self.chromosome, self.strand, self.rthick, self.rthick]

    def get_exon_maps(self):
        return [[self.chromosome, self.strand, str(self.starts[i]),
            str(self.sizes[i] + self.starts[i] - 1)]
            for i in range(0, self.exon_count)]

class IndexItem:
    def __init__(self, reference):
        self.chr_id = reference.chr_id
        self.chr_pos = reference.chr_pos
        self.plus_strand = reference.plus_strand
        self.bases_count = reference.bases_count