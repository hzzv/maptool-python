# -*- coding: utf-8 -*-
# Copyright 2013 by Petra Kubincova

import sys
import read_header_file
import bgzf_tool as bgzftool
import seq_structures as seq


messages = {"inf_preceed":\
    "In informant: one sequence does not preceed the next one",\
    "inf_chr":\
    "In informant: sequences are from different source chromosomes",\
    "inf_strand":\
    "In informant: sequences are from different strands",\
    "inf_overlap":\
    "In informant: two sequences overlap",\
    "inf_gap":\
    "In informant: there are gaps of width ",\
    "ref_strand":\
    "In reference: sequences are from different strands",\
    "ref_left":\
    "Reference does not contain left ending point of the interval",\
    "ref_gap":\
    "In reference: there are gaps of width ",\
    "ref_cont":\
    "Interval is not contained whole in reference sequence",\
    "pos_to_gap":\
    "Position maps to gap",\
    "int_to_gap":\
    "Interval maps to gap",\
    "discontinuous_int":\
    "Interval maps to discontinuous intervals",\
    "invalid_int":\
    "Input interval is not valid",\
    "invalid_chr":\
    "Invalid chromosome name",\
    "incorrect_int":\
    "Interval maps to incorrect intervals",
    "thick":\
    "Mapping of the thick interval is incorrect",\
    "no_mapping":\
    "There is no mapping of the interval",\
    "inner":\
    "One or both of the ending points were mapped in the inner way",\
    "outer":\
    "One or both of the ending points were mapped in the outer way"}

terminal_errors = {"ref_strand": False, "ref_left": False,\
    "ref_gap": False, "ref_cont": False, "discontinuous_int": False,\
    "invalid_int": False, "invalid_chr": False, "incorect_int": False,\
    "thick": False, "no_mapping": False}
errors = {"inf_preceed": False, "inf_chr": False, "inf_strand": False,\
    "inf_overlap": False, "inf_gap": False, "pos_to_gap": False,\
    "int_to_gap": False}
warnings = {"inner": False, "outer": False}
inf_gaps = []
ref_gaps = []
all_errors_count = terminal_errors_count = 0

# Clears all errors
def clear_errors():
    global all_errors_count, terminal_errors_count 
    for key in terminal_errors:
        terminal_errors[key] = False
    for key in errors:
        errors[key] = False
    for key in warnings:
        warnings[key] = False
    inf_gaps = []
    ref_gaps = []
    all_errors_count = terminal_errors_count = 0

# Adds error with id 'error' to current errors ('terminal_errors'
# or 'erorrs' or 'warnings')
def add_error(error, gap=0):
    global all_errors_count, terminal_errors_count
    all_errors_count += 1
    if error in terminal_errors:
        terminal_errors[error] = True
        terminal_errors_count += 1
    elif error in errors:
        errors[error] = True
    elif error in warnings:
        warnings[error] = True
    if error == "ref_gap" and gap > 0:
        ref_gaps.append(str(gap))
    elif error == "inf_gap" and gap > 0:
        inf_gaps.append(str(gap))

# Prints all errors, then clears them
def print_errors(name, always_map):
    no_t_errors = no_errors = True
    for key in terminal_errors:
        if terminal_errors[key]:
            no_t_errors = False
            ending = ""
            if key == "ref_gap":
                ending = "\t" + ",".join(ref_gaps)
            print >> sys.stderr, name + "\tno\t" + key + "\t" +\
                messages[key] + ending
    for key in errors:
        if errors[key]:
            no_errors = False
            ending = ""
            if key == "inf_gap":
                ending = "\t" + ",".join(inf_gaps)
            if no_t_errors and always_map:
                print >> sys.stderr, name + "\tw/error\t" + key + "\t" +\
                    messages[key] + ending
            else:
                print >> sys.stderr, name + "\tno\t" + key + "\t" +\
                    messages[key] + ending
    for key in warnings:
        if warnings[key]:
            no_errors = False
            if no_t_errors and always_map:
                print >> sys.stderr, name + "\tw/error\t" + key + "\t" +\
                    messages[key]
            else:
                print >> sys.stderr, name + "\tno\t" + key + "\t" +\
                    messages[key]
    if no_t_errors and no_errors:
        print >> sys.stderr, name + "\t" + "mapped"
    clear_errors()

# When 'chr_map' maps chromosome name to its [id, length], this function
# returns map which maps chromosome id to its [name, length]
def inverse_to(chr_map):
    inverse = {}
    for key in chr_map:
        inverse[chr_map[key][0]] = [key, chr_map[key][1]]
    return inverse

# Binary search: returns tuple from 'reference_index' and index of it. Tuple
# represents chunk of sequence of reference genome which contains position
# 'number', or, if 'number' is not contained in any chunk of sequence, it
# returns the closest one to the right (and index of it)
def get_reference(reference_index, number, where_from=0):
    if number < reference_index[where_from][0]:
        return (reference_index[where_from], where_from)
    elif number >= reference_index[-1][0] + reference_index[-1][1]:
        return (None, -1)
    lo, hi = where_from, len(reference_index)
    while True:
        middle = (lo+hi)/2
        if (number >= reference_index[middle][0] and
            number < reference_index[middle][0] + reference_index[middle][1]) \
            or (middle != 0 and number < reference_index[middle][0] and
            number >=\
            reference_index[middle-1][0] + reference_index[middle-1][1]):
            return (reference_index[middle], middle)
        elif middle == lo or lo >= hi:
            return None
        elif number < reference_index[middle][0]:
            hi = middle
        elif number >= reference_index[middle][0] + reference_index[middle][1]:
            lo = middle

# Are Informant-s 'inf1', 'inf2' following each other?
# Returns [bool, list], where 'bool' is True if 'inf2' is following 'inf1' close
# enough, on the right strand and in the right chromosome, and 'list' is list
# of conditions that were not met
def following_informants(inf1, inf2, allowed_gap):
    chromosome = inf1.chr_id == inf2.chr_id
    strand = inf1.plus_strand == inf2.plus_strand
    overlap = inf1.chr_pos + inf1.bases_count <= inf2.chr_pos or\
        inf2.chr_pos + inf2.bases_count <= inf1.chr_pos
    follows = not strand or inf1.chr_pos < inf2.chr_pos
    gap = allowed_gap == -1 or\
        inf1.chr_pos + inf1.bases_count + allowed_gap >= inf2.chr_pos
    ret = []
    if not follows:
        ret.append(["inf_preceed"])
    if not chromosome:
        ret.append(["inf_chr"])
    if not strand:
        ret.append(["inf_strand"])
    if not overlap:
        ret.append(["inf_overlap"])
    if not gap:
        ret.append(["inf_gap", inf2.chr_pos - inf1.chr_pos - inf1.bases_count])
    return [follows and chromosome and strand and overlap and gap, ret]

# Try to find interval in index of references, read all the references
# 'ref_chr_id' is integer id of reference chromosome from which is interval
# 'mapping' is Mapping structure containing information about desired mapping
# 'index' si list of tuples containing information about Reference-s written
# in 'source'
# 'inf_id' is integer id of informant to which is interval mapped
# 'reader' is ReadBgzf instance
# Returns list of references containing the interval, list of occured errors
# and list of continuous intervals
def get_references(ref_chr_id, mapping, index, inf_id, reader):
    left, right = mapping.left, mapping.right
    ref_information = []
    references = []
    i = 0
    last_informant = None
    inf_errors = []
    continuous = [] # List of continuous intervals of informants
    while left <= right:
        ref_info, i = get_reference(index[ref_chr_id], left, i)
        # Is the reference sequence on correct strand?
        if ref_info != None and ref_info[3] != mapping.strand:
            # If no, search for another
            continue
        if ref_info != None and len(ref_information) > 0 and\
            ref_information[-1] != None and\
            ref_info[3] != ref_information[-1][3]:
            add_error("ref_strand")
            if not mapping.always_map:
                sys.exit(0)
        # Is there gap at most 'reference_gap' long?
        if ref_info == None and right - left <= mapping.ref_maxgap:
            break
        # Now the 'ref_info' can be added to already collected 'ref_information'
        ref_information.append(ref_info)
        if ref_information[-1] == None and right - left > mapping.ref_maxgap:
            if left == mapping.left:
                add_error("ref_left")
                if mapping.always_map:
                    left = right
                    continue
                else:
                    sys.exit(0)
            else:
                if mapping.always_map:
                    add_error("ref_gap", right-left)
                    return (references, inf_errors, continuous)
                else:
                    add_error("ref_cont")
                    sys.exit(0)
        if ref_information[-1] != None and left < ref_information[-1][0] \
            and ref_information[-1][0] - left > mapping.ref_maxgap:
            if mapping.left == left:
                if right < ref_information[-1][0] and\
                    not mapping.always_map:
                    add_error("ref_cont")
                    sys.exit(0)
                else:
                    add_error("ref_left")
                    if not mapping.always_map:
                        sys.exit(0)
            else:
                add_error("ref_gap", ref_information[-1][0]-left)
                if not mapping.always_map:
                    sys.exit(0)
        # Read reference corresponding to ref_info from .bgzf file
        references.append(reader.read_reference(ref_info[2]))
        references[-1].bases_count = ref_info[1]
        references[-1].chr_id = ref_chr_id
        references[-1].chr_pos = ref_info[0]
        references[-1].plus_strand = ref_info[3]
        if not inf_id in references[-1].informants:
            references[-1].informants[inf_id] = ()
        if len(continuous) == 0 and len(references[-1].informants[inf_id]) > 0:
            continuous.append([references[-1].informants[inf_id][0],
                references[-1].informants[inf_id][0]])
        # Read through chunks of informant sequence aligned to this reference:
        # Can the first informant of this reference be following last informant
        # of the previous reference?
        if last_informant != None and \
            len(references[-1].informants[inf_id]) > 0:
            ret = following_informants(references[-2].informants[inf_id][-1],
            references[-1].informants[inf_id][0], mapping.inf_maxgap)
            if not ret[0]:
                inf_errors.append((references[-2].informants[inf_id][-1],
                    references[-1].informants[inf_id][0], ret[1]))
                continuous[-1][1] = references[-2].informants[inf_id][-1]
                continuous.append([references[-1].informants[inf_id][0], 0])
        # For every 'j': can the 'j'-th informant be following the 'j-1'-th?
        for j in range(1, len(references[-1].informants[inf_id])):
            ret = following_informants(references[-1].
                informants[inf_id][j-1],
                references[-1].informants[inf_id][j], mapping.inf_maxgap)
            if not ret[0]:
                inf_errors.append((references[-1].informants[inf_id][j-1],
                    references[-1].informants[inf_id][j], ret[1]))
                continuous[-1][1] = references[-1].informants[inf_id][j-1]
                continuous.append([references[-1].informants[inf_id][j],
                    0])
        # Shift: left ending point ('left') of reference is now position just
        # after the end of the reference found in this iteration
        # And the 'last_informant' is now the rightmost informant
        # of the reference found in this iteration
        left = ref_information[-1][0] + ref_information[-1][1]
        if len(references[-1].informants[inf_id]) > 0:
            last_informant = references[-1].informants[inf_id][-1]
            continuous[-1][1] = references[-1].informants[inf_id][-1]
    return (references, inf_errors, continuous)


# Find the informant aligned to position 'position' or if there is no such,
# find the closest one (to the right if 'to_the_right' is True, to the left
# if not), find the aligned (or closest-to-aligned) position in it
# 'references' is list of references containing the position
# 'position' is integer position to align
# 'ref_begin' is index to 'references' where to start search
# 'to_the_right' is boolean - True when aligned position is to be seeked
# to the right, False otherwise
# Return tuple: index of the referenece in 'ref_index', index of the informant
# in 'inf_index', and position in the informant in 'inf_position'
def find_aligned_position(references, position, ref_begin, to_the_right, inf_id):
    inf_index = -1
    inf_position = -1
    ref_index = 0
    if to_the_right:
        # Search through all references (from 'ref_begin' to the right)
        # and through all informants belonging to them
        for r in range(ref_begin, len(references)):
            for j in range(0, len(references[r].informants[inf_id])):
                # Is the end of the informant sequence beyond 'position'?
                if position < references[r].informants[inf_id][j].seq_pos +\
                    len(references[r].informants[inf_id][j].get_sequence()):
                    inf_index = j
                    ref_index = r
                    # Is the beginning of the informant sequence before
                    # 'position'?
                    if position >= references[r].informants[inf_id][j].seq_pos:
                        inf_position = position -\
                            references[r].informants[inf_id][j].seq_pos
                    else:
                        inf_position = 0
                    break
            if inf_index != -1:
                break
        if inf_index == -1:
            return (-1, -1, -1)
    else:
        # Search through all references (from 'ref_begin' to the left)
        # and through all informants belonging to them
        for r in range(ref_begin, -1, -1):
            for j in range(0, len(references[r].informants[inf_id])):
                # Is the end of the informant sequence beyond 'position'?
                if position < references[r].informants[inf_id][j].seq_pos +\
                    len(references[r].informants[inf_id][j].get_sequence()):
                    ref_index = r
                    # Is the beginning of the informant sequence before
                    # 'position'?
                    if position >= references[r].informants[inf_id][j].seq_pos:
                        inf_index = j
                        inf_position = position -\
                            references[r].informants[inf_id][j].seq_pos
                    elif j > 0:
                        inf_index = j-1
                        inf_position =\
                            len(references[r].informants[inf_id][j-1].
                            get_sequence()) - 1
                    break
                # If all informant sequences were before 'position'
                elif j == len(references[r].informants[inf_id]) - 1:
                    ref_index = r
                    inf_index = j
                    inf_position = len(references[r].informants[inf_id][j].
                        get_sequence()) - 1
            if inf_index != -1:
                break
        if inf_index == -1:
            return (-1, -1, -1)
    return (ref_index, inf_index, inf_position)

# Now find base in informant closest to 'inf_position', but aligned
# to some base in reference
# 'references' is list of references containing the position
# 'ref_index' is integer index to 'references' where to start the search
# 'inf_index' is integer index to 'references[ref_index].informants[inf_id]'
# where to start the search
# 'inf_position' is integer index to sequence of informant where to start the
# search
# 'inf_id' is id of informant which is to be searched
# 'to_the_right' is boolean - True when aligned position is to be seeked
# to the right, False otherwise
# Return tuple of altered (ref_index, inf_index, inf_position)
def find_aligned_base(references, ref_index, inf_index, inf_position, inf_id,
    to_the_right):
    difference = [-1, 1][to_the_right]
    ref_sequence = references[ref_index].get_sequence()
    while True:
        inf_sequence =\
            references[ref_index].informants[inf_id][inf_index].get_sequence()
        # Iterate through 'inf_sequence'
        while len(inf_sequence) > inf_position >= 0:
            if inf_sequence[inf_position] == '1' \
                and ref_sequence[\
                references[ref_index].informants[inf_id][inf_index].seq_pos +\
                inf_position] == '1':
                return (ref_index, inf_index, inf_position)
            inf_position += difference
        # When the end of 'inf_sequence' is reached, there are 3 possibilities:
        # Are there more informants aligned to this reference?
        if 0 <= inf_index + difference <\
            len(references[ref_index].informants[inf_id]):
            inf_index += difference
            inf_position = 0
            if not to_the_right:
                inf_position = len(references[ref_index].\
                    informants[inf_id][inf_index].get_sequence()) - 1
        # Are there more references?
        elif 0 <= ref_index + difference < len(references):
            ref_index += difference
            ref_sequence = references[ref_index].get_sequence()
            inf_index = 0
            if not to_the_right:
                inf_index = len(references[ref_index].informants[inf_id]) - 1
            inf_position = 0
            if not to_the_right:
                inf_position =\
                    len(references[ref_index].informants[inf_id][inf_index].
                    get_sequence()) - 1
        # We reached the end of whole sequence: there is no base on which
        # can be the 'position' mapped
        else:
            inf_position = -1
            break
    return (ref_index, inf_index, inf_position)

# Return genome position of 'inf_position' on 'informants' sequence,
# its chromosome name (using 'id_map'), strand ('+' or '-')
# and informant, if 'inf_position' == -1 return last position
# does not convert from '+' coordinates to '-'
def get_pos_chr_strand_inf(informant, inf_position, id_map):
    if inf_position != -1:
        counter = 0
        inf_sequence = informant.get_sequence()
        for i in range(0, inf_position):
            if inf_sequence[i] == '1':
                counter += 1
        position = counter + informant.chr_pos
    else:
        position = informant.chr_pos + informant.bases_count - 1
    # Is the position on + or - strand?
    strand = ['-', '+'][informant.plus_strand]
    # Conversion of position for '-' strand: not needed
    #position = id_map[informant.chr_id][1] - position
    return (position, id_map[informant.chr_id][0], strand, informant)


# Mapping: linear algorithm
# 'references' is list of all references containing sequence from the interval
# 'position' is integer position to be mapped
# 'to_the_right' is boolean - whether to map to the right or to the left
# 'in_this_reference' is integer index of reference in 'references' where
# is the 'position' contained
# 'inf_id' is integer id of informant to which is the position mapped
# Returns mapped position in informant genome or -1 in case of impossible
# mapping
def map_position(references, position, to_the_right, in_this_reference, inf_id,
    id_map):
    # Position in sequence = position in whole genome - start of sequence
    in_sequence = position - references[in_this_reference].chr_pos
    if in_sequence < 0:
        in_sequence = 0
    elif in_sequence >= references[in_this_reference].chr_pos +\
        references[in_this_reference].bases_count:
        in_sequence = references[in_this_reference].bases_count - 1
    # Find index of 'in_sequence'-th 1, store it in 'i'
    ref_sequence = references[in_this_reference].get_sequence()
    for i in range(0, len(ref_sequence)):
        if in_sequence == 0 and ref_sequence[i] == '1':
            break
        elif ref_sequence[i] == '1':
            in_sequence -= 1
    # Get position in informant aligned to 'i', get index of the informant
    # and index of the reference where it occurs
    ref_index, inf_index, inf_position = find_aligned_position(references,
        i, in_this_reference, to_the_right, inf_id)
    if ref_index == -1:
        return -1
    # Get closest base to 'inf_position' which is aligned to some base
    # in reference, store it in inf_position
    ref_index, inf_index, inf_position = find_aligned_base(references,
        ref_index, inf_index, inf_position, inf_id, to_the_right)
    if inf_position == -1:
        return -1
    return get_pos_chr_strand_inf(references[ref_index].informants[inf_id]\
        [inf_index], inf_position, id_map)

# Makes one interval in format [chr_name, strand, start, end] from left
# and right ending points ('start' and 'end') in format
# [position, chr_name, strand, informant] (output format of
# get_pos_chr_strand_inf())
def make_interval(start, end):
    if start[2] == True:
        start[2] = '+'
    elif start[2] == False:
        start[2] = '-'
    return [start[1], start[2], start[0], end[0]]

# Read header, get name of reference genome. Will be needed in any case
def get_header_info(filename):
    genome_map, chr_maps, index = read_header_file.read(open(filename, "rb"))
    ref_name = ""
    for key in genome_map:
        if genome_map[key] == 0:
            ref_name = key
    return (genome_map, chr_maps, index, ref_name)

# Print informants contained in 'genome_map'
def print_informants(genome_map, ref_name):
    print "Informants aligned to {0}:".format(ref_name)
    for key in genome_map:
        if key != ref_name:
            print key

# Print chromosomes of reference genome
def print_chromosomes(chr_maps, ref_name):
    print "Chromosomes of reference genome ({0}):".format(ref_name)
    for key in chr_maps[0]:
        print key


# If this is mapping of position (= interval of 1 base),
# check if it is correct
# Uses integer positions 'left' and 'right', output of map_position() 'left_map'
# and 'right_map', flags 'inner' and 'always_map'
# Returns information if line "Mapping might be incorrect..." was printed
def interval_of_one(left, right, left_map, right_map, inner, always_map):
    if left == right and (left_map == -1 or right_map == -1 or\
        left_map[0] != right_map[0] or left_map[1] != right_map[1] or\
        left_map[2] != left_map[2]):
        if inner and not always_map:
            add_error("pos_to_gap")
            sys.exit(0)
        else:
            add_error("pos_to_gap")
    if left_map != -1 and right_map != -1 and\
        left_map[1] == right_map[1] and left_map[2] == right_map[2] and\
        left_map[0] > right_map[0]:
        if inner and not always_map:
            add_error("int_to_gap")
            sys.exit(0)
        else:
            add_error("int_to_gap")

# Check if 'left_map' and 'right_map' (outputs of map_position()) are correct,
# if not (and 'always_map' is True), try to map 'mapping.left'
# or 'mapping.right' the other way (if 'mapping.inner' is True, then outer,
# else inner way)
# Uses 'references': list of seq.Reference-s;  information from Mapping
# structure 'mapping'; id of informant 'inf_id'; 'id_map'
# Returns not necessarily altered 'left_map' and 'right_map'
def map_the_other_way(left_map, right_map, mapping, references, inf_id, id_map):
    if left_map == -1 or right_map == -1:
        way = ("outer", "inner")
        # Mapping of the left ending point
        if mapping.always_map and left_map == -1:
            left_map = map_position(references, mapping.left, not mapping.inner,
                0, inf_id, id_map)
            add_error(way[not mapping.inner])
        # Mapping of the right ending point
        if mapping.always_map and right_map == -1:
            right_map = map_position(references, mapping.right, mapping.inner,
                len(references)-1, inf_id, id_map)
            add_error(way[not mapping.inner])
        # If mapping is not forced or is not possible even if forced, quit
        if not mapping.always_map or left_map == -1 or right_map == -1:
            add_error("no_mapping")
            sys.exit(0)
    return (left_map, right_map)

# Error handling after mapping (when 'left_map' and 'right_map' are correct),
# checks if there are errors between relevant informant sequences.
# Computes indices of first and last relevant intervals in 'continuous'.
# Uses 'left_map' and 'right_map': output of map_position(); 'continuous': list
# of continuous intervals of informant sequence; 'inf_errors': list of errors in
# informant genome; flag 'always_map'
# Returns 'first' and 'last' - indices of first and last relevant interval
def mapping_errors(left_map, right_map, continuous, inf_errors, always_map,
    chr_map):
    global all_errors_count, terminal_errors_count
    if always_map or terminal_errors_count == 0:
        # find first and last relevant interval of informants,
        first_end = last_start = None
        first = last = -1
        for i in range(0, len(continuous)):
            if continuous[i][0].chr_pos <= left_map[0] and \
                left_map[0] < continuous[i][1].chr_pos +\
                continuous[i][1].bases_count and\
                ((left_map[2] == '+' and continuous[i][0].plus_strand) or\
                (left_map[2] == '-' and not continuous[i][0].plus_strand)) and\
                continuous[i][0].chr_id == chr_map[left_map[1]][0]:
                first = i
                first_end = continuous[i][1]
                if last != -1:
                    break
            if continuous[i][0].chr_pos <= right_map[0] and \
                right_map[0] < continuous[i][1].chr_pos +\
                continuous[i][1].bases_count and\
                ((right_map[2] == '+' and continuous[i][0].plus_strand) or\
                (right_map[2] == '-' and not continuous[i][0].plus_strand)) and\
                continuous[i][0].chr_id == chr_map[right_map[1]][0]:
                last = i
                last_start = continuous[i][0]
                if first != -1:
                    break
        if first == -1 or first > last:
            if always_map:
                add_error("int_to_gap")
                if first == -1:
                    first = last
            else:
                add_error("int_to_gap")
                sys.exit(0)
        # then search through errors and save only those in relevant intervals,
        relevant_errors = False
        if first != last:
            for e in inf_errors:
                if e[0] == first_end:
                    relevant_errors = True
                if relevant_errors:
                    for inf_e in e[2]:
                        if len(inf_e) == 1:
                            add_error(inf_e[0])
                        else:
                            add_error(inf_e[0], inf_e[1])
                if e[1] == last_start:
                    break
    if (not always_map and all_errors_count > 0) or\
        terminal_errors_count > 0:
        sys.exit(0)
    return (first, last)


# Checks last two intervals in list 'intervals' if they are following
def check_following_intervals(intervals):
    if intervals[-1] == None or len(intervals) < 2:
        return True
    if intervals[-1][0] != intervals[-2][0] or\
        intervals[-1][1] != intervals[-2][1]:
        return False
    else:
        return True


# Checks all pairs of following exons if they really are following
# one another
def check_exons(exons):
    for i in range(0, len(exons)-1):
        if exons[i][0] != exons[i+1][0] or\
            exons[i][1] != exons[i+1][1] or\
            exons[i][3] >= exons[i+1][2]:
                return False
    return True


# Returns all relevant intervals from 'continuous' (beginning with 'first',
# ending with 'last') regarding 'left_map' and 'right_map' (mapping of leftmost
# and rightmost positions), using 'id_map'
def get_mapped_intervals(left_map, right_map, first, last, continuous,
    id_map):
    intervals = []
    if first < last:
        intervals.append(make_interval(left_map,
            get_pos_chr_strand_inf(continuous[first][1], -1, id_map)))
        if not check_following_intervals(intervals):
            return -1
    for i in range(first+1, last):
        intervals.append(make_interval(
            get_pos_chr_strand_inf(continuous[i][0], 0, id_map),
            get_pos_chr_strand_inf(continuous[i][1],-1, id_map)))
        if not check_following_intervals(intervals):
            return -1
    if first < last:
        intervals.append(make_interval(
            get_pos_chr_strand_inf(continuous[last][0], 0, id_map),
            right_map))
        if not check_following_intervals(intervals):
            return -1
    elif first == last:
        intervals.append(make_interval(left_map, right_map))
        if not check_following_intervals(intervals):
            return -1
    else:
        intervals.append(make_interval(left_map, left_map))
        intervals.append(make_interval(right_map, right_map))
        if not check_following_intervals(intervals):
            return -1
    return intervals

# Prints BED format of intervals (in one line)
def print_bed(intervals, thick_interval, bed, chr_map):
    if intervals == -1 or len(intervals) == 0 or thick_interval == -1 or\
        len(thick_interval) > 1 or len(thick_interval) == 0:
        add_error("discontinuous_int")
        return
    if intervals[0][1] == "-":
        for interval in intervals:
            interval[2], interval[3] =\
                chr_map[interval[0]][1] - interval[3] - 1,\
                chr_map[interval[0]][1] - interval[2] - 1
        intervals.reverse()
    if thick_interval[0][1] == "-":
        if thick_interval[0][2] != -1 and thick_interval[0][3] != -1:
            thick_interval[0][2], thick_interval[0][3] =\
                chr_map[thick_interval[0][0]][1] - thick_interval[0][3] - 1,\
                chr_map[thick_interval[0][0]][1] - thick_interval[0][2] - 1
    if (bed.strand == '+' and intervals[0][1] == '+') or\
        (bed.strand == '-' and intervals[0][1] == '-'):
        strand = '+'
    else:
        strand = '-'
    to_print = intervals[0][0] + "\t" + str(intervals[0][2]) +\
        "\t" + str(intervals[-1][3] + 1)
    if bed.length > 3:
        to_print += "\t" + bed.name
    if bed.length > 4:
        to_print += "\t" + str(bed.score)
    if bed.length > 5:
        to_print += "\t" + strand
    if thick_interval[0][2] != -1 or bed.length > 7 or bed.exon_count != 0:
        to_print += "\t" + str(thick_interval[0][2]) + "\t" +\
            str(thick_interval[0][3] + 1)
    if bed.length > 8 or bed.exon_count != 0:
        to_print += "\t" + bed.rgb
    if bed.exon_count != 0:
        to_print += "\t" + str(len(intervals)) + "\t"
        mapped_starts = []
        mapped_sizes = []
        for interval in intervals:
            mapped_starts.append(str(interval[2] - intervals[0][2]))
            mapped_sizes.append(str(interval[3] - interval[2] + 1))
        to_print += ",".join(mapped_sizes) + "\t" + ",".join(mapped_starts)
    print to_print

# Map interval ['mapping.left', 'mapping.right'] from 'mapping.strand'
# of 'mapping.chromosome' of reference genome to 'mapping.inf_name' genome
# with maximal gap of 'mapping.inf_maxgap' between two informants, if
# 'mapping.inner' then inner mapping (else outer), if 'mapping.always_map'
# force mapping and if 'mapping.all_errors', show all errors
# Uses 'genome_map', 'chr_maps', 'index' and 'ref_name'
def map_interval(mapping, genome_map, chr_maps, index, ref_name, reader):
    if mapping.left > mapping.right:
        add_error("invalid_int")
        sys.exit(0)
    ref_chr_id, ref_chr_len = chr_maps[0][mapping.chromosome]
    inf_id = genome_map[mapping.inf_name]
    # Flags
    reference_gap = 0
    informant_gap = mapping.inf_maxgap
    inner = mapping.inner
    always_map = mapping.always_map
    full = True

    # Get list of references containing given interval, check for errors
    references, inf_errors, continuous =\
        get_references(ref_chr_id, mapping, index, inf_id, reader)
    
    if len(references) == 0:
        add_error("ref_cont")
        sys.exit(0)
    
    # Dictionary that maps chromosome id to its length
    id_map = inverse_to(chr_maps[inf_id])
    # Map left and right ending points of the interval
    left_map = map_position(references, mapping.left, inner, 0, inf_id,
        id_map)
    right_map = map_position(references, mapping.right, not inner,
        len(references)-1, inf_id, id_map)

    # Is this correct mapping of interval of 1 base?
    interval_of_one(mapping.left, mapping.right, left_map,
        right_map, inner, always_map)

    # If there is no inner/outer mapping and -alwaysmap flag is present,
    # try it the other way
    left_map, right_map = map_the_other_way(left_map,
        right_map, mapping, references, inf_id, id_map)
    
    # Error handling after mapping (when 'left_map' and 'right_map' are correct)
    first, last = mapping_errors(left_map, right_map, continuous,
        inf_errors, always_map, chr_maps[inf_id])
    
    global all_errors_count
    if all_errors_count > 0 and mapping.all_errors:
        sys.exit(0)
    # Return all the intervals that input maps to
    return get_mapped_intervals(left_map, right_map, first, last, continuous,
        id_map)

# Here comes the "main function"
if __name__ == "__main__":
    # Usage
    if len(sys.argv) < 2 or (len(sys.argv) < 3 and sys.argv[1] != "help")\
        or (len(sys.argv) >= 3 and len(sys.argv) < 8 \
            and sys.argv[1] != "informants" and sys.argv[1] != "chromosomes") \
        or (len(sys.argv) >= 10 and sys.argv[1] != "map"):
        print "Usage:"
        print "To map interval [L,R] from reference genome to informant:"
        print "\tpython {0} map <header.bin> <from.bgzf>".format(sys.argv[0]),
        print "<chromosome> {+,-} <L> <R> <informant> -maxgap n",
        print "-{inner,outer} [-{allerrors,alwaysmap}]"
        print "To show list of informants:"
        print "\tpython {0} informants <header.bin>".format(sys.argv[0])
        print "To show list of reference chromosomes:"
        print "\tpython {0}".format(sys.argv[0]),
        print "chromosomes <header.bin>"
        sys.exit(0)

    # Display informants aligned to reference
    if sys.argv[1] == "informants":
        genome_map, chr_maps, index, ref_name = get_header_info(sys.argv[2])
        print_informants(genome_map, ref_name)

    # Display chromosomes of reference genome
    if sys.argv[1] == "chromosomes":
        genome_map, chr_maps, index, ref_name = get_header_info(sys.argv[2])
        print_chromosomes(chr_maps, ref_name)
