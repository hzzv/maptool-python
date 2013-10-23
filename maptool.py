# -*- coding: utf-8 -*-
# Copyright 2013 by Petra Kubincova

import read_maf as preprocessing
import seq_structures as seq
import bgzf_tool as bgzftool
import mapping
import sys
import os


def parse_mapping_input(args):
    try:
        import optparse
        parser = optparse.OptionParser()
        parser.add_option("--maxgap", "-m", dest="maxgap", metavar="N",
            default="0")
        parser.add_option("--outer", "-o", dest="inner",
            action="store_const", const="-outer", default="-inner")
        parser.add_option("--alwaysmap", dest="alwaysmap",
            action="store_true", default=False)
        parser.add_option("--allerrors", dest="allerrors",
            action="store_true", default=False)
        parsed, redundant = parser.parse_args(args)
    except ImportError as e:
        try:
            import argparse
            parser = argparse.ArgumentParser()
            parser.add_argument("--maxgap", "-m", dest="maxgap", metavar="N",
                default="0")
            parser.add_argument("--outer", "-o", dest="inner",
                action="store_const", const="-outer", default="-inner")
            parser.add_argument("--alwaysmap", dest="alwaysmap",
                action="store_true", default=False)
            parser.add_argument("--allerrors", dest="allerrors",
                action="store_true", default=False)
            parsed = parser.parse_args(args)
        except ImportError as e:
            print >> sys.stderr, "Error: neither of optparse and argparse " +\
                "modules is available."
            sys.exit(0)
    last_option = [""]
    if parsed.alwaysmap:
        last_option = ["-alwaysmap"]
    elif parsed.allerrors:
        last_option = ["-allerrors"]
    return [parsed.maxgap, parsed.inner] + last_option

def check_files(file1, file2=False):
    nonexistent = []
    if not os.path.isfile(file1):
        nonexistent.append(file1)
    if file2 and not os.path.isfile(file2):
        nonexistent.append(file2)
    if len(nonexistent) > 0:
        print >> sys.stderr, "Error: File(s) " + ", ".join(nonexistent) +\
            " is (are) nonexistent"
        sys.exit(0)

def check_inf_chr(chromosome, informant, genome_map, chr_maps):
    error = False
    if chromosome and not chromosome in chr_maps[0]:
        print >> sys.stderr, "Error:", chromosome, "is not a known chromosome"
        error = True
    if informant and not informant in genome_map:
        print >> sys.stderr, "Error:", informant, "is not a known informant"
        error = True
    if error:
        sys.exit(0)

# Displays help about 'command'
def help(command, program_name):
    if command == "preprocess":
        print "Usage:",
        print "python {0} preprocess <from.maf> <header.bin> <to.bgzf>".\
            format(program_name)
        print "Creates binary header and compressed file with alignments",
        print "from source .maf. Needed to be run only once for each .maf."
        print "<from.maf> is source .maf file to be preprocessed"
        print "<to.bgzf> is path where to write compressed output"
        print "<header.bin> is path where to write header information",
        print "(about data in compressed file)"
    elif command == "informants":
        print "Usage:",
        print "python {0} informants <header.bin>".format(program_name)
        print "Prints list of informants contained in .maf alignment to",
        print "which was by 'preprocess' command created <header.bin>."
        print "<header.bin> is path to header file"
    elif command == "chromosomes":
        print "Usage:",
        print "python {0} chromosomes <header.bin>".format(program_name)
        print "Prints list of chromosomes of reference genome contained",
        print "in .maf alignment to which was by 'preprocess' command",
        print "created <header.bin>."
        print "<header.bin> is path to header file"
    elif command == "bed":
        print "Usage:",
        print "python {0} bed <header.bin>".format(program_name),
        print "<from.bgzf> <informant> [--maxgap N] [--outer] [--allerrors]",
        print "[--alwaysmap]"
        print "Reads from stdin, writes to stdout and stderr.",
        print "Expects and produces",
        print "BED format. For each line of input prints one line",
        print "of output (if mapping is possible),",
        print "which may contain more intervals (if interval in",
        print "reference genome maps to more intervals in informant genome,",
        print "or if on input were more intervals (exons)).",
        print "Also for each line of input prints at least one line to",
        print "stderr: information whether was mapping successful or which",
        print "errors it did encounter."
        print "If N (maxgap) is -1, gaps of arbitrary length are allowed."
        print "If --outer flag is present, left ending points",
        print "of the input intervals will be mapped to the left and right",
        print "ending points to the right. Otherwise left will be mapped",
        print "to the right and right to the left."
        print "If --allerrors",
        print "flag is present, to stderr will be printed as many errors as",
        print "were encountered during the whole mapping process. If",
        print "--alwaysmap flag is present, all encountered errors will be",
        print "printed on stderr, and if none is terminal, mapping will be",
        print "carried out. If none of these two flags is present, after",
        print "encountering some errors mapping will be stopped (not all",
        print "errors will be printed on stderr)."
        print "Printout on stderr will be in table format with columns",
        print "'Name, Whether or not was that line mapped, [Error_id],",
        print "[Error description], [Optional error description]', separated",
        print "by tabs."
    elif command == "help":
        print "Usage:",
        print "python {0} help [<command>]".format(program_name)
        print "Prints help about <command>."
    else:
        print "'{0}' is not a valid command.".format(command)

if __name__ == "__main__":
    # Usage
    if len(sys.argv) <= 1 or (len(sys.argv) == 2 and sys.argv[1] != "help" and\
        sys.argv[1] != "-h" and sys.argv[1] != "--help")\
        or (len(sys.argv) == 3 and sys.argv[1] != "informants" and\
            sys.argv[1] != "chromosomes" and sys.argv[1] != "help") \
        or len(sys.argv) == 4 \
        or (len(sys.argv) == 5 and sys.argv[1] != "preprocess" and\
            sys.argv[1] != "bed") \
        or (10 >= len(sys.argv) >= 6 and sys.argv[1] != "bed") \
        or len(sys.argv) > 10:
        print "Usage:"
        print "To preprocess .maf file:"
        print "\tpython {0} preprocess <from.maf> <header.bin> <to.bgzf>".\
            format(sys.argv[0])
        print "To map set of intervals from input .bed file to set",
        print "of intervals to output .bed file:"
        print "\tpython {0} bed <header.bin> <from.bgzf>".format(sys.argv[0]),
        print "<informant> [--maxgap N] [--outer] [--allerrors] [--alwaysmap]"
        print "To show list of informants:"
        print "\tpython {0} informants <header.bin>".format(sys.argv[0])
        print "To show list of reference chromosomes:"
        print "\tpython {0} chromosomes <header.bin>".format(sys.argv[0])
        print "To show help:"
        print "\tpython {0} help [<command>]".format(sys.argv[0])
        sys.exit(0)
    
    # Help
    if sys.argv[1] == "help" or sys.argv[1] == "-h" or sys.argv[1] == "--help" \
        and 2 <= len(sys.argv) <= 3:
        if len(sys.argv) == 2:
            print "Maptool, version 1. Created by Petra Kubincova."
            print "Needed: Python2.x, x >= 6, and Biopython module bgzf.py"
        else:
            help(sys.argv[2], sys.argv[0])

    # Preprocessing
    elif sys.argv[1] == "preprocess" and len(sys.argv) == 5:
        check_files(sys.argv[2])
        preprocessing.do_all_the_magic(["preprocess"] + [sys.argv[2]] +\
            [sys.argv[4]] + [sys.argv[3]])

    # Map or display chromosomes or informants
    else:
        sys.argv.append(False)
        check_files(sys.argv[2], sys.argv[3])
        sys.argv.pop(-1)
        # Read header, get name of reference genome. Will be needed in any case
        try:
            genome_map, chr_maps, index, ref_name =\
                mapping.get_header_info(sys.argv[2])
        except IOError as e:
            print >> sys.stderr, "Wrong filename:", e.filename
            sys.exit(0)

        # Display informants aligned to reference
        if sys.argv[1] == "informants":
            mapping.print_informants(genome_map, ref_name)

        # Dislplay chromosomes of reference genome
        elif sys.argv[1] == "chromosomes":
            mapping.print_chromosomes(chr_maps, ref_name)
        
        # Map intervals in BED format
        elif sys.argv[1] == "bed":
            check_inf_chr(False, sys.argv[4], genome_map, chr_maps)
            try:
                reader = bgzftool.ReadBgzf(sys.argv[3])
            except IOError as e:
                print >> sys.stderr, "Wrong filename:", e.filename
                sys.exit(0)
            options = parse_mapping_input(sys.argv[5:])
            lines = sys.stdin.readlines()
            line_counter = 0
            for line in lines:
                name = str(line_counter + 1)
                # Parse input
                user_input = [x for x in line.strip().split() if x != ""]
                if len(user_input) < 3 or user_input[0] == "browser" or\
                    user_input[0] == "track" or user_input[0][0] == "#":
                    continue
                if not user_input[0] in chr_maps[0]:
                    mapping.add_error("invalid_chr")
                    line_counter += 1
                    mapping.print_errors(str(line_counter),
                        options[2] == "-alwaysmap")
                    continue
                try:
                    line_counter += 1
                    bed = seq.BedQuery(user_input, line_counter)
                    name = bed.name
                    # If line contains only one interval to map
                    if not bed.exon_count:
                        exon_maps = [bed.get_map()]
                    # If line contains more intervals (exons)
                    else:
                        exon_maps = bed.get_exon_maps()
                    intervals = []
                    for to_map in exon_maps:
                        new_intervals = mapping.map_interval(seq.Mapping(
                            to_map, sys.argv[4], options,\
                            (1, 2)[bed.exon_count > 0]), genome_map,\
                            chr_maps, index, ref_name, reader)
                        if new_intervals == -1 or len(new_intervals) == 0:
                            mapping.add_error("incorrect_int")
                            sys.exit(0)
                        else:
                            intervals.extend(new_intervals)
                    if not mapping.check_exons(intervals):
                        mapping.add_error("incorrect_int")
                        sys.exit(0)
                    # Map thick interval if is one on the input
                    thick_interval = [(-1, -1, -1, -1)]
                    if int(bed.rthick) >= 0:
                        lthick_map = seq.Mapping(bed.get_left_thick(),\
                            sys.argv[4], ["-1", "-outer"] + options[2:], 1)
                        rthick_map = seq.Mapping(bed.get_right_thick(),\
                            sys.argv[4], ["-1", "-outer"] + options[2:], 1)
                        lthick = mapping.map_interval(lthick_map,
                            genome_map, chr_maps, index, ref_name, reader)
                        rthick = mapping.map_interval(rthick_map,
                            genome_map, chr_maps, index, ref_name, reader)
                        # Is thick interval mapped correctly?
                        if lthick == -1 or rthick == -1 or\
                            lthick[0][0:2] != rthick[0][0:2] or\
                            lthick[0][2] >= rthick[0][3]:
                            mapping.add_error("thick")
                            sys.exit(0)
                        thick_interval = [lthick[0][0:3] + rthick[0][3:4]]
                    mapping.print_bed(intervals, thick_interval,
                        bed, chr_maps[genome_map[sys.argv[4]]])
                    mapping.print_errors(name, options[2] == "-alwaysmap")
                except SystemExit as e:
                    mapping.print_errors(name, options[2] == "-alwaysmap")
                except Exception as e:
                    print >> sys.stderr, name + "\tno\tanother\t" + str(e)
