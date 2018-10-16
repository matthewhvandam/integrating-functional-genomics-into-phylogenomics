#!/usr/bin/env python
# coding: utf-8
# read from stdin and write to stdout inserting introns between exon lines of a gff3 format file
# written by JB Henderson

import sys

def add_introns(infile, outfile):
    def to_int(digits_str):  # do not throw a ValueError exception if input happens to not be a string of just digits
        return int(digits_str) if digits_str.isdigit() else 0

    cur_flds = []
    for gff_ln in infile:
        if len(gff_ln) < 3 or gff_ln[0] == "#":  # output empty or comment lines
            outfile.write(gff_ln)
            continue
        
        lst_flds = list(cur_flds)  # copy the last non-comment line fields
        cur_flds = gff_ln.split("\t")
        
        too_short = len(cur_flds) < 5 or len(lst_flds) < 5
        if too_short or cur_flds[2].lower() != "exon" or lst_flds[2].lower() != "exon":
            outfile.write(gff_ln)
            continue
        
        # last line and current line both are exon lines, this is where to insert the intron line
        assert cur_flds[2].lower() == "exon" and lst_flds[2].lower() == "exon"
        
        cur_beg = to_int(cur_flds[3]); cur_end = to_int(cur_flds[4])
        lst_beg = to_int(lst_flds[3]); lst_end = to_int(lst_flds[4])
        
        if cur_beg==0 or cur_end==0 or lst_beg==0 or lst_end==0:
            # we expected positive integers but didn't get them
            outfile.write(gff_ln)
            continue
        
        if lst_end < cur_beg:  # positive strand
            int_beg = lst_end + 1
            int_end = cur_beg - 1
        else:  # negative strand
            int_beg = cur_end + 1
            int_end = lst_beg - 1
        
        intron = "\t".join([ lst_flds[0], "Insert","intron", str(int_beg),str(int_end), lst_flds[5],lst_flds[6],lst_flds[7] ]) + "\n"
        outfile.write(intron)  # output intron line
        outfile.write(gff_ln)  # output current exon line

if __name__ == '__main__':
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) # don't complain if pipe closes output (head or less commands will cause this)
    try:
        add_introns(sys.stdin, sys.stdout)
    except:
        pass
