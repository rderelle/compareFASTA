import argparse
import os
import sys
import collections


def parse_args():
    # define and parse command-line arguments
    parser = argparse.ArgumentParser(description='            compareALI v1.0', add_help=False, formatter_class=argparse.RawTextHelpFormatter, epilog=' \n')
    
    common = parser.add_argument_group()
    common.add_argument('-r',            help='reference alignment', metavar='')    
    common.add_argument('-f',            help='alignment to be tested', metavar='')
    common.add_argument('-h','-help',    action="help", help="show this help message and exit")

    args = parser.parse_args()
    
    return args.r, args.f


def pre_checking(dir_, ):
    l = list()
    for file in os.listdir(dir_):
        if file.endswith('.fas') or file.endswith('.fasta') or file.endswith('.fna'):
            l.append(file)
    return l 


def read_fasta(fasta_content):
    name, seq = None, []
    for line in fasta_content:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name.replace('>',''), ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name.replace('>',''), ''.join(seq))


def build_columns(d_seq, len_seq):
    # get sorted list of names
    l_names = list(d_seq)
    l_names.sort()
    
    # build SNP columns
    test_counts = dict()
    for pos in range(len_seq):
        tmp_l = list()
        for name in l_names:
            nucl = d_seq[name][pos]
            if nucl in {'A','C','G','T'}:
                tmp_l.append(nucl)
            else:
                tmp_l.append('-')
        # save it in str format
        str_col = ''.join(tmp_l)
        if str_col in test_counts:
            test_counts[str_col] += 1
        else:
            test_counts[str_col] = 1
    return l_names, test_counts


def complement(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
    return "".join(nn[n] for n in st)

def compare_dna(seq1, seq2):
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            if a != b:
                return False  # Exit immediately on mismatch
    return True  # Return True only if no mismatches were found



if __name__ == "__main__":
    
    ## get all arguments
    ref_ali_file, ali_test = parse_args()
    
    ## check arguments
    if ref_ali_file is None or ali_test is None :
        sys.exit('\n            ERROR: you need to specify an input directory and reference alignment (see -help)\n\n')
    
    ## get reference alignment
    ref_seq = dict()
    with open(ref_ali_file) as fasta_content:
        for name, seq in read_fasta(fasta_content):
            ref_seq[name] = list(seq)
            len_ref_seq = len(seq)
    
    ## get ali to test
    all_seq = dict()
    with open(ali_test) as fasta_content:
        for name, seq in read_fasta(fasta_content):
            all_seq[name] = list(seq)
            len_seq = len(seq)
    
    ## build columns to test
    l_names, test_counts = build_columns(all_seq, len_seq)
    
    ## build reference columns
    l_ref_names, ref_counts = build_columns(ref_seq, len_ref_seq)
    
    ## get identical SNP_columns
    TP = 0
    FP = 0
    FN = 0
    found_identical = 0
    all_SNPs = list(test_counts)
    for SNP in all_SNPs:
        nb = test_counts[SNP]
        if SNP in test_counts and SNP in ref_counts:
            if nb == ref_counts[SNP]:
                TP += nb
                found_identical += nb
                del ref_counts[SNP]
                del test_counts[SNP]
            elif nb > ref_counts[SNP]:
                TP += ref_counts[SNP]
                found_identical += ref_counts[SNP]
                test_counts[SNP] -= ref_counts[SNP]
                del ref_counts[SNP]
            elif nb < ref_counts[SNP]:
                TP += nb
                found_identical += nb
                ref_counts[SNP] -= nb
                del test_counts[SNP]
    
    ## get identical SNP_columns in complement
    all_SNPs = list(test_counts)
    for SNP in all_SNPs:
        compl = complement(SNP)
        nb = test_counts[SNP]
        if compl in ref_counts:
            if nb == ref_counts[compl]:
                TP += nb
                found_identical += nb
                del ref_counts[compl]
                del test_counts[SNP]
            elif nb > ref_counts[compl]:
                TP += ref_counts[compl]
                found_identical += ref_counts[compl]
                test_counts[SNP] -= ref_counts[compl]
                del ref_counts[compl]
            else:
                TP += nb
                found_identical += nb
                ref_counts[compl] -= nb
                del test_counts[SNP]
   
    ## check conflict between SNPs and combine them
    d_no_conflict = dict()
    d_ref_included = dict()
    for SNP in test_counts:
        compl = complement(SNP)
        tmp_l = list()
        for ref_SNP in ref_counts:
            no_conflict = compare_dna(SNP, ref_SNP)
            no_conflict2 = compare_dna(compl, ref_SNP)
            if no_conflict or no_conflict2:
                tmp_l.append(ref_SNP)
                if ref_SNP in  d_ref_included:
                    d_ref_included[ref_SNP] += 1
                else:
                    d_ref_included[ref_SNP] = 1
        
        if tmp_l:
            d_no_conflict[SNP] = tmp_l
    
    
    l_test = list(test_counts)
    for SNP in l_test:
        if not SNP in d_no_conflict:
            FP += test_counts[SNP]
            del test_counts[SNP]
    
    ref_seen = 0
    for SNP, l in d_no_conflict.items():
        for ref_SNP in l:
            ref_seen += ref_counts[ref_SNP] / d_ref_included[ref_SNP]
    
    remaning_test = sum(test_counts.values())
    remaning_ref  = sum(ref_counts.values())
    
    TP += min(ref_seen,remaning_test)
    new_FP = 0
    if ref_seen > remaning_ref:
        new_FP = ref_seen - remaning_ref
        FP += new_FP
    
    FN += (remaning_ref + new_FP - remaning_test)
        
    sensitivity = TP / (TP + FN)
    precision = TP / (TP + FP)
    
    # output values
    print('reference_alignment	tested_alignment	TP	FN	FP	sensitivity	precision')
    print(ref_ali_file + '	' + ali_test + '	' + str(TP)+ '	' + str(FN) + '	' + str(FP) + '	' + str(round(sensitivity,4)) + '	' + str(round(precision,4)))
