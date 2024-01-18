import argparse
import os
import sys


def parse_args():
    # define and parse command-line arguments
    parser = argparse.ArgumentParser(description='            compareALI', add_help=False, formatter_class=argparse.RawTextHelpFormatter, epilog=' \n')
    
    common = parser.add_argument_group()
    common.add_argument('-r',            help='reference alignment', metavar='')    
    common.add_argument('-d',            help='directory containing alignments to be tested', metavar='')
    common.add_argument('-c',            help='turn off complement testing', action="store_false")
    common.add_argument('-h','-help',    action="help", help="show this help message and exit")

    args = parser.parse_args()
    
    return args.r, args.d, args.c


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
    SNP_columns = list()
    for pos in range(len_seq):
        tmp_l = list()
        for name in l_names:
            nucl = d_seq[name][pos]
            if nucl in {'A','C','G','T'}:
                tmp_l.append(nucl)
            else:
                tmp_l.append('-')
        # save it in str format
        SNP_columns.append(''.join(tmp_l))
    return l_names, SNP_columns


def complement(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
    return "".join(nn[n] for n in st)




if __name__ == "__main__":
    
    ## get all arguments
    ref_ali_file, dir_test, complement_testing = parse_args()
    
    ## check arguments
    if ref_ali_file is None or dir_test is None :
        sys.exit('\n            ERROR: you need to specify an input directory and reference alignment (see -help)\n\n')
    
    ## get reference alignment
    ref_seq = dict()
    with open(ref_ali_file) as fasta_content:
        for name, seq in read_fasta(fasta_content):
            ref_seq[name] = list(seq)
            len_ref_seq = len(seq)
    
    ## create output file
    out = open('output_compareALI.csv', 'w+')
    out.write('file,TP,FN,FP\n')
    
    ## get list of files to test
    l_ali_files = pre_checking(dir_test)
    
    ## compare each file 1 by 1 against the reference file
    for file in l_ali_files:

        ## get ali to test
        all_seq = dict()
        with open(dir_test + '/' + file) as fasta_content:
            for name, seq in read_fasta(fasta_content):
                all_seq[name] = list(seq)
                len_seq = len(seq)
        
        ## build columns to test
        l_names, SNP_columns = build_columns(all_seq, len_seq)
        
        ## build reference columns
        l_ref_names, ref2_SNP_columns = build_columns(ref_seq, len_ref_seq)
        
        ## get identical SNP_columns
        TP = 0
        other_SNP_columns = list()
        for col in SNP_columns:
            found = False
            compl_col = complement(col)
            for i,col2 in enumerate(ref2_SNP_columns):
                if col == col2:
                    found = True
                    index_ = i
                    break
                # try complement
                elif complement_testing and compl_col == col2:
                    found = True
                    index_ = i
                    break
            
            if found:
                TP += 1
                del ref2_SNP_columns[index_]
            else:
                other_SNP_columns.append(col)
        
        
        # get SNP columns that only differ by missing characters
        left_over = list()
        for col in other_SNP_columns:
            found = False
            compl_col = complement(col)
            for i,col2 in enumerate(ref2_SNP_columns):
                same = True
                for n in range(len(col)):
                    if col[n] != '-' and col2[n] != '-' and col[n] != col2[n]:
                         same = False
                         break
                
                if same:
                    found = True
                    index_ = i
                    break
                
                # try complement
                elif complement_testing:
                    same = True
                    for n in range(len(compl_col)):
                        if compl_col[n] != '-' and col2[n] != '-' and compl_col[n] != col2[n]:
                             same = False
                             break
                    if same:
                        found = True
                        index_ = i
                        break
                
            if found:
                TP += 1
                del ref2_SNP_columns[index_]
            else:
                left_over.append(col)
        
        FN = len(ref2_SNP_columns)    
        FP = len(left_over)    
        
        out.write(file + ',' + str(TP)+ ',' + str(FN) + ',' + str(FP) + '\n')

out.close()
    
















