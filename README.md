<p>This script is designed to compare a set of alignments to an reference alignment, all in FASTA format,
   in order to identify true positive (TP), false negative (FN) and false positive (FP). It works as follows:
   
- 1: after sorting the sequences by names, it extracts all positions from the reference and tested alignments 
- 2: it identifies strictly identical positions between reference and tested alignments as true positives
- 3: for the remaining positions of the tested alignments, it then identifies as true positives positions that only differ by the presence of missing data between reference and tested alignments (by considering all non-ATGC characters as missing data)

<p>Steps 2 and 3 are each repeated with the complement of the positions, but this behaviour can be turned off ('-c option). The output of the script should be treated with caution if the alignments contain high levels of missing data per position.</p>

### requirements
Python 3+ (tested with Python 3.10.9)

### usage
```
python3 compareFASTA -r ref_example.fas -d examples
```
This test run will produce a CSV file called 'output_compareFASTA.csv' that should contain:
```
file,TP,FN,FP
SKA2_D39V__out2.fas,85,2,0
SKA2_D39V__out3.fas,82,5,0
SKA2_D39V__out1.fas,87,0,0
SKA2_D39V__out0.fas,83,4,0
SKA2_D39V__out4.fas,87,0,0
snippy_D39V__out0.fas,82,5,105
snippy_D39V__out1.fas,79,8,114
snippy_D39V__out3.fas,78,9,109
snippy_D39V__out2.fas,83,4,75
snippy_D39V__out4.fas,81,6,87
```

### input files
The script requieres a reference alignment to be compared with ('-r' option) and the name of the directory containing all alignments to test ('-d' option). All alignments should be in FASTA format, with extension '.fas', '.fasta' or '.fna'.

### what the script does not do
- it does not check if the sequence names are identical between the reference alignment and tested alignments
- it does not check if the alignments are true alignments (i.e. sequences of same length) 




