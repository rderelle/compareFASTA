<p>This script is designed to compare a set of alignments to a reference alignment, all in FASTA format,
   in order to identify true positive (TP), false negative (FN) and false positive (FP). It works as follows:
   
- 1: after sorting the sequences by names, it extracts all positions from the reference and tested alignments 
- 2: it identifies strictly identical positions between reference and tested alignments as true positives
- 3: for the remaining positions of the tested alignments, it then identifies as true positives positions that only differ by the presence of missing data between reference and tested alignments (by considering all non-ATGC characters as missing data)

<p>Steps 2 and 3 are each repeated with the complement of the positions. The output of the script should be treated with caution if the alignments contain high levels of missing data per position.</p>

### requirements
Python 3+ (tested with Python 3.10.9)

### Example usage
```
python3 compareALI.py -r ref_example.fas -f SKA2_D39V__out0.fas
```
This command line will output:
```
reference_alignment	tested_alignment	TP	FN	FP	sensitivity	precision
ref_example.fas	SKA2_D39V__out0.fas	83	4	0	0.954	1.0
```

### input files
The script requieres a reference alignment to be compared with ('-r' option) and an alignment to be tested ('-f option). All alignments should be in FASTA format.

### what the script does not do
- it does not check if the sequence names are identical between the reference alignment and tested alignments
- it does not check if the alignments are true alignments (i.e. sequences of same length) 




