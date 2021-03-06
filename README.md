# Bioinformatics

#### Steps followed
1. Human genome file is huge around 3gb. To extract the sequence by reading this file would be time consuming.
2. Extracted the sequence for a specific chromosome in to temporary fasta file.
3. For the given chromosome, below steps are followed to extract the protein sequence.
4. Extracted sequences for every exon frame, excep for frame = -1. 
5. For Negative strand - 
<br>a. Read the frames, start and end positions in the reverse order.
<br>b. Using exon start and end positions, extract the sequence.
<br>c. Gather sequences for each exon frame in to a list and join in to a single sequence.
<br>d. Compute the starting position of neucleotide sequence using codon end position (cdsEnd from annotation table).
<br>e. Substring the sequence from this new start position.
6. For Positive strand -
<br>a. Read the frames, start and end positions in the normal order.
<br>b. Using exon start and end positions, extract the sequence.
<br>c. Gather sequences for each exon frame in to a list and join in to a single sequence.
<br>d. Compute the starting position of neucleotide sequence using codon start position (cdsSt from annotation table).
<br>e. Substring the sequence from this new start position.    
7. Translate the above extracted squence in to protein sequence using codon table. Stop the translation once a stop codon is read. 
8. Write the protein sequence in to an output file in the below format.

Output Format:
>name:name2 <br>
Protein sequence