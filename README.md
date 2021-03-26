# Bioinformatics

### Project 1
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

### Project 2
#### Steps followed
1. We have fragments (reads) of a dna sequence in a fasta file. Read the fasta file and collect all the reads in to a list.
2. First derive the individual reads that have overlaps between them.
3. Gather such reads together to create a consensus out of them. In this step, we create the shortest path for all the reads that have overlaps between them.
4. Once we have the consensus for the reads from step2, we do this recursively so that we get a final contig which is the shortest path encompassing all the reads without any repetition.
5. We use Pydna for this task. Below features of Pydna were used
<br>a. To derive overlaps between reads (multiple)
<br>b. To assemble reads and find the shortest path.
6. We use below cases to test the reconstructed sequence
<br>a. Reconstructed sequence should cover every read from the input file. We validate every read, if it is a substring of the reconstructed sequence.
<br>b. Reconstructed sequence should not have repetitions. We find how many times a read exists in the reconstructed sequence. Each read should exist only once.
7. Write the reconstructed sequence in to an output file in the below format.

">reconstructed genome sequence <br>
"ouput assembled sequence
