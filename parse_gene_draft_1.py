"""
This script is meant to take fasta files for DNA sequences and screen for potential Guide RNA sequences" 
gRNAs always have an NGG PAM (Protospacer Adjacent Motif) which is upstream of the 20 nucleotide sequence. 
- First section imports BioPython and the SeqIO modulem, used to read and Parse sequences, 
Different parts of SeqRecord: 
 record.id → the ID from the FASTA header

record.description → longer description from the header

record.seq → the actual DNA sequence (string-like)

record.annotations → extra info (for formats like GenBank)   
"""

import os
from Bio import SeqIO # Importing the sub module from BioPython 

#os.chdir ("/Users/jaichowdhary/gRNA_design") # uses the os module to to always change the directory in this script to "gRNA_design".
"""
TP53_file = SeqIO.read("TP53.fna","fasta")      # using the SeqIO.read() function to read the FASTA file
print (TP53_file.id)
print (TP53_file.seq[:50])     # printing the actual sequence but using the python slice notation to print only the first 50 characters. 
This code was meant for .fna file with only one record in it. 
"""

#for sequences in SeqIO.parse("TP53.fna","fasta"):
    #print (sequences.id)

sequence_list = list(SeqIO.parse("TP53.fna","fasta")) # turns the output of the func SeqIO.parse() into a list
print (sequence_list[0].seq[:50])  # prints the 0th index position of the list i.e. the first sequence (.seq), and the only the 1st 50 characters. 

print (len(sequence_list[0].seq)) #prints the length of the sequence 
print (sequence_list[0].seq[630:700]) # prints the sequences between 630 and 700

print(sequence_list[0].seq[-50:]) #prints the last 50 nucleotides/ Note: for slicing a string = [beginning:ending]

sequence_slice = sequence_list[0].seq[50:80] # assigns this variable to the 30 nucleotides that have been sliced out of the sequence
print (sequence_slice)

