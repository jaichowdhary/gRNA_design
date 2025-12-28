
"""
gRNA Design Project Draft 2 

Updated Draft: 
    This draft works on the previous draft by still parsing through the FASTA 
    file. However it now preserves CRISPR logic by normalizing the strand 
    first. 
    
    Important additions include looping through record to identify PAM sites 
    and slicing out the appropriate 20 basepair sequences that can be used as 
    Guide RNA sequences. 
"""

from Bio import SeqIO

#using the SeqIO parse function, we're able to iterate through the two records
#in the file and put them into a list format underneath. 
tp53 = list(SeqIO.parse("TP53.fna","fasta"))

print(tp53[0].id)
print(tp53[0].seq[0:50])

print(tp53[1].id)
print(tp53[1].seq[0:50])


#Record of choice is NC_000017.11
record = tp53[0]

#the working record is the reverse complement of the minus strand to make sure 
#that CRISPR strand reading logic is preserved. 

working_record = record.reverse_complement()
print(working_record)

#using the reverse complement function creates a whole new record(object)
#which needs to be processed similarly to the original record.

working_seq = working_record.seq 

#iterate through the sequence and find "NGG" PAM Sites
#spCas9 PAM = NGG, where "N" can be any nucleotide. 

pam_size = 3

for counter in range(0,len(working_seq)-pam_size) :
   
    if (working_seq[counter+1] == "G") and (working_seq[counter+2]=="G"):
        
        PamIndexPos = counter
        
        if counter - 20 >= 0:       #gRNA length is 20 nt. 
            
            guide_seq = working_seq[counter-20:counter]
            print("PAM Position: ", PamIndexPos)
            print ("Associated Guide RNA: ", guide_seq)
        
        

    
        
