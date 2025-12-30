"""
gRNA Design Project Draft 6

    This script scans the selected record from the TP53 FASTA file, and 
    identifies candidate guide RNA's for the CRISPR SpCas9 complex.
    
    Important Biological Assumptions: 
        - Cas Variant: SpCas9
        - PAM Sequence: NGG 
        - gRNA Length: 20 nucleotides
        - Normalization to plus strand: CRISPR complex is direction specific
        
    Additions in Draft 6: 
        - Added histogram of GC% values for filtered gRNAs
            
    Current Limitations/Future Additions
    -Currently only works for selected TP53 FASTA file 

"""

from Bio import SeqIO
import matplotlib.pyplot as plt 

guides_list = []


#using the SeqIO parse function, we're able to iterate through the two records
#in the file and put them into a list format underneath. 
tp53 = list(SeqIO.parse("TP53.fna","fasta"))


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
        
        pam_index = counter
        
        if counter - 20 >= 0:       #gRNA length is 20 nt. It will be found 
                                    #upstream of the PAM.
            
            guide_seq = working_seq[counter-20:counter]
            
            total_gc = 0
            start = 0
            
            #adding a GC% filter to measure binding stability of gRNA
            
            for i in range(start,len(guide_seq)):
                if (guide_seq[i] == "G") or (guide_seq[i] == "C"):
                    total_gc = total_gc+1
            gc_content = (total_gc/len(guide_seq))*100
            
            #the ideal range for binding efficiency is 40% to 60% 
            #thus we only select those that fall in this range 
            
            if (gc_content <= 60) and (gc_content >= 40):
                print("PAM Position: ", pam_index)
                print ("Associated Guide RNA: ", guide_seq)
                print("Total GC Percentage: ",gc_content,"%")
                print("")
                
                guides_list.append({
                    "pam_position":pam_index,
                    "guide_seq": str(guide_seq),
                    "gc_percent":gc_content
                    })
                
with open("gRNA_doc.tsv","w") as file_name:
    file_name.write("pam_position\tguide_sequence\tgc_percentage\n")
                   
    for guide in guides_list:
        file_name.write(f"{guide['pam_position']}\t{guide['guide_seq']}\t{guide['gc_percent']}\n")
                       
gc_values = [guide["gc_percent"] for guide in guides_list]

 
if (len(gc_values) == 0):
    print("No Guides passed through the filter, thus no plot was generated.")
else :
    plt.figure()
    plt.hist(gc_values, bins = 20)
    plt.axvline(40)
    plt.axvline(60)
    plt.xlabel("GC Percentage")
    plt.ylabel("Number of Guides")
    plt.title("GC% Distribution of Candidate TP53 Guide RNAs (within 40-60%)")
    plt.savefig("gc_distribution.png", dpi = 300)
    plt.show()             
        


        
