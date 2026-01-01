"""
gRNA Design Project Draft 7

    This script scans the selected record from the TP53 FASTA file, and 
    identifies candidate guide RNA's for the CRISPR SpCas9 complex.
    
    Important Biological Assumptions: 
        - Cas Variant: SpCas9
        - PAM Sequence: NGG 
        - gRNA Length: 20 nucleotides
        - Normalization to plus strand: CRISPR complex is direction specific
        
    Additions in Draft 7: 
        - Added summary statistics for total and passing guide counts
        - Improved output readability and organization
        - Updated histogram visualization, by including all candidate
            guide RNAs
        - Added biologically relevant GC% thresholds markers at 40 and 60%. 
            
    Current Limitations/Future Additions
    -Currently only works for selected TP53 FASTA file 

"""

from Bio import SeqIO
import matplotlib.pyplot as plt 

all_guides = []
passing_guides = []



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
            
        #all_guides contains all the possible candidates for gRNAs
        #regardless of GC%
            
            all_guides.append({
                "pam_position":pam_index,
                "guide_seq": str(guide_seq),
                "gc_percent":gc_content
                })
            
        #the ideal range for binding efficiency is 40% to 60% 
        #thus we only select those that fall in this range for passing_guides
            
            if (gc_content <= 60) and (gc_content >= 40):
                #print("PAM Position: ", pam_index)
                #print ("Associated Guide RNA: ", guide_seq)
                #print("Total GC Percentage: ",gc_content,"%")
                #print("")
                
                passing_guides.append({
                    "pam_position":pam_index,
                    "guide_seq": str(guide_seq),
                    "gc_percent":gc_content
                    })
                
with open("gRNA_doc.tsv","w") as file_name:
    file_name.write("pam_position\tguide_sequence\tgc_percentage\n")
                   
    for guide in passing_guides:
        file_name.write(f"{guide['pam_position']}\t{guide['guide_seq']}\t{guide['gc_percent']}\n")
                       
all_gc = [guide["gc_percent"] for guide in all_guides]


total_candidates = len(all_guides)
passing_candidates = len(passing_guides)
 


if (total_candidates == 0):
    print("No Candidate Guides were found, thus no plot was generated.")

else:
    
    #printing out the percentasge of viable candidates for the user
   
    print("Total candidates found: ", total_candidates)
    print("Viable candidates post filtering(40-60%): ", passing_candidates)
    passing_percentage = passing_candidates/total_candidates*100
    print("The percentage of appropriate candidates is: ", passing_percentage, "%")
    
   
    #creating a histogram for all the guides found, with a filter at the 40-60% 
    
    plt.figure()
    plt.hist(all_gc, bins = 20, edgecolor='black')
    
    plt.axvline(40, color='black', linestyle='--', label='Lower GC cutoff (40%)')
    plt.axvline(60, color='black', linestyle='--', label='Upper GC cutoff (60%)')

    plt.xlabel("GC Percentage")
    plt.ylabel("Number of Guides")
    plt.title("GC% Distribution of Candidate TP53 Guide RNAs (Filter Window: 40-60%)")
    plt.legend()
    plt.savefig("gc_distribution.png", dpi = 300)
    plt.show()                 


        
