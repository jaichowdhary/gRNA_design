"""
gRNA Design Project Draft 8

    This script scans the selected record from the TP53 FASTA file, and 
    identifies candidate guide RNA's for the CRISPR SpCas9 complex.
    
    Important Biological Assumptions: 
        - Cas Variant: SpCas9
        - PAM Sequence: NGG 
        - gRNA Length: 20 nucleotides
        - Normalization to plus strand: CRISPR complex is direction specific
        
    Additions in Draft 8: 
        - Refactored script into modular, reusable functions
        - Separated scanning, filtering, output, and visualization logic
        - Improved code readability and maintainability
            
    Current Limitations/Future Additions
    - Currently FASTA path is hardcoded

"""

from Bio import SeqIO
import matplotlib.pyplot as plt 



#CONSTANTS
PAM_SIZE = 3 
GUIDE_LEN = 20
MIN_GC = 40 
MAX_GC = 60 



def main():
    
    file_records = load_records("TP53.fna")
    
    selected_record = select_record(file_records)
    
    normalized_sequence = normalize_plus_strand(selected_record)
    
    all_guides = find_candidate_guides(normalized_sequence)

    viable_guides = filter_guides_by_gc(all_guides, MIN_GC, MAX_GC)
    
    guides_to_file(viable_guides)
    
    plot_histogram_and_summary(all_guides, viable_guides)

#using the SeqIO parse function, we're able to iterate through the two records
#in the file and put them into a list format underneath. 

def load_records(fasta_pathway):
    records = list(SeqIO.parse(fasta_pathway,"fasta"))
    return records


def select_record(records, record_index=0):
    
    #Record of choice is NC_000017.11
    record_of_choice = records[record_index]
    return record_of_choice

#the working record is the reverse complement of the minus strand to make sure 
#that CRISPR strand reading logic is preserved. 

def normalize_plus_strand(selected_record):
    
    #using the reverse complement function creates a whole new record(object)
    #which needs to be processed similarly to the original record.
    
    working_seq = selected_record.seq.reverse_complement()
    return working_seq



#iterate through the sequence and find "NGG" PAM Sites
#spCas9 PAM = NGG, where "N" can be any nucleotide. 


def find_candidate_guides(normalized_sequence):
    
    all_guides = []
    
    for counter in range(0,len(normalized_sequence)-PAM_SIZE+1) :
       
        if (normalized_sequence[counter+1] == "G") and (normalized_sequence[counter+2]=="G"):
            
            pam_index = counter
            
            if counter - GUIDE_LEN >= 0:         #gRNA length is 20 nt. It will be found 
                                                #upstream of the PAM.
                
                guide_seq = normalized_sequence[counter-GUIDE_LEN:counter]
                
                total_gc = 0
                
                
                #calculating the GC% percentage to measure binding stability of each gRNA    
                for i in range(len(guide_seq)):
                    if (guide_seq[i] == "G") or (guide_seq[i] == "C"):
                        total_gc = total_gc+1
                gc_content = (total_gc/len(guide_seq))*100
                
                all_guides.append({
                    "pam_position":pam_index,
                    "guide_seq": str(guide_seq),
                    "gc_percent":gc_content
                    })
    return all_guides
                
                
def filter_guides_by_gc(all_guides, min_gc, max_gc):            
                
#the ideal range for binding efficiency is 40% to 60%
#thus we only select those that fall in this range for passing_guides
    
    passing_guides = []
    
    for guide in all_guides:
                      
        if (guide["gc_percent"] <= max_gc) and (guide["gc_percent"] >= min_gc):
                    
            passing_guides.append(guide)
                
    return passing_guides
    

def guides_to_file(viable_guides):                  

    with open("gRNA_doc.tsv","w") as file_name:
    
        file_name.write("pam_position\tguide_sequence\tgc_percentage\n")
                   
        for guide in viable_guides:
        
            file_name.write(f"{guide['pam_position']}\t"
                            f"{guide['guide_seq']}\t"
                            f"{guide['gc_percent']:.2f}\n")
           


def plot_histogram_and_summary(all_guides, viable_guides):
            
    
    total_candidates = len(all_guides)
    passing_candidates = len(viable_guides)
 


    if (total_candidates == 0):
        print("No Candidate Guides were found, thus no plot was generated.")
        return 
    
    else:
    
        #extracting GC content for all selected guides
        all_gc = [guide["gc_percent"] for guide in all_guides]
       
        #printing out the percentasge of viable candidates for the user
     
        print("Total candidates found: ", total_candidates)
        print("Viable candidates post filtering(40-60%): ", passing_candidates)
        
        passing_percentage = passing_candidates/total_candidates*100
        print(f"The percentage of appropriate candidates is: {passing_percentage:.2f}%")
        
       
        #plot histogram for all the guides found, with a filter at the 40-60% 
      
        plt.figure()
        plt.hist(all_gc, bins = 20, edgecolor='black')
        
        plt.axvline(MIN_GC, color='black', linestyle='--', label='Lower GC cutoff (40%)')
        plt.axvline(MAX_GC, color='black', linestyle='--', label='Upper GC cutoff (60%)')
    
        plt.xlabel("GC Percentage")
        plt.ylabel("Number of Guides")
        plt.title("GC% Distribution of Candidate TP53 Guide RNAs (Filter Window: 40-60%)")
        plt.legend()
        plt.savefig("gc_distribution.png", dpi = 300)
        plt.show()                 


        
if __name__ == "__main__":
    main()
