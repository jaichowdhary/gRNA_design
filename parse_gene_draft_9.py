"""
gRNA Design Project Draft 9

This script scans a selected record from a FASTA file and identifies
candidate guide RNAs (gRNAs) compatible with the CRISPR SpCas9 system.

Biological Assumptions:
    - Cas variant: SpCas9
    - PAM sequence: NGG
    - gRNA length: 20 nucleotides
    - Strand normalization: sequences may be reverse-complemented to ensure
      consistent upstream PAM scanning logic

Pipeline Overview:
    1. Load FASTA records
    2. Select a target record
    3. Normalize sequence orientation
    4. Scan for NGG PAM sites and extract upstream gRNA candidates
    5. Compute GC% for each candidate guide
    6. Filter guides based on GC content (40–60%)
    7. Write viable guides to a TSV file
    8. Generate a GC% distribution histogram with cutoff annotations

Additions in Draft 9:
    - Refactored script into modular, reusable functions
    - Introduced a command-line interface (CLI) using argparse
    - Separated scanning, filtering, output, and visualization logic
    - Improved code readability, configurability, and reproducibility

Current Limitations / Future Improvements:
    - Sequence retrieval currently requires a local FASTA file
    - Record selection is index-based rather than accession-based
    - Only SpCas9 (NGG) PAMs are supported

"""


from Bio import SeqIO
import matplotlib.pyplot as plt 



#CONSTANTS
PAM_SIZE = 3 
GUIDE_LEN = 20
MIN_GC = 40 
MAX_GC = 60 



def main(fasta_path, record_index=0, reverse_complement=True, out_tsv="gRNA_doc.tsv", out_png="gc_distribution.png"):

    file_records = load_records(fasta_path)
    
    selected_record = select_record(file_records, record_index)
    
    normalized_sequence = normalize_plus_strand(selected_record, reverse_complement)

    all_guides = find_candidate_guides(normalized_sequence)
   
    viable_guides = filter_guides_by_gc(all_guides, MIN_GC, MAX_GC)

    guides_to_file(viable_guides, out_tsv)
   
    plot_histogram_and_summary(all_guides, viable_guides, out_png)


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

def normalize_plus_strand(selected_record, reverse_complement=True):
    
    seq = selected_record.seq
    
    #normalizes to the plus strand only if required
    #thus maintaining CRISPR logic
    
    if reverse_complement:
        return seq.reverse_complement()
    return seq



#iterate through the sequence and find "NGG" PAM Sites
#spCas9 PAM = NGG, where "N" can be any nucleotide. 


def find_candidate_guides(normalized_sequence):
    
    all_guides = []
    
    for counter in range(0,len(normalized_sequence)-PAM_SIZE+1) :
       
        if (normalized_sequence[counter+1] == "G") and (normalized_sequence[counter+2]=="G"):
            
            # pam_index corresponds to the first base (N) of the NGG PAM
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
    

def guides_to_file(viable_guides, out_tsv="gRNA_doc.tsv"):
    
    with open(out_tsv, "w") as file_name:
    
        file_name.write("pam_position\tguide_sequence\tgc_percentage\n")
        
        for guide in viable_guides:
        
            file_name.write(
                f"{guide['pam_position']}\t"
                f"{guide['guide_seq']}\t"
                f"{guide['gc_percent']:.2f}\n")

           


def plot_histogram_and_summary(all_guides, viable_guides, out_png="gc_distribution.png"):
   
    
    total_candidates = len(all_guides)
    passing_candidates = len(viable_guides)



    if (total_candidates == 0):
        print("No Candidate Guides were found, thus no plot was generated.")
        return

    else:
        
        #extracting GC content for all selected guides
        all_gc = [guide["gc_percent"] for guide in all_guides]
    
        #printing out the percentage of viable candidates for the user
        
        print("Total candidates found:", total_candidates)
        print(f"Viable candidates post filtering ({MIN_GC}-{MAX_GC}%):", passing_candidates)
       
        passing_percentage = (passing_candidates / total_candidates) * 100
        print(f"Percentage of appropriate candidates: {passing_percentage:.2f}%")
    
        
        #plot histogram for all the guides found, with a filter at the 40-60%
        
        plt.figure()
        plt.hist(all_gc, bins=20, edgecolor="black")
    
        plt.axvline(MIN_GC, color="black", linestyle="--", label=f"Lower GC cutoff ({MIN_GC}%)")
        plt.axvline(MAX_GC, color="black", linestyle="--", label=f"Upper GC cutoff ({MAX_GC}%)")
    
        plt.xlabel("GC Percentage")
        plt.ylabel("Number of Guides")
        plt.title(f"GC% Distribution of Candidate Guides (Filter Window: {MIN_GC}-{MAX_GC}%)")
        plt.legend()
        plt.savefig(out_png, dpi=300)
        plt.show()
                  


#command line interface (CLI): written for reproducible execution

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(
        description="Scan a FASTA record for SpCas9 (NGG) PAMs and extract 20-nt upstream gRNA candidates."
    )

    parser.add_argument("--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("--record-index", type=int, default=0, help="FASTA record index to scan (0-based)")
    parser.add_argument("--no-rc", action="store_true",
                        help="Do NOT reverse-complement the selected record before scanning")
    parser.add_argument("--out-tsv", default="gRNA_doc.tsv", help="Output TSV filename")
    parser.add_argument("--out-png", default="gc_distribution.png", help="Output PNG filename for GC histogram")

    args = parser.parse_args()

    main(
        fasta_path=args.fasta,
        record_index=args.record_index,
        reverse_complement=(not args.no_rc),
        out_tsv=args.out_tsv,
        out_png=args.out_png
    )

