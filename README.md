CRISPR SpCas9 gRNA Design Pipeline  
Author: Jai Singh Chowdhary


Overview
-  This project implements a modular Python pipeline for identifying candidate guide RNAs (gRNAs) compatible with the CRISPR SpCas9 system from a user-provided FASTA file.
-  The script scans a selected FASTA record for NGG PAM sites, extracts upstream 20-nt guide sequences, computes GC content for each candidate, filters guides based on 
   biologically motivated GC thresholds, and generates both tabular and visual outputs.
-  The pipeline is designed to be reproducible, configurable, and extensible, following common bioinformatics scripting practices.


Biological Assumptions
-  Cas variant: SpCas9
-  PAM sequence: NGG
-  gRNA length: 20 nucleotides
-  Strand normalization: Sequences may be reverse-complemented to ensure consistent upstream PAM scanning logic
-  GC content filter: 40–60% (proxy for gRNA binding stability)
Note: These assumptions reflect commonly documented and used CRISPR–SpCas9 guide design practices.

Pipeline Steps
-  Load FASTA records using Biopython
-  Select a target record (index-based)
-  Normalize sequence orientation (optional reverse complement)
-  Scan for NGG PAM sites
-  Extract upstream gRNA candidates
-  Compute GC% for each candidate
-  Filter guides by GC content
-  Write viable guides to a TSV file
-  Generate a GC% distribution histogram with cutoff annotations

Requirements
-  Python 3.8+
-  Biopython
-  matplotlib

Usage Instructions 
-  Basic run:
      python grna_design.py --fasta TP53.fna
-  Specify record index:
      python grna_design.py --fasta TP53.fna --record-index 1
-  Disable reverse complement:
      python grna_design.py --fasta TP53.fna --no-rc
-  Customize output filenames:
      python grna_design.py --fasta TP53.fna --out-tsv tp53_guides.tsv --out-png tp53_gc.png
-  View help menu:
      python grna_design.py --help


Outputs

- TSV file (gRNA_doc.tsv by default)
    Contains PAM position, gRNA sequence, and GC percentage for guides passing GC filtering.
- PNG histogram (gc_distribution.png by default)
    Displays GC% distribution of all candidate guides with threshold annotations.



Current Limitations

- FASTA records are selected by index rather than accession ID
- Sequence retrieval requires a local FASTA file
- Only SpCas9 (NGG) PAMs are currently supported

Future Improvements

- Accession-based record selection
- Automated sequence retrieval from public databases (e.g., NCBI, Ensembl)
- Support for additional Cas variants and PAM motifs
- Off-target analysis and additional gRNA quality metrics


Author Notes

This project was developed to explore CRISPR–Cas9 guide RNA design from a computational perspective and to demonstrate biological understanding informed by literature review alongside foundational programming skills. The gRNA design problem was chosen because it requires translating well-established biological constraints into explicit algorithmic logic, including sequence scanning, filtering, and quantitative feature calculation. From a software perspective, the pipeline emphasizes structured programming, modular function design, and clear data flow between analysis stages. The resulting script serves both as a functional gRNA candidate identification tool and as a demonstration of basic bioinformatics pipeline construction.

