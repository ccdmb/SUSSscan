# SUSSscan
Scan fungal proteome FASTA for identifiable effectors and predicted SUSS structural families based on RemEff clusters from the publication below:
  ```
  Remote homology clustering identifies lowly conserved families of effector proteins in plant-pathogenic fungi
  DAB Jones, PM Moolhuijzen, JK Hane
  Microbial Genomics 7 (9), 000637
  https://doi.org/10.6084/m9.figshare.32859416
  ```
Outputs: CSV and GFF3

SUSSscan.py enables rapid identification of known effectors and prediction of SUSS (sequence-unrelated structurally-similar) effector families, for further validation with more cpu-intensive structural bioinformatic methods.  The RemEff dataset contains 3 levels of protein clusters containing at least 1 known effector.  Level 1 = proteins clustered on similar sequence.  Level 2 = HMM-HMM clusters.  Level 3 = combination of Hidden-Markov, Greedy and Connected Component clusters of similar level 2 clusters.  SUSSscan adds Level 4 descriptive groupings of SUSS families.

Details of SUSS Level 4 groupings are contained in data/SUSS.csv
Details of selected Pfam and publication-derived HMM models relevant to known effectors are in data/Other.csv 


# Installation
Clone the repository:
   ```bash
   git clone https://github.com/ccdmb/SUSSscan.git
   cd SUSSscan
  ```

run install.sh to download Remeff and Other (selected Pfam and publication-derived HMM) databases from https://doi.org/10.6084/m9.figshare.32859416

# Basic Scan (All Hits):
  ```bash
  python SUSSscan.py -i input_proteins.fasta -o output_results.tsv
  ```


# Filter to Top 5 Hits Per Sequence:
  ```bash
  python SUSSscan.py -i input_proteins.fasta -o output_results.tsv --top-hits 5
  ```


# Run using multiple CPU threads for speed:
  ```bash
  python SUSSscan.py -i input_proteins.fasta -o output_results.tsv --threads 8
  ```



# Usage:
  ```python SUSSscan.py -h
  usage: SUSSscan.py [-h] -i INPUT -o OUTPUT [-d DATABASE] [--other-database OTHER_DATABASE] [-c CSV]
                     [--other-csv OTHER_CSV] [-e EVALUE] [--keep-tblout] [--hmmsearch HMMSEARCH] [--cpu CPU]
                     [--top-hits TOP_HITS]
  
  SUSSscan: Search protein sequences against SUSS/Remeff database and annotate hits.
  
  options:
    -h, --help            show this help message and exit
    -i, --input INPUT     Input protein FASTA file
    -o, --output OUTPUT   Output annotated CSV file
    -d, --database DATABASE
                          SUSS/Remeff hmmsearch database file, prefix, or glob pattern (default: data/Remeff.*.hmm)
    --other-database OTHER_DATABASE
                          Other hmmsearch database file, prefix, or glob pattern (default: data/Other.*.hmm)
    -c, --csv CSV         SUSS reference CSV (default: data/SUSS.csv)
    --other-csv OTHER_CSV
                          Other reference CSV (default: data/Other.csv)
    -e, --evalue EVALUE   E-value cutoff (default: 1e-3)
    --keep-tblout         Keep the intermediate .domtblout file
    --hmmsearch HMMSEARCH
                          Path to hmmsearch executable if not in PATH
    --cpu CPU             Number of CPUs to use for hmmsearch (default: all available CPUs)
    --top-hits TOP_HITS   Only output up to the top N hits per query sequence. Use 0 or negative to output all hits
  ```
