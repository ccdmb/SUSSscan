import argparse
import subprocess
import os
import pandas as pd
import sys
import tempfile
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="SUSSscan: Search protein sequences against SUSS/Remeff database and annotate hits.")
    parser.add_argument("-i", "--input", required=True, help="Input protein FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output annotated CSV file")
    parser.add_argument("-d", "--database", default=None, help="SUSS/Remeff hmmsearch database file, prefix, or glob pattern (default: data/Remeff.*.hmm)")
    parser.add_argument("--other-database", default=None, help="Other hmmsearch database file, prefix, or glob pattern (default: data/Other.*.hmm)")
    parser.add_argument("-c", "--csv", default=os.path.join("data", "SUSS.csv"), help="SUSS reference CSV (default: data/SUSS.csv)")
    parser.add_argument("--other-csv", default=os.path.join("data", "Other.csv"), help="Other reference CSV (default: data/Other.csv)")
    parser.add_argument("-e", "--evalue", type=float, default=1e-3, help="E-value cutoff (default: 1e-3)")
    parser.add_argument("--keep-tblout", action="store_true", help="Keep the intermediate .domtblout file")
    parser.add_argument("--hmmsearch", default="hmmsearch", help="Path to hmmsearch executable if not in PATH")
    parser.add_argument("--cpu", type=int, default=None, help="Number of CPUs to use for hmmsearch (default: all available CPUs)")
    parser.add_argument("--top-hits", type=int, default=0, help="Only output up to the top N hits per query sequence. Use 0 or negative to output all hits (default: 0).")
    return parser.parse_args()

def run_hmmsearch(input_fasta, database, hmmsearch_cmd, evalue_cutoff, cpu_count):
    # Use a temporary file for .domtblout output
    fd, tblout_path = tempfile.mkstemp(suffix=".domtblout")
    os.close(fd)
    
    cmd = [
        hmmsearch_cmd,
        "-E", str(evalue_cutoff),
        "--domE", str(evalue_cutoff),
        "--cpu", str(cpu_count),
        "--domtblout", tblout_path,
        database,
        input_fasta
    ]
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmsearch:\n{e.stderr.decode('utf-8')}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: hmmsearch executable '{hmmsearch_cmd}' not found. Please ensure it is installed and in your PATH.")
        sys.exit(1)
        
    return tblout_path

def parse_hmmsearch_domtblout(tblout_path):
    hits = []
    
    with open(tblout_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.split()
            if len(parts) >= 22:
                seq_id = parts[0]
                hit_id = parts[3]
                hit_acc = parts[4]
                e_val = parts[6]
                score = parts[13]
                hmm_start = parts[15]
                hmm_end = parts[16]
                start = parts[17]
                end = parts[18]
                
                hits.append({
                    "Query_ID": seq_id,
                    "Start": int(start),
                    "End": int(end),
                    "HMM_Start": int(hmm_start),
                    "HMM_End": int(hmm_end),
                    "Hit_ID": hit_id,
                    "Hit_Accession": hit_acc,
                    "E-value": float(e_val),
                    "Score": float(score)
                })
    
    return pd.DataFrame(hits)

def annotate_all_hits(suss_hits, other_hits, suss_csv_path, other_csv_path):
    # Load SUSS reference CSV
    print(f"Loading SUSS reference annotations from {suss_csv_path}...")
    try:
        suss_ref = pd.read_csv(suss_csv_path, low_memory=False)
    except Exception as e:
        print(f"Error reading SUSS reference CSV: {e}")
        sys.exit(1)
        
    # Load Other reference CSV
    print(f"Loading Other reference annotations from {other_csv_path}...")
    try:
        other_ref = pd.read_csv(other_csv_path, low_memory=False)
    except Exception as e:
        print(f"Error reading Other reference CSV: {e}")
        sys.exit(1)
        
    # Annotate SUSS hits
    suss_annotated = pd.DataFrame()
    if not suss_hits.empty:
        suss_annotated = suss_hits.merge(suss_ref, left_on="Hit_ID", right_on="Level_2_Cluster", how="left")
        suss_annotated["Level4ID"] = suss_annotated.get("Level_4_SUSS_Family")
        suss_annotated["Description"] = suss_annotated.get("Level_4_SUSS_Description")
        
    # Annotate Other hits
    other_annotated = pd.DataFrame()
    if not other_hits.empty:
        # Pfam accessions in domtblout are like PF03032.21, we strip the version.
        # Fall back to Hit_ID if accession is missing/'-'.
        other_hits["Accession_clean"] = other_hits.apply(
            lambda r: r["Hit_Accession"].split('.')[0] if (isinstance(r["Hit_Accession"], str) and r["Hit_Accession"] != "-") else r["Hit_ID"],
            axis=1
        )
        other_annotated = other_hits.merge(other_ref, left_on="Accession_clean", right_on="Accession", how="left")
        other_annotated["Level_4_SUSS_Family"] = other_annotated.get("Level4ID")
        other_annotated["Level_4_SUSS_Description"] = other_annotated.get("Description")
        # Clean up temporary column
        other_annotated = other_annotated.drop(columns=["Accession_clean"])
        
    # Combine the datasets
    if suss_annotated.empty and other_annotated.empty:
        return pd.DataFrame()
    elif suss_annotated.empty:
        combined_df = other_annotated
    elif other_annotated.empty:
        combined_df = suss_annotated
    else:
        combined_df = pd.concat([suss_annotated, other_annotated], ignore_index=True)
        
    # Reorder columns
    cols = combined_df.columns.tolist()
    important_cols = [
        "Query_ID", "Start", "End", "HMM_Start", "HMM_End", "Hit_ID", "E-value", "Score",
        "Level4ID", "Level_4_SUSS_Family", "Suggested_Level_4_SUSS_Family",
        "Level_3_Supercluster", "Description", "Level_4_SUSS_Description"
    ]
    
    ordered_cols = [c for c in important_cols if c in cols]
    for c in cols:
        if c not in ordered_cols:
            ordered_cols.append(c)
            
    return combined_df[ordered_cols]

def get_fasta_ids(fasta_path):
    ids = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                ids.append(seq_id)
    return ids

def write_gff3(df, gff_path):
    print(f"Writing GFF3 results to {gff_path}...")
    try:
        with open(gff_path, 'w', encoding='utf-8') as f:
            f.write("##gff-version 3\n")
            if df.empty:
                return
                
            for idx, row in df.reset_index(drop=True).iterrows():
                seqid = row['Query_ID']
                source = "SUSSscan"
                feature_type = "polypeptide_domain"
                start = row['Start']
                end = row['End']
                score = row['Score']
                strand = "."
                phase = "."
                
                # Construct attributes
                attrs = []
                attrs.append(f"ID=match_{idx+1}")
                attrs.append(f"Name={row['Hit_ID']}")
                
                hit_acc = row.get('Hit_Accession')
                if pd.notna(hit_acc) and hit_acc != "-":
                    clean_acc = hit_acc.split('.')[0]
                    attrs.append(f"Dbxref=Pfam:{clean_acc}")
                    
                hmm_start = row.get('HMM_Start')
                hmm_end = row.get('HMM_End')
                if pd.notna(hmm_start) and pd.notna(hmm_end):
                    attrs.append(f"Target={row['Hit_ID']} {int(hmm_start)} {int(hmm_end)}")
                    
                e_val = row.get('E-value')
                if pd.notna(e_val):
                    attrs.append(f"evalue={e_val}")
                    
                l4_id = row.get('Level4ID')
                if pd.notna(l4_id) and l4_id != "":
                    attrs.append(f"level4ID={l4_id}")
                    
                desc = row.get('Description')
                if pd.notna(desc) and desc != "":
                    clean_desc = str(desc).replace(';', ',').replace('=', ' ').replace('%', '%25')
                    attrs.append(f"Note={clean_desc}")
                    
                attr_str = ";".join(attrs)
                
                f.write(f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attr_str}\n")
    except Exception as e:
        print(f"Error writing GFF3 file: {e}")
        sys.exit(1)

def main():
    args = parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
        
    # Resolve split databases
    if args.database:
        if "*" in args.database or "?" in args.database:
            databases = sorted(glob.glob(args.database))
        elif os.path.isdir(args.database):
            databases = sorted(glob.glob(os.path.join(args.database, "*.hmm")))
        elif os.path.exists(args.database):
            databases = [args.database]
        else:
            # Try matching as prefix
            databases = sorted(glob.glob(args.database + "*.hmm"))
    else:
        databases = sorted(glob.glob(os.path.join("data", "Remeff.*.hmm")))
        if not databases:
            fallback = os.path.join("data", "SUSS.hmm")
            if os.path.exists(fallback):
                databases = [fallback]
                
    if not databases:
        print(f"Error: SUSS database files matching '{args.database or 'data/Remeff.*.hmm'}' not found.")
        sys.exit(1)
        
    # Resolve split Other databases
    if args.other_database:
        if "*" in args.other_database or "?" in args.other_database:
            other_databases = sorted(glob.glob(args.other_database))
        elif os.path.isdir(args.other_database):
            other_databases = sorted(glob.glob(os.path.join(args.other_database, "Other.*.hmm")))
            if not other_databases and os.path.exists(os.path.join(args.other_database, "Other.hmm")):
                other_databases = [os.path.join(args.other_database, "Other.hmm")]
        elif os.path.exists(args.other_database):
            other_databases = [args.other_database]
        else:
            other_databases = sorted(glob.glob(args.other_database + "*.hmm"))
    else:
        other_databases = sorted(glob.glob(os.path.join("data", "Other.*.hmm")))
        if not other_databases:
            fallback = os.path.join("data", "Other.hmm")
            if os.path.exists(fallback):
                other_databases = [fallback]
                
    if not other_databases:
        print(f"Error: Other database files matching '{args.other_database or 'data/Other.*.hmm'}' not found.")
        sys.exit(1)
        
    cpu_count = args.cpu
    if cpu_count is None:
        cpu_count = os.cpu_count()
        if cpu_count is None:
            cpu_count = 1
            
    print(f"Starting SUSSscan with evalue cutoff {args.evalue}...")
    
    print("Searching against SUSS database...")
    suss_hits_list = []
    temp_tblouts = []
    for db in databases:
        print(f"Searching against database partition {db}...")
        tblout = run_hmmsearch(args.input, db, args.hmmsearch, args.evalue, cpu_count)
        temp_tblouts.append(tblout)
        hits = parse_hmmsearch_domtblout(tblout)
        if not hits.empty:
            suss_hits_list.append(hits)
            
    if suss_hits_list:
        suss_hits = pd.concat(suss_hits_list, ignore_index=True)
        if args.top_hits > 0:
            print(f"Filtering SUSS hits to top {args.top_hits} per query sequence...")
            suss_hits = suss_hits.sort_values(['Query_ID', 'E-value'])
            suss_hits = suss_hits.groupby('Query_ID').head(args.top_hits)
    else:
        suss_hits = pd.DataFrame()
    
    print("Searching against Other database...")
    other_hits_list = []
    other_temp_tblouts = []
    for db in other_databases:
        print(f"Searching against database partition {db}...")
        tblout = run_hmmsearch(args.input, db, args.hmmsearch, args.evalue, cpu_count)
        other_temp_tblouts.append(tblout)
        hits = parse_hmmsearch_domtblout(tblout)
        if not hits.empty:
            other_hits_list.append(hits)
            
    if other_hits_list:
        other_hits = pd.concat(other_hits_list, ignore_index=True)
    else:
        other_hits = pd.DataFrame()
    
    # Annotate hits from both databases
    annotated_df = annotate_all_hits(suss_hits, other_hits, args.csv, args.other_csv)
    
    if annotated_df.empty:
        print("No significant hits found.")
    else:
        print("Sorting results...")
        fasta_ids = get_fasta_ids(args.input)
        unique_fasta_ids = list(dict.fromkeys(fasta_ids))
        
        # Query_ID sorting based on input fasta order
        annotated_df['Query_ID'] = pd.Categorical(annotated_df['Query_ID'], categories=unique_fasta_ids, ordered=True)
        
        # Temporary sorting column for Level4ID
        annotated_df['Level4ID_sort'] = annotated_df['Level4ID'].fillna('')
        
        # Sort values first by Query_ID, Level4ID_sort, and sequence coordinates to calculate consecutive overlapping groups
        annotated_df = annotated_df.sort_values(['Query_ID', 'Level4ID_sort', 'Start', 'End'])
        
        # Calculate locus groups (overlapping hits) for each Query_ID and Level4ID
        locus_ids = []
        current_key = None
        locus_counter = 0
        max_end = -1
        
        for idx, row in annotated_df.iterrows():
            key = (row['Query_ID'], row['Level4ID_sort'])
            start = row['Start']
            end = row['End']
            
            if key != current_key:
                current_key = key
                locus_counter = 0
                max_end = end
                locus_ids.append((idx, locus_counter))
            else:
                if start <= max_end:
                    max_end = max(max_end, end)
                    locus_ids.append((idx, locus_counter))
                else:
                    locus_counter += 1
                    max_end = end
                    locus_ids.append((idx, locus_counter))
                    
        locus_df = pd.DataFrame(locus_ids, columns=['index', 'Locus_ID']).set_index('index')
        annotated_df['Locus_ID'] = locus_df['Locus_ID']
        
        # Sort by:
        # 1. Query_ID (fasta order)
        # 2. Level4ID (alphabetical)
        # 3. Locus_ID (consecutive along sequence)
        # 4. E-value (lowest first)
        annotated_df = annotated_df.sort_values(['Query_ID', 'Level4ID_sort', 'Locus_ID', 'E-value'])
        
        # Drop temporary sorting helper columns
        annotated_df = annotated_df.drop(columns=['Level4ID_sort', 'Locus_ID'])
        
        # Reorder Hit_Accession if it exists to place it next to Hit_ID
        if 'Hit_Accession' in annotated_df.columns:
            cols = annotated_df.columns.tolist()
            cols.remove('Hit_Accession')
            hit_id_idx = cols.index('Hit_ID')
            cols.insert(hit_id_idx + 1, 'Hit_Accession')
            annotated_df = annotated_df[cols]
            
        annotated_df.to_csv(args.output, index=False)
        print(f"Success! Annotated results saved to {args.output}")
        
        base_path, _ = os.path.splitext(args.output)
        gff_path = base_path + ".gff3"
        write_gff3(annotated_df, gff_path)
        
    if not args.keep_tblout:
        for tblout in temp_tblouts:
            if os.path.exists(tblout):
                os.remove(tblout)
        for tblout in other_temp_tblouts:
            if os.path.exists(tblout):
                os.remove(tblout)
    else:
        print(f"Raw SUSS hmmsearch outputs kept at: {', '.join(temp_tblouts)}")
        print(f"Raw Other hmmsearch outputs kept at: {', '.join(other_temp_tblouts)}")

if __name__ == "__main__":
    main()
