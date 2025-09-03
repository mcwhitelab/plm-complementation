import pandas as pd
import argparse
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
import subprocess
import tempfile
import os
import multiprocessing as mp
from functools import partial
import tqdm

def load_sequences(fasta_path):
    """Load sequences and create mapping dictionary"""
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    
    # Create mapping from short ID to sequence
    short_to_seq = {}
    for record in sequences.values():
        short_id = record.id.split('|')[1]  # e.g., 'P33334' from 'sp|P33334|XXX_YEAST'
        short_to_seq[short_id] = str(record.seq)
    
    return short_to_seq

def calculate_alignment_scores(seq1, seq2):
    """Calculate sequence identity and similarity using MAFFT E-INS-i"""
    # Clean sequences
    seq1 = seq1.upper().strip()
    seq2 = seq2.upper().strip()
    
    # Create temporary FASTA file for MAFFT input
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_in:
        temp_in.write(f">seq1\n{seq1}\n>seq2\n{seq2}\n")
        temp_in_path = temp_in.name
    
    # Create temporary file for MAFFT output
    temp_out_path = temp_in_path + '.aligned'
    
    try:
        # Run MAFFT with E-INS-i mode
        cmd = f"mafft-einsi --quiet --genafpair --maxiterate 1000 --ep 0 --thread 4 {temp_in_path} > {temp_out_path}"
        subprocess.run(cmd, shell=True, check=True)
        
        # Read aligned sequences
        aligned_seqs = list(SeqIO.parse(temp_out_path, "fasta"))
        aligned_seq1 = str(aligned_seqs[0].seq)
        aligned_seq2 = str(aligned_seqs[1].seq)
        
        # Calculate identity and similarity including positions where only one sequence has a gap
        aligned_positions = [(a, b) for a, b in zip(aligned_seq1, aligned_seq2) 
                           if not (a == '-' and b == '-')]  # Only exclude positions where both sequences have gaps
        if not aligned_positions:
            return 0.0, 0.0
            
        matches = sum(a == b for a, b in aligned_positions)
        identity = matches / len(aligned_seq1)
        
        # Calculate normalized similarity using BLOSUM62 matrix
        matrix = substitution_matrices.load("BLOSUM62")
        #max_score = max(matrix[a, a] for a in matrix.alphabet)  # Maximum possible score
        
        # Calculate similarity scores, treating gaps as zero similarity
        #similarity_scores = []
        #for a, b in aligned_positions:
        #    if a == '-' or b == '-':
        #        similarity_scores.append(0)  # Gap penalty
        #    else:
        #        similarity_scores.append(matrix[a, b] / max_score)
                
        #similarity = sum(max(0, score) for score in similarity_scores) / len(aligned_positions)
       
        similar_count = 0
        total_positions = len(aligned_positions)

        for a, b in aligned_positions:
             if a != '-' and b != '-' and matrix[a, b] > 0:
                  similar_count += 1

        # Calculate percent similarity
        similarity = similar_count / len(aligned_seq1) 

 
        return identity, similarity
        
    finally:
        # Clean up temporary files
        os.unlink(temp_in_path)
        if os.path.exists(temp_out_path):
            os.unlink(temp_out_path)

def process_pair(short_to_seq, pair):
    """Process a single pair of sequences"""
    try:
        prot1_short, prot2_short = pair
        
        # Get sequences
        seq1 = short_to_seq[prot1_short]
        seq2 = short_to_seq[prot2_short]
        
        # Calculate identity and similarity
        identity, similarity = calculate_alignment_scores(seq1, seq2)
        
        return {
            'protein1_id': prot1_short,
            'protein2_id': prot2_short,
            'identity': identity,
            'similarity': similarity
        }
        
    except Exception as e:
        print(f"Error processing pair {prot1_short} vs {prot2_short}: {str(e)}")
        return {
            'protein1_id': prot1_short,
            'protein2_id': prot2_short,
            'identity': None,
            'similarity': None
        }

def calculate_similarities(pairs_file, fasta_path, output_path):
    # Load sequences
    print(fasta_path)
    short_to_seq = load_sequences(fasta_path)
    print(len(short_to_seq))
    
    # Read pairs file
    pairs_df = pd.read_csv(pairs_file)
    pairs = list(zip(pairs_df['protein1_id'], pairs_df['protein2_id']))
    
    # Find and report missing proteins
    missing_proteins = set()
    for p1, p2 in pairs:
        if p1 not in short_to_seq:
            missing_proteins.add(p1)
        if p2 not in short_to_seq:
            missing_proteins.add(p2)
    
    if missing_proteins:
        print("\nProteins not found in FASTA file:")
        for prot in sorted(missing_proteins):
            print(f"  {prot}")
    
    # Filter pairs to only include those where both proteins are in the FASTA file
    filtered_pairs = [(p1, p2) for p1, p2 in pairs if p1 in short_to_seq and p2 in short_to_seq]
    
    # Report on filtered pairs
    n_original = len(pairs)
    n_filtered = len(filtered_pairs)
    if n_filtered < n_original:
        print(f"\nFiltered out {n_original - n_filtered} pairs where one or both proteins were not found in the FASTA file")
        print(f"Processing {n_filtered} pairs")
    
    # Set up multiprocessing
    n_cores = max(1, mp.cpu_count() - 1)  # Leave one core free
    print(f"Using {n_cores} CPU cores")
    
    # Create pool and process pairs
    with mp.Pool(n_cores) as pool:
        process_func = partial(process_pair, short_to_seq)
        results = list(tqdm.tqdm(
            pool.imap(process_func, filtered_pairs),
            total=len(filtered_pairs),
            desc="Processing pairs"
        ))
    
    print("\nAlignment calculations complete!")
    
    # Create and save output DataFrame
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description='Calculate sequence similarities from pairs file using sequence alignment')
    parser.add_argument('-p', '--pairs-file', type=str, required=True,
                       help='Path to the pairs CSV file')
    parser.add_argument('-f', '--fasta-path', type=str, required=True,
                       help='Path to the FASTA file containing sequences')
    parser.add_argument('-o', '--output-path', type=str, required=True,
                       help='Path to save the output CSV')
    
    args = parser.parse_args()
    
    calculate_similarities(args.pairs_file, args.fasta_path, args.output_path)

if __name__ == '__main__':
    main()
