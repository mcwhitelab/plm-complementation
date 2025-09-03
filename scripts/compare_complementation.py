import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import sys
from scipy import stats
from io import StringIO

def parse_args():
    parser = argparse.ArgumentParser(
        description='Compare protein sequence identity scores with complementation data',
        epilog='''
Examples:
  # Compare mean embeddings vs SWE vs simple attention
  python compare_complementation.py -e mean_sims.csv swe_sims.csv simple_attn_sims.csv -a align_sims.csv -c paralogs.txt -o results/
  
  # Compare all attention types
  python compare_complementation.py -e mean_sims.csv swe_sims.csv simple_attn_sims.csv self_attn_sims.csv local_self_attn_sims.csv -a align_sims.csv -c paralogs.txt -o results/
  
  # Compare just two embedding types
  python compare_complementation.py -e embeddings.csv swe.csv -a align.csv -c paralogs.txt -o results/
        '''
    )
    parser.add_argument('--embedding_files', '-e', nargs='+', required=True,
                      help='Paths to embedding similarity files (can specify multiple)')
    parser.add_argument('--alignment_similarities', '-a', required=True,
                      help='Path to sequence alignment similarities file')
    parser.add_argument('--complements', '-c', required=True,
                      help='Path to yeast_paralog_complement.txt file')
    parser.add_argument('--output', '-o', required=True,
                      help='Output directory for figures')
    parser.add_argument('--log', '-l', default=None,
                      help='Path to output log file (default: output_dir/analysis_log.txt)')
    return parser.parse_args()

# Create a custom logging class to capture output
class Logger:
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, 'w')
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        
    def flush(self):
        self.terminal.flush()
        self.log.flush()
        
    def close(self):
        self.log.close()

def process_similarity_file(similarity_file, complements):
    similarities = pd.read_csv(similarity_file)
    results = []

    # Determine which column name to use based on the file type
    # For sequence alignment file, use 'identity', for others use 'similarity'
    score_column = 'identity' if 'identity' in similarities.columns else 'similarity'

    # Group by yeast protein
    for _, group in complements.groupby('ScUniProt'):
        if len(group) < 2:  # Skip if less than 2 human paralogs
            continue
        
        yeast_id = group['ScUniProt'].iloc[0]
        
        # Get all human proteins and their complementation status
        for _, row in group.iterrows():
            human_id = row['HsUniProt']
            
            # Find identity/similarity score in similarities file
            match = similarities[
                ((similarities['protein1_id'] == yeast_id) & (similarities['protein2_id'] == human_id)) |
                ((similarities['protein2_id'] == yeast_id) & (similarities['protein1_id'] == human_id))
            ]
            
            if len(match) == 0:
                continue
            
            results.append({
                'yeast_id': yeast_id,
                'human_id': human_id,
                'seq_identity': match[score_column].iloc[0],
                'complements': row['CompStatus'] == 'Complement',
                'rank': row['HsOrthoRank']
            })
    
    return pd.DataFrame(results)

def analyze_results(results_df, method_name):

    identity_col = 'seq_identity'
    
    complement_scores = results_df[results_df['complements']][identity_col]
    non_complement_scores = results_df[~results_df['complements']][identity_col]
    
    print(f"\nSummary Statistics for {method_name}:")
    print(f"Average sequence identity for complementing pairs: {complement_scores.mean():.4f}")
    print(f"Average sequence identity for non-complementing pairs: {non_complement_scores.mean():.4f}")
    
    # Count cases where higher sequence identity doesn't predict complementation
    pairs = results_df.groupby('yeast_id').agg({
        identity_col: list,
        'complements': list,
        'human_id': list
    }).reset_index()
    
    counter_intuitive = 0
    complement_higher = 0
    non_complement_higher = 0
    total_pairs = 0
    
    # Create a list to store pair-wise comparison data for tidy output
    tidy_data = []
    
    for _, row in pairs.iterrows():
        # Add the number of paralogs in the group
        num_paralogs = len(row[identity_col])
        
        if len(row[identity_col]) == 2:  # Only consider complete pairs
            total_pairs += 1
            
            # Get indices for complement and non-complement
            comp_idx = 0 if row['complements'][0] else 1
            non_comp_idx = 1 if row['complements'][0] else 0
            
            comp_id = row[identity_col][comp_idx]
            non_comp_id = row[identity_col][non_comp_idx]
            comp_human_id = row['human_id'][comp_idx]
            non_comp_human_id = row['human_id'][non_comp_idx]
            
            # Store data for tidy output
            tidy_data.append({
                'yeast_id': row['yeast_id'],
                'complementing_human_id': comp_human_id,
                'noncomplementing_human_id': non_comp_human_id,
                'complement_seq_identity': comp_id,
                'noncomplement_seq_identity': non_comp_id,
                'complement_higher': comp_id > non_comp_id,
                'method': method_name,
                'num_paralogs': num_paralogs  # Add number of paralogs
            })
            
            if row[identity_col][0] > row[identity_col][1]:
                if row['complements'][0]:
                    complement_higher += 1
                else:
                    non_complement_higher += 1
                    counter_intuitive += 1
            elif row[identity_col][1] > row[identity_col][0]:
                if row['complements'][1]:
                    complement_higher += 1
                else:
                    non_complement_higher += 1
                    counter_intuitive += 1
    
    # Create a DataFrame from the tidy data
    tidy_df = pd.DataFrame(tidy_data)
    
    print(f"Pairs where complement has higher sequence identity: {complement_higher}/{total_pairs}")
    print(f"Pairs where non-complement has higher sequence identity: {non_complement_higher}/{total_pairs}")
    print(f"Cases where higher sequence identity doesn't predict complementation: {counter_intuitive}/{total_pairs}")
    
    print(f"\nDetailed Analysis for {method_name}:")
    print(f"Total number of pairs: {total_pairs}")
    print(f"Number of cases where complement has higher sequence identity: {complement_higher} ({complement_higher/total_pairs*100:.1f}%)")
    print(f"Number of cases where non-complement has higher sequence identity: {non_complement_higher}/{total_pairs} ({non_complement_higher/total_pairs*100:.1f}%)")
    print(f"Average difference (complement - non-complement): {(complement_scores.mean() - non_complement_scores.mean()):.4f}")
    
    # Add statistical significance test
    t_stat, p_value = stats.ttest_ind(complement_scores, non_complement_scores)
    print(f"T-test p-value: {p_value:.4f}")
    
    print("\nCases where higher sequence identity fails to predict complementation:")
    print(f"{'Yeast_ID':<12} {'Human_Complement':<20} {'Human_Non_Complement':<20} {'Complement_Seq_Identity':<20} {'Non_Complement_Seq_Identity':<20} {'Difference':<10}")
    print("-" * 102)
    
    for _, row in pairs.iterrows():
        if len(row[identity_col]) == 2:  # Only consider complete pairs
            # Get indices for complement and non-complement
            comp_idx = 0 if row['complements'][0] else 1
            non_comp_idx = 1 if row['complements'][0] else 0
            
            comp_id = row[identity_col][comp_idx]
            non_comp_id = row[identity_col][non_comp_idx]
            
            if non_comp_id > comp_id:  # Failed prediction case
                diff = non_comp_id - comp_id
                print(f"{row['yeast_id']:<12} {row['human_id'][comp_idx]:<20} {row['human_id'][non_comp_idx]:<20} {comp_id:>19.4f} {non_comp_id:>19.4f} {diff:>9.4f}")
    
    # Create a list to store paralog group analysis data
    paralog_group_data = []
    
    print("\nDetailed Analysis of Paralog Groups:")
    for _, row in pairs.iterrows():
        if len(row[identity_col]) >= 2:  # Consider groups with 2 or more paralogs
            print(f"\nYeast protein: {row['yeast_id']}")
            print(f"Number of paralogs: {len(row[identity_col])}")
            
            # Sort paralogs by sequence identity
            sorted_indices = sorted(range(len(row[identity_col])), 
                                 key=lambda k: row[identity_col][k],
                                 reverse=True)
            
            print(f"{'Human_ID':<20} {'Complements':<12} {'Seq_Identity':<10}")
            print("-" * 42)
            
            # Store the group info for output
            yeast_id = row['yeast_id']
            
            for rank, idx in enumerate(sorted_indices, 1):
                human_id = row['human_id'][idx]
                complements = row['complements'][idx]
                identity = row[identity_col][idx]
                
                print(f"{human_id:<20} {str(complements):<12} {identity:.4f}")
                
                # Add this paralog to our data list
                paralog_group_data.append({
                    'yeast_id': yeast_id,
                    'human_paralog_id': human_id,
                    'complements': complements,
                    'seq_identity_score': identity,
                    'seq_identity_rank': rank,  # 1-based rank (highest identity = 1)
                    'method': method_name,
                    'num_paralogs_in_group': len(row[identity_col])
                })
            
            # Calculate statistics for this group
            comp_ids = [row[identity_col][i] for i in range(len(row[identity_col])) 
                      if row['complements'][i]]
            non_comp_ids = [row[identity_col][i] for i in range(len(row[identity_col])) 
                          if not row['complements'][i]]
            
            if comp_ids and non_comp_ids:
                print(f"\nAvg complement sequence identity: {sum(comp_ids)/len(comp_ids):.4f}")
                print(f"Avg non-complement sequence identity: {sum(non_comp_ids)/len(non_comp_ids):.4f}")
                print(f"Difference (complement - non-complement): {(sum(comp_ids)/len(comp_ids) - sum(non_comp_ids)/len(non_comp_ids)):.4f}")
    
    # Convert the list to a DataFrame
    paralog_groups_df = pd.DataFrame(paralog_group_data) if paralog_group_data else None
    
    return counter_intuitive/total_pairs if total_pairs > 0 else 0, tidy_df, paralog_groups_df

def main():
    args = parse_args()
    
    # Set up logging
    if args.log is None:
        args.log = os.path.join(args.output, "analysis_log.txt")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Set up the logger
    logger = Logger(args.log)
    sys.stdout = logger
    
    print(f"Analysis started - Output will be saved to {args.log}")
    print(f"Command line arguments: {args}")
    
    # Read complement file
    complements = pd.read_table(args.complements, encoding='latin-1')
    
    # Process each similarity file
    all_results = {}
    error_rates = {}
    all_tidy_data = []
    all_paralog_groups = []  # New list to collect paralog group data
    
    # Process the similarity files
    similarity_files = {}
    
    # Add embedding files with descriptive names
    for i, file_path in enumerate(args.embedding_files):
        # Extract method name from filename
        filename = os.path.basename(file_path)
        if 'swe' in filename.lower():
            method_name = 'SWE Embeddings'
        elif 'simpleattention' in filename.lower() or 'simple_attention' in filename.lower():
            method_name = 'Simple Attention Embeddings'
        elif 'localselfattention' in filename.lower() or 'local_self_attention' in filename.lower():
            method_name = 'Local Self Attention Embeddings'
     
        elif 'selfattention' in filename.lower() or 'self_attention' in filename.lower():
            method_name = 'Self Attention Embeddings'

        elif 'meansig' in filename.lower() or 'sigma' in filename.lower():
            method_name = 'Mean+Sigma Embeddings'
        else:
            method_name = f'Embeddings {i+1}'
        
        # Ensure unique names
        counter = 1
        original_name = method_name
        while method_name in similarity_files:
            method_name = f"{original_name} {counter}"
            counter += 1
        
        similarity_files[method_name] = file_path
    
    # Always add alignment similarities
    similarity_files['Sequence Alignment'] = args.alignment_similarities
    
    print(f"\nProcessing {len(similarity_files)} similarity files:")
    for method_name, file_path in similarity_files.items():
        print(f"  {method_name}: {file_path}")
    print()
    
    for method_name, file_path in similarity_files.items():
        if os.path.exists(file_path):
            results_df = process_similarity_file(file_path, complements)
            all_results[method_name] = results_df
            error_rate, tidy_df, paralog_groups_df = analyze_results(results_df, method_name)
            error_rates[method_name] = error_rate
            
            # Add embedding filename to the tidy data
            tidy_df['embedding_filename'] = os.path.basename(file_path)
            all_tidy_data.append(tidy_df)
            
            # Add paralog group data if available
            if paralog_groups_df is not None:
                paralog_groups_df['embedding_filename'] = os.path.basename(file_path)
                all_paralog_groups.append(paralog_groups_df)
        else:
            print(f"Warning: File not found - {file_path}")
    
    # Combine all tidy data and save to file
    if all_tidy_data:
        combined_tidy_data = pd.concat(all_tidy_data)
        
        # Save the tidy data to a CSV file
        tidy_output_path = os.path.join(args.output, "complementation_comparison_tidy.csv")
        combined_tidy_data.to_csv(tidy_output_path, index=False)
        print(f"\nTidy data saved to: {tidy_output_path}")
        
        # Also print the data in a format similar to your example
        for method_name in similarity_files.keys():
            if method_name in all_results:
                method_data = combined_tidy_data[combined_tidy_data['method'] == method_name]
                print(f"\nRaw data pairs for {method_name}:")
                print(f"{'Yeast_ID':<12}{'Complementing_Human_ID':<20}{'Noncomplementing_Human_ID':<25}{'Complement_seq_identity':<20}{'Noncomplement_seq_identity':<20}{'Complement_Higher':<15}{'Num_Paralogs':<15}{'Embedding_File'}")
                print("-" * 130)
                
                for _, row in method_data.iterrows():
                    print(f"{row['yeast_id']:<12}{row['complementing_human_id']:<20}{row['noncomplementing_human_id']:<25}{row['complement_seq_identity']:<20.4f}{row['noncomplement_seq_identity']:<20.4f}{str(row['complement_higher']):<15}{row['num_paralogs']:<15}{row['embedding_filename']}")
    
    # Combine all paralog group data and save to file
    if all_paralog_groups:
        combined_paralog_groups = pd.concat(all_paralog_groups)
        
        # Save the paralog group data to a CSV file
        paralog_groups_output_path = os.path.join(args.output, "paralog_groups_data.csv")
        combined_paralog_groups.to_csv(paralog_groups_output_path, index=False)
        print(f"\nParalog group data saved to: {paralog_groups_output_path}")
        
        # Print summary statistics about the paralog groups
        print("\nSummary of Paralog Groups Data:")
        print(f"Total number of yeast genes: {combined_paralog_groups['yeast_id'].nunique()}")
        print(f"Total number of human paralogs: {len(combined_paralog_groups)}")
        print(f"Average number of paralogs per yeast gene: {len(combined_paralog_groups) / combined_paralog_groups['yeast_id'].nunique():.2f}")
        
        # Calculate how often the complementing paralog is ranked highest
        complementing_paralogs = combined_paralog_groups[combined_paralog_groups['complements']]
        top_ranked_complements = complementing_paralogs[complementing_paralogs['seq_identity_rank'] == 1]
        
        print(f"\nCases where a complementing paralog is ranked #1 by sequence identity: {len(top_ranked_complements)}/{len(complementing_paralogs.groupby('yeast_id'))} " + 
              f"({len(top_ranked_complements)/len(complementing_paralogs.groupby('yeast_id'))*100:.1f}%)")
    
    # Create scatter plots for each method
    fig, axes = plt.subplots(1, len(all_results), figsize=(15, 5))
    if len(all_results) == 1:
        axes = [axes]
    
    for (name, df), ax in zip(all_results.items(), axes):
        # Create scatter plot comparing complementing vs non-complementing
        sns.scatterplot(data=df, x='seq_identity', y='rank', 
                       hue='complements', alpha=0.6, ax=ax)
        ax.set_title(f"{name}\nParalog Sequence Identity vs Rank")
        ax.set_xlabel('Sequence Identity Score')
        ax.set_ylabel('Ortholog Rank')
    
    plt.tight_layout()
    plt.savefig(f"{args.output}/seq_identity_scatter_multi.png")
    plt.close()
    
    # Create paired points plot for each method
    fig, axes = plt.subplots(1, len(all_results), figsize=(15, 5))
    if len(all_results) == 1:
        axes = [axes]
        
    for (name, df), ax in zip(all_results.items(), axes):
        print(f"\nRaw data pairs for {name}:")
        print("Yeast_ID\tComplement_seq_identity\tNon-complement_seq_identity\tComplement_Higher")
        
        identity_col = 'seq_identity'
        
        pairs = df.groupby('yeast_id').agg({
            identity_col: list,
            'complements': list
        }).reset_index()
        
        for _, row in pairs.iterrows():
            if len(row[identity_col]) == 2 and len(row['complements']) == 2:
                # Get the identities in the right order
                complement_id = row[identity_col][0] if row['complements'][0] else row[identity_col][1]
                non_complement_id = row[identity_col][1] if row['complements'][0] else row[identity_col][0]
                complement_higher = complement_id > non_complement_id
                print(f"{row['yeast_id']}\t{complement_id:.4f}\t{non_complement_id:.4f}\t{complement_higher}")
                
                # Plot points and connecting line
                x = [0, 1]
                y = [complement_id, non_complement_id]  # Now ordered as complement first, non-complement second
                ax.plot(x, y, 'gray', alpha=0.3)
                ax.scatter([0], [y[0]], color='blue', alpha=0.5, label='Complement' if _ == 0 else "")
                ax.scatter([1], [y[1]], color='red', alpha=0.5, label='Non-complement' if _ == 0 else "")
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Complement', 'Non-complement'])
        ax.set_title(name)
        if _ == 0:  # Only add legend for first plot
            ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{args.output}/paired_seq_identity.png")
    plt.close()
    
    # Create violin plot comparing all methods
    plt.figure(figsize=(15, 6))
    data = []
    labels = []
    for name, df in all_results.items():
        data.extend(df['seq_identity'])
        labels.extend([f"{name}\n{['Non-complement', 'Complement'][int(c)]}" 
                      for c in df['complements']])
    
    sns.violinplot(x=labels, y=data)
    plt.xticks(rotation=45)
    plt.title('Distribution of Sequence Identity Scores by Method and Complementation Status')
    plt.xlabel('')
    plt.ylabel('Sequence Identity Score')
    plt.tight_layout()
    plt.savefig(f"{args.output}/seq_identity_comparison.png")
    plt.close()
    
    # Create bar plot of error rates
    plt.figure(figsize=(10, 6))
    plt.bar(error_rates.keys(), error_rates.values())
    plt.title('Error Rate by Method\n(Cases where higher sequence identity fails to predict complementation)')
    plt.xticks(rotation=45)
    plt.ylabel('Error Rate')
    plt.tight_layout()
    plt.savefig(f"{args.output}/error_rates.png")
    plt.close()

    # Create integrated comparison table
    print("\nIntegrated comparison across methods:")
    integrated_results = {}
    
    # First collect all yeast IDs that have complete pairs in any method
    all_yeast_ids = set()
    for df in all_results.values():
        pairs = df.groupby('yeast_id').filter(lambda x: len(x) == 2)
        all_yeast_ids.update(pairs['yeast_id'].unique())
    
    # Create comparison data for each yeast ID
    print(f"{'Yeast_ID':<12}", end='')
    for method in all_results.keys():
        print(f"{method:<20}", end='')
    print()
    print("-" * (12 + 20 * len(all_results)))
    
    for yeast_id in sorted(all_yeast_ids):
        print(f"{yeast_id:<12}", end='')
        for method, df in all_results.items():
            pair = df[df['yeast_id'] == yeast_id]
            if len(pair) == 2:
                # Get complement and non-complement sequence identities
                comp_id = pair[pair['complements']]['seq_identity'].iloc[0]
                non_comp_id = pair[~pair['complements']]['seq_identity'].iloc[0]
                complement_higher = comp_id > non_comp_id
                print(f"{str(complement_higher):<20}", end='')
            else:
                print(f"{'NA':<20}", end='')
        print()
    
    print(f"\nAnalysis complete. Log saved to: {args.log}")
    
    # Close the logger
    sys.stdout = sys.__stdout__
    logger.close()

if __name__ == "__main__":
    main()

