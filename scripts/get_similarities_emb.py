import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import torch
import torch.nn.functional as F
from enum import Enum
import os

# Try to import FAISS, fall back to numpy if not available
try:
    import faiss
    FAISS_AVAILABLE = True
    print("✅ FAISS available - using FAISS for similarity calculations")
except ImportError:
    FAISS_AVAILABLE = False
    print("⚠️  FAISS not available - falling back to numpy arrays")
except Exception as e:
    FAISS_AVAILABLE = False
    print(f"⚠️  FAISS import failed ({e}) - falling back to numpy arrays")

class SimilarityMethod(Enum):
    COSINE = 'cosine'
    ATTENTION = 'attention'
    DIRECT_TEMP = 'direct_temp'
    TOPK = 'topk'
    GAUSSIAN = 'gaussian'
    MULTISCALE = 'multiscale'
    WEIGHTED = 'weighted'
    PCA_COSINE = 'pca_cosine'

def load_seqnames(index_path):
    """Load seqnames and create mapping dictionaries"""
    seqnames_path = Path(str(index_path) + '.seqnames')
    with open(seqnames_path, 'r') as f:
        seqnames = [line.strip() for line in f]
    
    # Create mapping from short ID to full name and vice versa
    short_to_full = {}
    full_to_idx = {}
    for idx, name in enumerate(seqnames):
        short_id = name.split('|')[1]  # e.g., 'P33334' from 'sp|P33334|XXX_YEAST'
        short_to_full[short_id] = name
        full_to_idx[name] = idx
    
    return short_to_full, full_to_idx

def load_embeddings_numpy(index_path):
    """Load embeddings as numpy array (fallback when FAISS not available)"""
    if str(index_path).endswith('.npy'):
        embeddings = np.load(index_path)
    else:
        # Try to load as .npy file
        npy_path = Path(str(index_path) + '.npy')
        if npy_path.exists():
            embeddings = np.load(npy_path)
        else:
            raise FileNotFoundError(f"Could not find embeddings file: {index_path} or {npy_path}")
    
    # Embeddings are already normalized when added to the index
    # No need to normalize again
    return embeddings

def get_embedding_numpy(embeddings, idx):
    """Get embedding at index from numpy array"""
    return embeddings[idx].reshape(1, -1)

def attention_similarity(emb1, emb2, temperature=1.0):
    """Original attention-based similarity."""
    raw_scores = torch.matmul(emb1, emb2.transpose(-2, -1))
    attention_scores = raw_scores / temperature
    attention_weights = F.softmax(attention_scores, dim=-1)
    similarity = torch.sum(attention_weights * raw_scores)
    return similarity.item()

def direct_temp_similarity(emb1, emb2, temperature=1.0):
    """Higher temperature = more emphasis on larger differences"""
    raw_scores = torch.matmul(emb1, emb2.transpose(-2, -1))
    return torch.mean(torch.pow(raw_scores, 1/temperature)).item()

def topk_similarity(emb1, emb2, k=10):
    """Only consider top K strongest interactions"""
    raw_scores = torch.matmul(emb1, emb2.transpose(-2, -1))
    top_k_values, _ = torch.topk(raw_scores.flatten(), min(k, raw_scores.numel()))
    return torch.mean(top_k_values).item()

def gaussian_similarity(emb1, emb2, sigma1=None, sigma2=None, default_sigma=1.0):
    """Use Gaussian kernel with controllable width using per-protein sigma values
    Args:
        emb1: First embedding
        emb2: Second embedding
        sigma1: Sigma values for first protein (optional)
        sigma2: Sigma values for second protein (optional)
        default_sigma: Default sigma to use if specific values aren't provided
    """
    diff = emb1 - emb2
    squared_diff = diff * diff
    
    if sigma1 is not None and sigma2 is not None:
        # Use a proper multivariate Gaussian kernel where each dimension 
        # is scaled by its corresponding sigma values
        # For numerical stability, use max to prevent division by very small numbers
        combined_sigma_squared = torch.max(sigma1 * sigma2, torch.ones_like(sigma1) * 1e-8)
        scaled_squared_diff = squared_diff / combined_sigma_squared
        similarity = torch.exp(-0.5 * torch.sum(scaled_squared_diff))
    else:
        # Fall back to scalar sigma
        sigma_squared = default_sigma * default_sigma
        similarity = torch.exp(-torch.sum(squared_diff) / (2 * sigma_squared))
        
    return similarity.item()

def multiscale_similarity(emb1, emb2, scales=[0.1, 1.0, 10.0]):
    """Compare at multiple scales/temperatures"""
    similarities = []
    for scale in scales:
        diff = emb1 - emb2
        sim = torch.exp(-torch.sum(diff * diff) / (2 * scale * scale))
        similarities.append(sim)
    return torch.mean(torch.stack(similarities)).item()

def weighted_similarity(emb1, emb2, temperature=1.0):
    """Weight dimensions by their magnitude"""
    weights = torch.mean(torch.abs(torch.cat([emb1, emb2])), dim=0)
    weights = F.softmax(weights / temperature, dim=0)
    return torch.sum(weights * (emb1 * emb2)).item()

def calculate_similarities(pairs_file, index_path, output_path, method=SimilarityMethod.COSINE, 
                         temperature=1.0, k=10, sigma=1.0, scales=[0.1, 1.0, 10.0], pca_components=50, sigma_index_path=None):
    
    # Load embeddings based on availability
    if FAISS_AVAILABLE:
        # Use FAISS
        index = faiss.read_index(str(index_path))
        print(f"Loaded FAISS index with {index.ntotal} embeddings")
        
        # Load sigma index if provided and method is Gaussian
        sigma_index = None
        if method == SimilarityMethod.GAUSSIAN and sigma_index_path:
            print("Loading sigma index...")
            sigma_index = faiss.read_index(str(sigma_index_path))
            if sigma_index.ntotal != index.ntotal:
                raise ValueError("Sigma index must have same number of entries as embedding index")
        
        # If using PCA, fit PCA on all embeddings first
        if method == SimilarityMethod.PCA_COSINE:
            print("Fitting PCA...")
            n = index.ntotal
            all_embeddings = np.vstack([index.reconstruct(i) for i in range(n)])
            pca = faiss.PCAMatrix(all_embeddings.shape[1], pca_components)
            pca.train(all_embeddings)
    else:
        # Use numpy fallback
        embeddings = load_embeddings_numpy(index_path)
        print(f"Loaded numpy embeddings with {embeddings.shape[0]} embeddings")
        
        # Load sigma embeddings if provided and method is Gaussian
        sigma_embeddings = None
        if method == SimilarityMethod.GAUSSIAN and sigma_index_path:
            print("Loading sigma embeddings...")
            sigma_embeddings = load_embeddings_numpy(sigma_index_path)
            if sigma_embeddings.shape[0] != embeddings.shape[0]:
                raise ValueError("Sigma embeddings must have same number of entries as embedding index")
        
        # If using PCA, fit PCA on all embeddings first
        if method == SimilarityMethod.PCA_COSINE:
            print("Fitting PCA...")
            from sklearn.decomposition import PCA
            pca = PCA(n_components=pca_components)
            pca.fit(embeddings)
    
    # Load mappings
    short_to_full, full_to_idx = load_seqnames(index_path)
    
    # Read pairs file
    pairs_df = pd.read_csv(pairs_file)
    
    # Prepare output data
    results = []
    skipped_count = 0
    
    for _, row in pairs_df.iterrows():
        prot1_short = row['protein1_id']
        prot2_short = row['protein2_id']
        
        try:
            # Get full names and embeddings
            prot1_full = short_to_full[prot1_short]
            prot2_full = short_to_full[prot2_short]
            
            idx1 = full_to_idx[prot1_full]
            idx2 = full_to_idx[prot2_full]
            
            # Get embeddings based on availability
            if FAISS_AVAILABLE:
                emb1 = torch.tensor(index.reconstruct(idx1).reshape(1, -1))
                emb2 = torch.tensor(index.reconstruct(idx2).reshape(1, -1))
            else:
                emb1 = torch.tensor(get_embedding_numpy(embeddings, idx1))
                emb2 = torch.tensor(get_embedding_numpy(embeddings, idx2))
            
            # Calculate similarity based on method
            if method == SimilarityMethod.PCA_COSINE:
                # Apply PCA transformation
                if FAISS_AVAILABLE:
                    emb1_pca = torch.tensor(pca.apply_py(emb1.numpy()))
                    emb2_pca = torch.tensor(pca.apply_py(emb2.numpy()))
                else:
                    emb1_pca = torch.tensor(pca.transform(emb1.numpy()))
                    emb2_pca = torch.tensor(pca.transform(emb2.numpy()))
                similarity = torch.nn.functional.cosine_similarity(emb1_pca, emb2_pca).item()
            elif method == SimilarityMethod.COSINE:
                similarity = torch.nn.functional.cosine_similarity(emb1, emb2).item()
            elif method == SimilarityMethod.ATTENTION:
                similarity = attention_similarity(emb1, emb2, temperature)
            elif method == SimilarityMethod.DIRECT_TEMP:
                similarity = direct_temp_similarity(emb1, emb2, temperature)
            elif method == SimilarityMethod.TOPK:
                similarity = topk_similarity(emb1, emb2, k)
            elif method == SimilarityMethod.GAUSSIAN:
                if FAISS_AVAILABLE and sigma_index:
                    # Get sigma values for both proteins (FAISS)
                    sigma1 = torch.tensor(sigma_index.reconstruct(idx1).reshape(1, -1)).squeeze()
                    sigma2 = torch.tensor(sigma_index.reconstruct(idx2).reshape(1, -1)).squeeze()
                    similarity = gaussian_similarity(emb1.squeeze(), emb2.squeeze(), sigma1, sigma2, default_sigma=sigma)
                elif not FAISS_AVAILABLE and sigma_embeddings is not None:
                    # Get sigma values for both proteins (numpy)
                    sigma1 = torch.tensor(sigma_embeddings[idx1].reshape(1, -1)).squeeze()
                    sigma2 = torch.tensor(sigma_embeddings[idx2].reshape(1, -1)).squeeze()
                    similarity = gaussian_similarity(emb1.squeeze(), emb2.squeeze(), sigma1, sigma2, default_sigma=sigma)
                else:
                    # Fall back to global sigma if no index provided
                    similarity = gaussian_similarity(emb1.squeeze(), emb2.squeeze(), default_sigma=sigma)
            elif method == SimilarityMethod.MULTISCALE:
                similarity = multiscale_similarity(emb1, emb2, scales)
            elif method == SimilarityMethod.WEIGHTED:
                similarity = weighted_similarity(emb1, emb2, temperature)
            
            results.append({
                'protein1_id': prot1_short,
                'protein2_id': prot2_short,
                'protein1_full': prot1_full,
                'protein2_full': prot2_full,
                'similarity': similarity
            })
        except KeyError:
            skipped_count += 1
            continue
    
    print(f"Skipped {skipped_count} pairs due to missing proteins in the index")
    
    # Create and save output DataFrame
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description='Calculate sequence similarities from pairs file using FAISS index')
    parser.add_argument('-p', '--pairs-file', type=str, required=True,
                       help='Path to the pairs CSV file')
    parser.add_argument('-i', '--index-path', type=str, required=True,
                       help='Path to the FAISS index')
    parser.add_argument('-o', '--output-path', type=str, required=True,
                       help='Path to save the output CSV')
    parser.add_argument('-m', '--method', type=str, default='cosine',
                       choices=[m.value for m in SimilarityMethod],
                       help='Similarity method to use')
    parser.add_argument('-t', '--temperature', type=float, default=1.0,
                       help='Temperature parameter for applicable methods')
    parser.add_argument('-k', '--topk', type=int, default=10,
                       help='K value for topk similarity method')
    parser.add_argument('-s', '--sigma', type=float, default=1.0,
                       help='Sigma value for Gaussian similarity')
    parser.add_argument('--scales', type=float, nargs='+', default=[0.1, 1.0, 10.0],
                       help='Scales for multiscale similarity')
    parser.add_argument('--pca-components', type=int, default=50,
                       help='Number of PCA components to use for pca_cosine method')
    parser.add_argument('--sigma-index-path', type=str,
                       help='Path to the sigma index')
    
    args = parser.parse_args()
    
    method = SimilarityMethod(args.method)
    
    calculate_similarities(
        pairs_file=args.pairs_file,
        index_path=args.index_path,
        output_path=args.output_path,
        method=method,
        temperature=args.temperature,
        k=args.topk,
        sigma=args.sigma,
        scales=args.scales,
        pca_components=args.pca_components,
        sigma_index_path=args.sigma_index_path
    )

if __name__ == '__main__':
    main()
