#!/bin/sh
# -------------
#SBATCH --job-name=embed_9606
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --gres=gpu:2
#SBATCH --mail-type=ALL

source ~/.bashrc
#conda activate /groups/clairemcwhite/envs/core_pkgs
#
f=/groups/clairemcwhite/claire_workspace/github/seqsim_project/static/proteomes/small.fasta

#f=/groups/clairemcwhite/claire_workspace/github/seqsim_project/static/proteomes/smaller.fasta


#o=${f}.prott5xl.4l.sweavg.pkl
#o=${f}.esm2_150M.4l.sweavg.pkl
m1=/groups/clairemcwhite/models/prot_t5_xl_uniref50/
m2=/groups/clairemcwhite/models/esm2_t30_150M_UR50D/ 
m3=/groups/clairemcwhite/models/esm2_t33_650M_UR50D/  
m4=/groups/clairemcwhite/models/esm2_t6_8M_UR50D/
#m5=/groups/clairemcwhite/models/esm2_t36_3B_UR50D/
m6=/groups/clairemcwhite/models/esm2_t12_35M_UR50D/
m7=/groups/clairemcwhite/models/ESMplusplus_small
m8=/groups/clairemcwhite/models/ESMplusplus_large
m9=/groups/clairemcwhite/models/prot_bert_bfd
m10=/groups/clairemcwhite/models/ESMC-300M-Protein-Function-Embedding-Only

#m9=../function_classification/models/go_split_v100_improved_fastlearn
#m10=../function_classification/models/go_split_v100_improved



folder=/groups/clairemcwhite/claire_workspace/github/seqsim_project/

# Configuration: Set the layer number for embeddings
# Change this value to use different layers (e.g., "4l", "8l", "12l", "last")
LAYER_NUM="8l"

#for m in $m1 $m2 $m3 $m4 $m5; do
#for m in  $m7 $m8 $m9 $m10; do
for m in $m7; do
    model_name=$(basename $m)
    output_dir=${folder}/data/embeddings/${model_name}
    mkdir -p $output_dir
    
    o=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.pkl
    
    # First generate the embeddings


    python /groups/clairemcwhite/mcwlab_utils/hf_embed_new.py -m $m -f $f -o $o -ss swe meansig simple_attention  self_attention local_self_attention  -s -l -1 -2 -3 -4 -5 -6 -7 -8   -b  1
    # Check if embedding generation was successful
    if [ ! -f "$o" ]; then
        echo "Error: Embedding file $o was not created"
        continue
    fi

    # Process each embedding type separately
    for key in sequence_embeddings sequence_embeddings_sigma sequence_embeddings_swe sequence_embeddings_simpleattention sequence_embeddings_localselfattention; do
        index_file=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.${key}.npy

        rm -f $index_file
        rm -f ${index_file}.seqnames
        rm -f ${index_file}.sources
        python /groups/clairemcwhite/mcwlab_utils/build_incremental_index.py -p $o -o $index_file -k $key

        # Check if index was created successfully
        if [ ! -f "$index_file" ]; then
            echo "Error: Index file $index_file was not created"
            continue
        fi

        # Calculate similarities immediately after building each index
        pairs=/groups/clairemcwhite/claire_workspace/github/seqsim_project/static/orthologies/inp9_559292_9606_prot.pairs 
        similarities=${index_file}.similarities
        #for method in cosine attention direct_temp topk gaussian multiscale weighted; do
        for method in cosine; do

            echo "METHOD: $method"

            python /groups/clairemcwhite/claire_workspace/github/seqsim_project/scripts/get_similarities_emb.py \
                -p "$pairs" \
                -i "$index_file" \
                -o "${similarities}_${method}" \
                -m "$method" \
                --pca-components 50
            echo $similarities
        done
    done

    # Now do the complementation comparison after all similarities are calculated
    #for method in cosine attention direct_temp topk gaussian multiscale weighted; do
    for method in cosine; do
        # Define similarity file paths
        emb_similarities=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.sequence_embeddings.npy.similarities_${method}
        swe_similarities=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.sequence_embeddings_swe.npy.similarities_${method}
        simple_attn_similarities=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.sequence_embeddings_simpleattention.npy.similarities_${method}
        self_attn_similarities=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.sequence_embeddings_selfattention.npy.similarities_${method}
        local_self_attn_similarities=${output_dir}/small.fasta.${LAYER_NUM}.${model_name}.sequence_embeddings_localselfattention.npy.similarities_${method}
        align_similarities=${folder}/data/similarities/${pairs##*/}.align.similarities
        echo "METHOD: $method"
        # Check if required files exist
        if [ ! -f "$emb_similarities" ] || [ ! -f "$swe_similarities" ] || [ ! -f "$simple_attn_similarities" ]  || [ ! -f "$local_self_attn_similarities" ] || [ ! -f "$align_similarities" ]; then
            echo "Error: Missing required similarity files for complementation comparison"
            echo "Looking for:"
            echo $emb_similarities
            echo $swe_similarities
            echo $simple_attn_similarities
            echo $self_attn_similarities
            echo $local_self_attn_similarities
            echo $align_similarities
            continue
        fi

        mkdir -p ${output_dir}/complementation_${model_name}

        # Run complementation comparison for all embedding types
        python /groups/clairemcwhite/claire_workspace/github/seqsim_project/scripts/compare_complementation.py \
            -e "$emb_similarities" "$swe_similarities" "$simple_attn_similarities" "$local_self_attn_similarities" \
            -a "$align_similarities" \
            -c "${folder}/static/external/laurent_garge/combined_paralogs.txt" \
            -o "${output_dir}/complementation_${model_name}"

            #-c "${folder}/static/external/yeast_paralog_complement.txt" \
    done

done
