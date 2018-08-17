#!/usr/bin/env bash

data=$PWD/data
intermediate=$PWD/intermediate
results=$PWD/results

num_permutations=100

# Hierarchical HotNet is parallelizable, but this script runs each Hierarchical
# HotNet sequentially.  Please see the example_commands_parallel.sh script for a
# parallelized example.

# Compile Fortran module.
cd ../src
f2py -c fortran_module.f95 -m fortran_module > /dev/null
cd ..

################################################################################
#
#   Prepare data.
#
################################################################################

# Create data, intermediate data and results, and results directories.
mkdir -p $data
mkdir -p $intermediate
mkdir -p $results

for network in network_1
do
    mkdir -p $intermediate/"$network"
done

for network in network_1
do
    for score in scores_1 scores_2
    do
        mkdir -p $intermediate/"$network"_"$score"
    done
done

################################################################################
#
#   Construct similarity matrices.
#
################################################################################

# Construct similarity matrices.
echo "Construct similarity matrices..."

for network in network_1
do
    python src/construct_similarity_matrix.py \
        -i   $data/"$network"_edge_list.tsv \
        -o   $intermediate/"$network"/similarity_matrix.h5 \
        -bof $intermediate/"$network"/beta.txt
done

################################################################################
#
#   Permute data.
#
################################################################################

# Permute networks.  We do not use permuted networks in this example, but we
# generate them here to test the network permutation scripts.
echo "Permuting networks..."

for network in network_1
do
    cp $data/"$network"_index_gene.tsv $intermediate/"$network"/index_gene_0.tsv
    cp $data/"$network"_edge_list.tsv $intermediate/"$network"/edge_list_0.tsv

    # Preserve connectivity of the observed graph.
    for i in `seq 1 4`
    do
        python src/permute_network.py \
            -i $intermediate/"$network"/edge_list_0.tsv \
            -s "$i" \
            -c \
            -o $intermediate/"$network"/edge_list_"$i".tsv
    done

    # Do not preserve connectivity of the observed graph.
    for i in `seq 5 8`
    do
        python src/permute_network.py \
            -i $intermediate/"$network"/edge_list_0.tsv \
            -s "$i" \
            -o $intermediate/"$network"/edge_list_"$i".tsv
    done
done

# Permute scores.  The permuted scores only exchange scores between vertices in
# the network.
echo "Permuting scores..."

for network in network_1
do
    for score in scores_1 scores_2
    do
        cp $data/"$score".tsv $intermediate/"$network"_"$score"/scores_0.tsv

        python src/find_permutation_bins.py \
            -gsf $intermediate/"$network"_"$score"/scores_0.tsv \
            -igf $data/"$network"_index_gene.tsv \
            -elf $data/"$network"_edge_list.tsv \
            -ms  1000 \
            -o   $intermediate/"$network"_"$score"/score_bins.tsv

        for i in `seq $num_permutations`
        do
            python src/permute_scores.py \
                -i  $intermediate/"$network"_"$score"/scores_0.tsv \
                -bf $intermediate/"$network"_"$score"/score_bins.tsv \
                -o  $intermediate/"$network"_"$score"/scores_"$i".tsv
        done
    done
done

################################################################################
#
#   Construct, summarize, and cut hierarchies.
#
################################################################################

# Construct hierarchies.
echo "Constructing hierarchies..."

for network in network_1
do
    for score in scores_1 scores_2
    do
        for i in `seq 0 $num_permutations`
        do
            python src/construct_hierarchy.py \
                -smf  $intermediate/"$network"/similarity_matrix.h5 \
                -igf  $data/"$network"_index_gene.tsv \
                -gsf  $intermediate/"$network"_"$score"/scores_"$i".tsv \
                -helf $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv \
                -higf $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv
        done
    done
done

# Summarize hierarchies.
echo "Summarizing hierarchies..."

for network in network_1
do
    for score in scores_1 scores_2
    do
        for i in `seq 0 $num_permutations`
        do
            python src/find_cluster_sizes.py \
                -elf $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv \
                -igf $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv \
                -csf $intermediate/"$network"_"$score"/sizes_"$i".txt
        done
    done
done

for network in network_1
do
    for score in scores_1 scores_2
    do
        python src/summarize_cluster_sizes.py \
            -i     $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/sizes_"$i".txt "; done) \
            -esf   $intermediate/"$network"_"$score"/sizes_expected.txt \
            -minsf $intermediate/"$network"_"$score"/sizes_min.txt \
            -maxsf $intermediate/"$network"_"$score"/sizes_max.txt
    done
done

# Cut hierarchies.
echo "Cutting hierarchies..."

for network in network_1
do
    for score in scores_1 scores_2
    do
        for i in `seq 0 $num_permutations`
        do
            python src/find_cut.py \
                -osf $intermediate/"$network"_"$score"/sizes_"$i".txt \
                -esf $intermediate/"$network"_"$score"/sizes_expected.txt \
                -lsb 1 \
                -hf $intermediate/"$network"_"$score"/height_"$i".txt \
                -rf $intermediate/"$network"_"$score"/ratio_"$i".txt
        done
    done
done

################################################################################
#
#   Process results.
#
################################################################################

# Plot cluster sizes.
echo "Plotting cluster sizes..."

for network in network_1
do
    for score in scores_1 scores_2
    do

        height=`cat $intermediate/"$network"_"$score"/height_0.txt`

        python src/plot_hierarchy_statistics.py \
            -osf   $intermediate/"$network"_"$score"/sizes_0.txt \
            -psf   $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/sizes_"$i".txt "; done) \
            -esf   $intermediate/"$network"_"$score"/sizes_expected.txt \
            -minsf $intermediate/"$network"_"$score"/sizes_min.txt \
            -maxsf $intermediate/"$network"_"$score"/sizes_max.txt \
            -ch    $height \
            -l     $network $score \
            -o    $results/sizes_"$network"_"$score".pdf
    done
done

# Identify clusters.
echo "Identifying clusters..."

for network in network_1
do
    for score in scores_1 scores_2
    do
        height=`cat $intermediate/"$network"_"$score"/height_0.txt`

        python src/cut_hierarchy.py \
            -elf $intermediate/"$network"_"$score"/hierarchy_edge_list_0.tsv \
            -igf $intermediate/"$network"_"$score"/hierarchy_index_gene_0.tsv \
            -cc  height \
            -ct  $height \
            -o   $results/clusters_"$network"_"$score".tsv
    done
done

# Evaluate statistical significance.
echo "Evaluating statistical significance..."

for network in network_1
do
    for score in scores_1 scores_2
    do
        python src/compute_p_value.py \
            -osf $intermediate/"$network"_"$score"/ratio_0.txt \
            -psf $(for i in `seq $num_permutations`; do echo -n " $intermediate/"$network"_"$score"/ratio_"$i".txt "; done) \
            -o   $results/p_value_"$network"_"$score".txt
    done
done

# Perform consensus.
echo "Performing consensus..."

python src/perform_consensus.py \
    -cf  $results/clusters_network_1_scores_1.tsv $results/clusters_network_1_scores_2.tsv \
    -igf $data/network_1_index_gene.tsv $data/network_1_index_gene.tsv \
    -elf $data/network_1_edge_list.tsv $data/network_1_edge_list.tsv \
    -n   network_1 network_1 \
    -s   scores_1 scores_2 \
    -t   2 \
    -cnf $results/consensus_nodes.tsv \
    -cef $results/consensus_edges.tsv
