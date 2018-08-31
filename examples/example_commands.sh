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

# This example does not use permuted networks, but these commands show how to
# generate them.
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
                -s  "$i" \
                -o  $intermediate/"$network"_"$score"/scores_"$i".tsv
        done
    done
done

################################################################################
#
#   Construct hierarchies.
#
################################################################################

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

################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

# This example uses -lsb/--lower_size_bound 1 because it is a small toy example
# with 25 vertices.  Use larger value (default is 10) for larger graphs.
for network in network_1
do
    for score in scores_1 scores_2
    do
        python src/process_hierarchies.py \
            -oelf $intermediate/"$network"_"$score"/hierarchy_edge_list_0.tsv \
            -oigf $intermediate/"$network"_"$score"/hierarchy_index_gene_0.tsv \
            -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv "; done) \
            -lsb  1 \
            -cf   $results/clusters_"$network"_"$score".tsv \
            -pl   $network $score \
            -pf   $results/sizes_"$network"_"$score".pdf
    done
done

################################################################################
#
#   Perform consensus.
#
################################################################################

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
