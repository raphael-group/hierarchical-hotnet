#!/usr/bin/env bash

data=$PWD/data
results=$PWD/results

num_permutations=100

# Compile Fortran routines.
cd ../src
f2py -c fortran_routines.f95 -m fortran_routines > /dev/null
cd ../examples

################################################################################
#
#   Prepare data
#
################################################################################

# Create data and results directories.
for network in network_1
do
    for score in score_1 score_2
    do
        mkdir -p $data/scores/$score/permuted/$network
        mkdir -p $data/networks/$network/permuted
        mkdir -p $data/hierarchies/"$network"_"$score"
        mkdir -p $data/hierarchies/"$network"_"$score"/permuted
        mkdir -p $results
    done
done

cd ../src

# Choose beta parameter.
echo "Choosing beta parameter..."

for network in network_1
do
    python choose_beta.py \
        -i $data/networks/$network/edge_list.tsv \
        -o $data/networks/$network/beta.txt
done

# Create similarity matrix.
echo "Creating similarity matrix..."

for network in network_1
do
    beta=`cat $data/networks/$network/beta.txt`

    python create_similarity_matrix.py \
        -i $data/networks/$network/edge_list.tsv \
        -b $beta \
        -o $data/networks/$network/similarity_matrix.h5
done

# Create permuted networks.
# We do not use permuted networks in this example, but we generate them here to
# test the scripts.
echo "Creating permuted networks..."

for network in network_1
do
    # Preserve connectivity of the original graph.
    for i in `seq 1 4`
    do
        python permute_network.py \
            -i $data/networks/$network/edge_list.tsv \
            -s $i \
            -c \
            -o $data/networks/$network/permuted/edge_list_$i.tsv
    done

    # Do not preserve connectivity of the original graph.
    for i in `seq 1 4`
    do
        python permute_network.py \
            -i $data/networks/$network/edge_list.tsv \
            -s $i \
            -o $data/networks/$network/permuted/edge_list_$i.tsv
    done
done

# Create permuted scores.
# The permuted scores only permute scores between vertices in the network.
echo "Creating permuted scores..."

for network in network_1
do
    for score in score_1 score_2
    do
        for i in `seq $num_permutations`
        do
            python permute_scores.py \
                -i   $data/scores/$score/gene_score.tsv \
                -igf $data/networks/$network/index_gene.tsv \
                -s   $i \
                -o   $data/scores/$score/permuted/$network/gene_score_$i.tsv
        done
    done
done

# Constructing hierarchies.
echo "Constructing hierarchies..."

for network in network_1
do
    for score in score_1 score_2
    do
        python construct_hierarchy.py \
            -smf  $data/networks/$network/similarity_matrix.h5 \
            -igf  $data/networks/$network/index_gene.tsv \
            -gsf  $data/scores/$score/gene_score.tsv \
            -df   $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -digf $data/hierarchies/"$network"_"$score"/index_gene.tsv

        for i in `seq $num_permutations`
        do
            python construct_hierarchy.py \
                -smf  $data/networks/$network/similarity_matrix.h5 \
                -igf  $data/networks/$network/index_gene.tsv \
                -gsf  $data/scores/$score/permuted/$network/gene_score_$i.tsv \
                -df   $data/hierarchies/"$network"_"$score"/permuted/edge_list_$i.tsv \
                -digf $data/hierarchies/"$network"_"$score"/permuted/index_gene_$i.tsv
        done
    done
done

################################################################################
#
#   Process results
#
################################################################################

# Plot statistic on hierarchy.
echo "Plotting statistic on hierarchy..."

for network in network_1
do
    for score in score_1 score_2
    do
        python plot_hierarchy_statistic.py \
            -odf  $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -oigf $data/hierarchies/"$network"_"$score"/index_gene.tsv \
            -pdf  $(for i in `seq $num_permutations`; do echo -n "$data/hierarchies/"$network"_"$score"/permuted/edge_list_$i.tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo -n "$data/hierarchies/"$network"_"$score"/permuted/index_gene_$i.tsv "; done) \
            -nc   1 \
            -l    $network $score \
            -o    $results/statistic_"$network"_"$score".pdf
    done
done

# Cut hierarchy.
echo "Cutting hierarchy..."

for network in network_1
do
    for score in score_1 score_2
    do
        python cut_hierarchy.py \
            -odf  $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -oigf $data/hierarchies/"$network"_"$score"/index_gene.tsv \
            -pdf  $(for i in `seq $num_permutations`; do echo -n "$data/hierarchies/"$network"_"$score"/permuted/edge_list_$i.tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo -n "$data/hierarchies/"$network"_"$score"/permuted/index_gene_$i.tsv "; done) \
            -nc   1 \
            -o    $results/clusters_"$network"_"$score".tsv
    done
done

# Perform consensus.
echo "Performing consensus..."

python perform_consensus.py \
    -cf  $results/clusters_network_1_score_1.tsv $results/clusters_network_1_score_2.tsv \
    -igf $data/networks/network_1/index_gene.tsv $data/networks/network_1/index_gene.tsv \
    -elf $data/networks/network_1/edge_list.tsv $data/networks/network_1/edge_list.tsv \
    -n   network_1 network_1 \
    -s   score_1 score_2 \
    -t   2 \
    -o   $results/consensus.tsv
