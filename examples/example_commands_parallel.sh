#!/usr/bin/env bash

data=$PWD/data
results=$PWD/results

num_permutations=100

# This script uses gnu parallel to parallizes some of the scripts.  If gnu
# parallel is not already installed on your system, then you can install it from
# https://www.gnu.org/software/parallel/.  You can set the number of cores for
# your system by changing the following num_cores variable.

num_cores=8

# Compile Fortran module.
cd ../src
f2py -c fortran_module.f95 -m fortran_module > /dev/null
cd ..

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

# Choose beta parameter.
echo "Choosing beta parameter..."

for network in network_1
do
    python src/choose_beta.py \
        -i $data/networks/$network/edge_list.tsv \
        -o $data/networks/$network/beta.txt
done

# Create similarity matrix.
echo "Creating similarity matrix..."

for network in network_1
do
    beta=`cat $data/networks/$network/beta.txt`

    python src/create_similarity_matrix.py \
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
    parallel -u -j $num_cores --bar \
        python src/permute_network.py \
            -i $data/networks/$network/edge_list.tsv \
            -s {} \
            -c \
            -o $data/networks/$network/permuted/edge_list_{}.tsv \
        ::: `seq 1 4`

    # Do not preserve connectivity of the original graph.
    parallel -u -j $num_cores --bar \
        python src/permute_network.py \
            -i $data/networks/$network/edge_list.tsv \
            -s {} \
            -o $data/networks/$network/permuted/edge_list_{}.tsv \
        ::: `seq 5 8`
done

# Create permuted scores.
# The permuted scores only permute scores between vertices in the network.
echo "Creating permuted scores..."

for network in network_1
do
    for score in score_1 score_2
    do
        python src/find_permutation_bins.py \
            -gsf $data/scores/$score/gene_score.tsv \
            -igf $data/networks/$network/index_gene.tsv \
            -elf $data/networks/$network/edge_list.tsv \
            -ms  100 \
            -o   $data/scores/$score/permuted/$network/gene_score_bin.tsv

        parallel -u -j $num_cores --bar \
            python src/permute_scores.py \
                -i  $data/scores/$score/gene_score.tsv \
                -bf $data/scores/$score/permuted/$network/gene_score_bin.tsv \
                -o  $data/scores/$score/permuted/$network/gene_score_{}.tsv \
            ::: `seq $num_permutations`
    done
done

# Constructing hierarchies.
echo "Constructing hierarchies..."

for network in network_1
do
    for score in score_1 score_2
    do
        python src/construct_hierarchy.py \
            -smf  $data/networks/$network/similarity_matrix.h5 \
            -igf  $data/networks/$network/index_gene.tsv \
            -gsf  $data/scores/$score/gene_score.tsv \
            -helf $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -higf $data/hierarchies/"$network"_"$score"/index_gene.tsv

        parallel -u -j $num_cores --bar \
            python src/construct_hierarchy.py \
                -smf  $data/networks/$network/similarity_matrix.h5 \
                -igf  $data/networks/$network/index_gene.tsv \
                -gsf  $data/scores/$score/permuted/$network/gene_score_{}.tsv \
                -helf $data/hierarchies/"$network"_"$score"/permuted/edge_list_{}.tsv \
                -higf $data/hierarchies/"$network"_"$score"/permuted/index_gene_{}.tsv \
            ::: `seq $num_permutations`
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
        python src/plot_hierarchy_statistic.py \
            -oelf $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -oigf $data/hierarchies/"$network"_"$score"/index_gene.tsv \
            -pelf $(for i in `seq $num_permutations`; do echo -n " $data/hierarchies/"$network"_"$score"/permuted/edge_list_$i.tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo -n " $data/hierarchies/"$network"_"$score"/permuted/index_gene_$i.tsv "; done) \
            -nc   $num_cores \
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

        echo $(for i in `seq $num_permutations`; do echo " $data/hierarchies/"$network"_"$score"/permuted/edge_list_$i.tsv"; done) > $data/hierarchies/"$network"_"$score"/permuted_edge_list_file.txt
        echo $(for i in `seq $num_permutations`; do echo " $data/hierarchies/"$network"_"$score"/permuted/index_gene_$i.tsv"; done) > $data/hierarchies/"$network"_"$score"/permuted_index_gene_file.txt

        python src/choose_cut_file.py \
            -oelf $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -oigf $data/hierarchies/"$network"_"$score"/index_gene.tsv \
            -pelf $data/hierarchies/"$network"_"$score"/permuted_edge_list_file.txt \
            -pigf $data/hierarchies/"$network"_"$score"/permuted_index_gene_file.txt \
            -hf   $data/hierarchies/"$network"_"$score"/height.txt \
            -sf   $data/hierarchies/"$network"_"$score"/statistic.txt \
            -rf   $data/hierarchies/"$network"_"$score"/ratio.txt

        parallel -u -j $num_cores --bar \
            python src/choose_cut_file.py \
                -oelf $data/hierarchies/"$network"_"$score"/permuted/edge_list_{}.tsv \
                -oigf $data/hierarchies/"$network"_"$score"/permuted/index_gene_{}.tsv \
                -pelf $data/hierarchies/"$network"_"$score"/permuted_edge_list_file.txt \
                -pigf $data/hierarchies/"$network"_"$score"/permuted_index_gene_file.txt \
                -hf   $data/hierarchies/"$network"_"$score"/permuted/height_{}.txt \
                -sf   $data/hierarchies/"$network"_"$score"/permuted/statistic_{}.txt \
                -rf   $data/hierarchies/"$network"_"$score"/permuted/ratio_{}.txt \
            ::: `seq $num_permutations`

    done
done

# Summarize results.
echo "Summarizing results..."
for network in network_1
do
    for score in score_1 score_2
    do

        height=`cat $data/hierarchies/"$network"_"$score"/height.txt`

        python src/cut_hierarchy.py \
            -elf $data/hierarchies/"$network"_"$score"/edge_list.tsv \
            -igf $data/hierarchies/"$network"_"$score"/index_gene.tsv \
            -cc  height \
            -ct  $height \
            -o   $results/"$network"_"$score"_clusters.tsv \

        python src/compute_p_value.py \
            -osf $data/hierarchies/"$network"_"$score"/ratio.txt \
            -psf $(for i in `seq $num_permutations`; do echo -n " $data/hierarchies/"$network"_"$score"/permuted/ratio_$i.txt "; done) \
            -o   $results/"$network"_"$score"_p_value.txt

    done
done

# Perform consensus.
echo "Performing consensus..."

python src/perform_consensus.py \
    -cf  $results/network_1_score_1_clusters.tsv $results/network_1_score_2_clusters.tsv \
    -igf $data/networks/network_1/index_gene.tsv $data/networks/network_1/index_gene.tsv \
    -elf $data/networks/network_1/edge_list.tsv $data/networks/network_1/edge_list.tsv \
    -n   network_1 network_1 \
    -s   score_1 score_2 \
    -t   2 \
    -o   $results/consensus_nodes.tsv \
    -oo  $results/consensus_edges.tsv
