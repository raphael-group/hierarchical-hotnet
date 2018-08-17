Hierarchical HotNet
=======================

Hierarchical HotNet is an algorithm for finding hierarchies of altered subnetworks.  While originally developed for use with cancer mutation data on protein-protein interaction networks, Hierarchical HotNet supports any application in which scores may be associated with the nodes of a network, i.e., a vertex-weighted graph.

Setup
------------------------
The setup process for Hierarchical HotNet requires the following steps:

### Download
Download Hierarchical HotNet.  The following command clones the current Hierarchical HotNet repository from GitHub:

    git clone https://github.com/raphael-group/hierarchical-hotnet.git

### Installation
The following software is either required for Hierarchical HotNet or optional but recommended for better performance or for visualizations.

##### Required
* Linux/Unix
* [Python (2.7 or 3.5)](http://python.org/)
* [NumPy (1.14)](http://www.numpy.org/)
* [SciPy (0.19)](http://www.scipy.org/)
* [NetworkX (1.11)](http://networkx.github.io/)
* [h5py (2.7)](http://www.h5py.org/)

##### Optional, but recommended
* A Fortran compiler, e.g., [gfortran 5.4](https://gcc.gnu.org/wiki/GFortran)
* [virtualenv](https://virtualenv.pypa.io/en/stable/)
* [GNU parallel](https://www.gnu.org/software/parallel/)
* [Matplotlib (2.1)](http://matplotlib.org/)

Most likely, Hierarchical HotNet will work with other versions of the above software.

In particular, both [virtualenv](https://virtualenv.pypa.io/en/stable/) and [GNU parallel](https://www.gnu.org/software/parallel/) are recommended in practice.  Virtualenv provides a virtual environment that allows Python packages to be installed or updated independently of the system packages.  GNU parallel facilitates running many scripts in parallel.  We highly recommend running Hierarchical HotNet in parallel.

### Compilation
Install a Fortran compiler, such as [gfortran](https://gcc.gnu.org/wiki/GFortran), for better performance.  The following command compiles the optional Fortan code used in Hierarchical HotNet:

    cd src
    f2py -c fortran_module.f95 -m fortran_module > /dev/null

We highly recommend using the Fortran code for better performance.  However, Hierarchical HotNet will transparently fall back to a Python-only implementation if a Fortran compiler is unavailable or if compilation is unsuccessful.

### Testing
To test Hierarchical HotNet on an example network with two sets of example scores, please run the following script:

    sh examples/example_commands.sh

This script illustrates the full Hierarchical HotNet pipeline.  It should require less than a minute or two of CPU time, 100 MB of RAM, and 1 MB of storage space.  If this script runs successfully, then Hierarchical HotNet is ready to use.

Alternatively, to run Hierarchical HotNet in parallel on the sample example data, please run the following script:

    sh examples/example_commands_parallel.sh

It is straightforward to modify the above scripts to run Hierarchical HotNet on a compute cluster.

Use
----------------
Hierarchical HotNet requires the use of several scripts on a few input files.

### Input
There are three input files for Hierarchical HotNet that together define a network with scores on the nodes of the network.  For example, the following example defines a network with an edge between the nodes ABC and DEF, which have scores 0.5 and 0.2, respectively.  For convenience, these files use the same format as the input files for HotNet2.

##### Index-to-gene file
This file associates each gene with an index, which we use for the edge list as well as a similarity matrix:

    1   ABC
    2   DEF

##### Edge list file
This file defines a network using the indices in the index-to-gene file:

    1    2

##### Gene-to-score file
This file associates each gene with a score:

    ABC 0.5
    DEF 0.2

### Running
We provide a script, `hierarchical_hotnet.py`, for running the entire Hierarchical HotNet pipeline.  **This script will be added shortly.**  This script combines the following steps, which can also be run separately (see `examples/example_commands.sh` and `examples/example_commands_parallel.sh` for examples):

1. Create a similarity matrix for each network by running the `src/create_similarity_matrix.py` script.

2. Depending on your choice of null model, generate either permuted scores for each set of scores and network by running the `src/find_permutation_bins.py` and `src/permute_scores.py` scripts or permuted networks for each network by running the `src/permute_networks.py` script.  In general, it is faster to permute scores than networks.

3. Construct hierarchies on the observed network and gene scores as well as the permuted networks and gene scores by running the `src/construct_hierarchies.py` script.

4. Find cuts of the hierarchy by running the `src/find_cluster_sizes.py`, `src/summarize_cluster_sizes.py`, `src/find_cut.py`, and `src/cut_hierarchy.py` scripts, and evaluate the statistical significance of the results by running the `src/compute_p_value.py` script

5. Perform the consensus summarization procedure on the results by running the `src/perform_consensus.py` script.

See `examples/example_commands.sh` or `examples/example_commands_parallel.sh` for full minimal working examples of Hierarchical HotNet that illustrate the use of each of these scripts, including the inputs and outputs for the Hierarchical HotNet pipeline.

The `hierarchical_hotnet.py` script runs Hierarchical HotNet checks for intermediate results, which it will use if available.  For example, if the similarity matrices are already available, then the script can skip this step of Hierarchical HotNet.

### Output
Hierarchical HotNet identifies statistically significant regions of a hierarchical clustering of topologically close, high-scoring genes.  Hierarchical HotNet also performs a consensus across hierarchical clusterings from different networks and gene scores.

Additional information
----------------

### Examples
See the `examples` directory for example data, scripts, and output for Hierarchical HotNet.

### Support
For support with Hierarchical HotNet, please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users).  Please try one of the examples in the `examples` directory before running Hierarchical HotNet with your own data, and please provide any error messages encountered with these examples to expedite troubleshooting.

### License
See `LICENSE.txt` for license information.  **Add license.**

### Citation
If you use Hierarchical HotNet in your work, then please cite the following reference:

> M.A. Reyna, M.D.M. Leiserson, B.J. Raphael. Hierarchical HotNet: identifying hierarchies of altered subnetworks. *Bioinformatics*.  2018.
