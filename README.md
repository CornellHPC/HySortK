## README

HySortK is a high-performance distributed memory K-mer counter. 
Our paper "High-Performance Sorting-Based K-mer Counting in Distributed Memory with Flexible Hybrid Parallelism" was accepted to ICPP24 and will be available online soon!

### Build

To compile, load the MPI + OpenMP environment and use `make`. On NERSC's Perlmutter, the default environment should be enough. 

You can change the compile time parameters in the makefile or with the `make` command. Compile time parameters include $K$ value, upper bound and lower bound, and optimization strategies.

For example, use the following command to count $31-mers$ with frequency in [2, 50]

```sh
make K=31 M=17 L=2 U=50 LOG=2 -j8
```

### Usage

Note: The current IO module of HySortK has not been optimized. It may not perform well with large datasets. Work in progress for related optimization.

Currently, the IO module needs an index file for the dataset. If you do not have one, You use `samtools` to index the dataset before running the kmer counter. On Perlmutter, for example, run:

```sh
# load samtools
module load spack
spack load samtools
# index the dataset
smatools faidx ${PATH_TO_DATASET}
```

To run the kmer counter, use your MPI/slurm runner with general parallelization settings followed by `./hysortk $PATH_TO_DATASE`. For example, to run on perlmutter (in the interactive environment) with 4 cpu nodes, use

```sh
srun -N 4 -n 32 -c 64 --cpu_bind=cores ./hysortk ${PATH_TO_DATASET}
```

Under default log level, the histogram of counted k-mers will be printed at the end of k-mer counting.

### Copyright

The following projects have been integrated into or used by HySortK:
- ELBA. See https://github.com/PASSIONLab/ELBA
- PARADIS. 3rd-party implementation  https://github.com/odanivan/simple_paradis/tree/master
- RADULS. https://github.com/refresh-bio/RADULS
