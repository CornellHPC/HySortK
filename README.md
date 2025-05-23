## README

HySortK is a high-performance distributed memory K-mer counter. Our paper "High-Performance Sorting-Based K-mer Counting in Distributed Memory with Flexible Hybrid Parallelism" was accepted to ICPP24 and is now available at [https://dl.acm.org/doi/abs/10.1145/3673038.3673072](https://dl.acm.org/doi/abs/10.1145/3673038.3673072).

There're two options for using HySortK: using it as a standalone program or integrating into existing projects.

### Standalone $k$-mer counter

Here's how to build and use HySortK as a standalone program, like KMC3 or kmerind.

#### Build

To compile, load your MPI + OpenMP environment and use `make standalone`. On NERSC's Perlmutter, the default environment should be enough. 

You can change the compile time parameters in the makefile or with the `make` command. Compile time parameters include $K$ value, upper bound and lower bound, and optimization strategies.

For example, use the following command to count $31-mers$ with frequency in [2, 50]

```sh
make standalone K=31 M=17 L=2 U=50 LOG=2 -j8
```

For more information on compile time parameters, please refer to the makefile.

#### Usage

The default IO module of HySortK needs an index file for the dataset. If you do not have one, You use `samtools` to index the dataset before running the kmer counter. On Perlmutter, for example, run:

```sh
# load samtools
module load spack
spack load samtools
samtools faidx ${PATH_TO_DATASET}
```

To run the kmer counter, use your MPI/slurm runner with general parallelization settings followed by `./hysortk $PATH_TO_DATASE`. For example, to run on perlmutter (in the interactive environment) with 4 cpu nodes, use

```sh
srun -N 4 -n 32 -c 64 --cpu_bind=cores ./hysortk ${PATH_TO_DATASET} ${OUTPUT_DIRECTORY}
```

The last argument is **optional**. If it is provided, all $k$-mers within the frequency range will be written to the output directory in {kmer, frequency} format. Each process will write to a separate file.

When the program finishes, a histogram of the $k$-mer frequency will be printed to the console.

On recent hardware ( such as EPYC processors on Perlmutter ), the number of MPI ranks per node to be set to the number of NUMA domains, or the number of Memory Controllers per node, since kmer counting is highly memory and communication bound instead of compute bound.

### Integration into projects

HySortK is also designed to be used as a library for other projects. The flexibility of parallel settings for HySortK makes it an ideal choice for many genomic projects.

An example of integrating HySortK into a existing project will be uploaded to the ELBA project repository.

#### Usage

To use HySortK as a library, you need to include the header file `hysortk.hpp` in your project. The header file contains the following functions:

```cpp
namespace hysortk {
    std::shared_ptr<DnaBuffer> read_dna_buffer(const std::string& fasta_fname, MPI_Comm comm);
    std::unique_ptr<KmerListS> kmer_count(const DnaBuffer& mydna, MPI_Comm comm);
    void print_kmer_histogram(const KmerListS& kmerlist, MPI_Comm comm);
    void write_output_file(const KmerListS& kmerlist, const std::string& output_dir, MPI_Comm comm);
}
```

All functions and data structures are wrapped in the `hysortk` namespace. Since HySortK define its own data structures, you have two choices when integrating HySortK into your project:

1. Use the data structures defined in HySortK throughout your project
2. Convert the data structures to your project's data structures. If you use the identical data structure representations but wrapped in your project's namespace, you can use `reinterpret_cast` to convert between the two data structures.

The first function reads the DNA buffer from a fasta file. The second function `kmer_count`, is the core function of HySortK. It returns a `KmerListS` object in each process, which contains the $k$-mers and their frequencies. The third function prints the histogram of the $k$-mer frequency. The last function writes the $k$-mers to the output directory.  `main.cpp` in the `standalone` directory provides a miminal example of how to use the provided function. You can choose which functions to use in your project, as the functions are designed to be used independently.

The extension information (ReadId and PosInRead) of each $k$-mer are provided when the `EXTENSION` flag is on. Please refer to the following section if they're required in your project.

#### Build

1. Include the `include` directory of HySortK when compiling your project. More specifically, add `-I${PATH_TO_HYSORTK}/include` to your compile command.

2. Currently HySortK will compile to a single `.o` object file. You can build HySortK by calling `make` in its directory with the chosen compiling parameters. For example, include the following in your project's makefile:

```makefile
YOUR_PROJECT:
    $(MAKE) -C ${PATH_TO_HYSORTK} K=31 M=17 L=2 U=50 LOG=2 -j8
```

3. When linking your project, include the `libhysortk.o` object file in the `obj` directory.

If your project needs the EXTENSION information, you need to pass `EXT=1` when building HySortK with Make. As part of the code in the header file is controlled by the flags, you also need to define the `EXTENSION` flag when compiling **all related files** in your project, even if it is not part of HySortK. 

### Citation

If you find this repo helpful to your work, please consider citing our article:

```bibtex
@inproceedings{10.1145/3673038.3673072,
author = {Li, Yifan and Guidi, Giulia},
title = {High-Performance Sorting-Based K-mer Counting in Distributed Memory with Flexible Hybrid Parallelism},
year = {2024},
isbn = {9798400717932},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3673038.3673072},
doi = {10.1145/3673038.3673072},
abstract = {In generating large quantities of DNA data, high-throughput sequencing technologies require advanced bioinformatics infrastructures for efficient data analysis. k-mer counting, the process of quantifying the frequency of fixed-length k DNA subsequences, is a fundamental step in various bioinformatics pipelines, including genome assembly and protein prediction. Due to the growing volume of data, the scaling of the counting process is critical. In the literature, distributed memory software uses hash tables, which exhibit poor cache friendliness and consume excessive memory. They often also lack support for flexible parallelism, which makes integration into existing bioinformatics pipelines difficult. In this work, we propose HySortK, a highly efficient sorting-based distributed memory k-mer counter. HySortK reduces the communication volume through a carefully designed communication scheme and domain-specific optimization strategies. Furthermore, we introduce an abstract task layer for flexible hybrid parallelism to address load imbalances in different scenarios. HySortK achieves a 2-10 \texttimes{} speedup compared to the GPU baseline on 4 and 8 nodes. Compared to state-of-the-art CPU software, HySortK achieves up to 2 \texttimes{} speedup while reducing peak memory usage by 30\% on 16 nodes. Finally, we integrated HySortK into an existing genome assembly pipeline and achieved up to 1.8 \texttimes{} speedup, proving its flexibility and practicality in real-world scenarios.},
booktitle = {Proceedings of the 53rd International Conference on Parallel Processing},
pages = {919–928},
numpages = {10},
keywords = {Computational Biology, Distributed Memory, Genome Analysis, Parallel Radix Sort, Performance Analysis, k-mer Counting},
location = {Gotland, Sweden},
series = {ICPP '24}
}
```


### Copyright

The following projects have been integrated into or used by HySortK:
- ELBA. See https://github.com/PASSIONLab/ELBA
- PARADIS. 3rd-party implementation  https://github.com/odanivan/simple_paradis/tree/master
- RADULS. https://github.com/refresh-bio/RADULS
