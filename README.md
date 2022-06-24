# HapDup

HapDup (haplotype duplicator) is a pipeline to convert a haploid long read assembly into a diploid assembly.
HapDup works with either Oxford Nanopore or PacBio (HiFi or CLR). The reconstructed haplotypes preserve heterozygous 
structural variants (in addition to small variants) and are locally phased.


## Version 0.8

Input requirements
------------------

HapDup takes as input a haploid long-read assembly, such as produced with [Flye](https://github.com/fenderglass/Flye) or 
[Shasta](https://github.com/chanzuckerberg/shasta). Currenty, ONT reads (Guppy 5+ recommended) and PacBio HiFi
reads are supported.

HapDup is currently designed for low-heterozygosity genomes (such as human). The expectation is that the assembly
has most of the diploid genome collapsed into a single haplotype. For assemblies with partially resolved haplotypes, alternative
alleles could be removed prior to running the pipeline using [purge_dups](https://github.com/dfguan/purge_dups).
We expect to add a better support of highly heterozygous genomes in the future.

The first stage is to realign the original long reads
on the assembly using [minimap2](https://github.com/lh3/minimap2). We recommend to use the latest minimap2 release.

```
minimap2 -ax map-ont -t 30 assembly.fasta reads.fastq | samtools sort -@ 4 -m 4G > lr_mapping.bam
samtools index -@ 4 assembly_lr_mapping.bam
```

Quick start using Docker
------------------------

HapDup is available on the [Docker Hub](https://hub.docker.com/repository/docker/mkolmogo/hapdup).

If Docker is not installed in your system, you need to set it up first following this [guide](https://docs.docker.com/engine/install/ubuntu/).

Next steps assume that your `assembly.fasta` and `lr_mapping.bam` are in the same directory,
which will also be used for HapDup output. If it is not the case, you might need to bind additional 
directories using the Docker's `-v / --volume` argument. The number of threads (`-t` argument)
should be adjusted according to the available resources. For PacBio HiFi input, use
`--rtype hifi` instead of `--rtype ont`.

```
cd directory_with_assembly_and_alignment
HD_DIR=`pwd`
docker run -v $HD_DIR:$HD_DIR -u `id -u`:`id -g` mkolmogo/hapdup:0.8 \
  hapdup --assembly $HD_DIR/assembly.fasta --bam $HD_DIR/lr_mapping.bam --out-dir $HD_DIR/hapdup -t 64 --rtype ont
```

Quick start using Singularity
-----------------------------

Alternatively, you can use [Singularity](https://sylabs.io/guides/3.5/user-guide/). First, you will need install
the client as descibed in the manual. One way to do it is through conda:

```
conda install singularity
```

Next steps assume that your `assembly.fasta` and `lr_mapping.bam` are in the same directory,
which will also be used for HapDup output. If it is not the case, you might need to bind additional 
directories using the `--bind` argument. The number of threads (`-t` argument)
should be adjusted according to the available resources. For PacBio HiFi input, use
`--rtype hifi` instead of `--rtype ont`.

```
singularity pull docker://mkolmogo/hapdup:0.8
HD_DIR=`pwd`
singularity exec --bind $HD_DIR hapdup_0.8.sif \
  hapdup --assembly $HD_DIR/assembly.fasta --bam $HD_DIR/lr_mapping.bam --out-dir $HD_DIR/hapdup -t 64 --rtype ont
```



Output files
------------

The output directory will contain:
* `hapdup_dual_{1,2}.fasta` - dual assembly
* `phased_blocks_hp{1,2}.bed` - phased blocks coordinates (in dual assmeblies)
* `hapdup_phased_{1,2}.fasta` - haplotype-resolved assmebly

Becuase the pipeline is only using long-read (ONT/HiFi) data, it typically does not achieve chromosome-level phasing,
as opposed to the approaches that use Hi-C or Strand-Seq. Hapdup represents the diploid assembly with in two modes: dual and phased.

Dual assembly will have the same contiguity as the original input, but may contain phase switches.
Location of the switches are given in the `bed` files. Dual assembly haplotypes will contain homozogous and heterozygous varinats 
(small and structural). Therefore this type of output is convenient for variant calling.
Read more about this assembly representation [here](http://lh3.github.io/2021/10/10/introducing-dual-assembly).

Alternatively, `hapdup_phased_*` files contain haplotype-resolved contigs. Phases are assigned arbitrary and
different contigs in one file may not be in synch. Because phasing N50 is typically lower than haploid assembly N50,
this representation will require original contigs to be broken and is therefore more fragmented.


Pipeline overview
-----------------

1. HapDup starts with filtering alignments that are likely originating from the unassembled parts of the genome.
Such alignments may later create false haplotypes if not removed (e.g. if reads from a segmental duplication with two copies
can create four haplotypes).

2. Afterwards, PEPPER is used to call SNPs from the filtered alignment file

3. Then we use Margin to phase SNPs and haplotype reads

4. We then use Flye to polish the initiall assembly with the reads from each of the two
haplotypes independently

5. Finally, we find (heterozygous) breakpoints in long-read alignments and apply
the corresponding structural changes to the corresponding polished haplotypes.
Currently, it allows to recover large heterozygous inversions.

Benchmarks
----------

We evaluated HapDup haplotypes in terms of reconstructed structural variants signatures (heterozygous & homozygous)
using the HG002 for which the [curated set of SVs](https://www.nature.com/articles/s41587-020-0538-8) 
is available. We used the [recent ONT data](https://s3-us-west-2.amazonaws.com/miten-hg002/index.html?prefix=guppy_5.0.7/) 
basecalled with Guppy 5.

Given HapDup haplotypes, we called SV using [dipdiff](https://github.com/fenderglass/dipdiff). We also compare SV
set against hifiasm assemblies, even though they were produced from HiFi, rather than ONT reads.
Evaluated using truvari with `-r 2000` option. GT refers to genotype-considered benchmarks.


| Method         | Precision | Recall | F1-score | GT Precision | GT Recall | GT F1-score |
|----------------|-----------|--------|----------|--------------|-----------|-------------|
| Shasta+HapDup  |  0.9500   | 0.9551 | 0.9525   | 0.934        | 0.9543    |  0.9405     |
| Sniffles       |  0.9294   | 0.9143 | 0.9219   | 0.8284       | 0.9051    |  0.8605     |
| CuteSV         |  0.9324   | 0.9428 | 0.9376   | 0.9119       | 0.9416    |  0.9265     |
| hifiasm        |  0.9512   | 0.9734 | 0.9622   | 0.9129       | 0.9723    |  0.9417     |

Yak k-mer based evaluations:

| Hap   |  QV  | Switch err | Hamming err |
|-------|------|------------|-------------|
|     1 |  35  |   0.0389   |   0.1862    |  
|     2 |  35  |   0.0385   |   0.1845    |

Given a minimap2 alignment, HapDup runs in ~400 CPUh and uses ~80 Gb of RAM.

Source installation
-------------------

If you prefer, you can install from source as follows:

```
#create a new conda environemnt and activate it
conda create -n hapdup python=3.8
conda activate hapdup

#get HapDup source
git clone https://github.com/fenderglass/hapdup
cd hapdup
git submodule update --init --recursive

#build and install Flye
pushd submodules/Flye/ && python setup.py install && popd

#build and install Margin
pushd submodules/margin/ && mkdir build && cd build && cmake .. && make && cp ./margin $CONDA_PREFIX/bin/ && popd

#build and install PEPPER and its dependencies
pushd submodules/pepper/ && python -m pip install . && popd
```

To run, ensure that the conda environemnt is activated and then execute:

```
conda activate hapdup
./hapdup.py --assembly assembly.fasta --bam lr_mapping.bam --out-dir hapdup -t 64 --rtype ont
```

Acknowledgements
----------------

The major parts of the HapDup pipeline are:

* [PEPPER](https://github.com/kishwarshafin/pepper)
* [Margin](https://github.com/UCSC-nanopore-cgl/margin)
* [Flye](https://github.com/fenderglass/Flye)


Authors
-------

The pipeline was developed at [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/), Benedict Paten's lab.

Pipeline code contributors:
* Mikhail Kolmogorov

PEPPER/Margin/Shasta support:
* Kishwar Shafin
* Trevor Pesout
* Paolo Carnevali

Citation
--------

If you use HapDup in your research, the most relevant papers to cite are:

Kishwar Shafin, Trevor Pesout, Pi-Chuan Chang, Maria Nattestad, Alexey Kolesnikov, Sidharth Goel, Gunjan Baid et al. 
"Haplotype-aware variant calling enables high accuracy in nanopore long-reads using deep neural networks." bioRxiv (2021).
[doi:10.1101/2021.03.04.433952](https://doi.org/10.1101/2021.03.04.433952)


Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using Repeat Graphs", Nature Biotechnology, 2019
[doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)

License
-------

HapDup is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.
Other software included in this discrubution is released under either MIT or BSD licenses.


How to get help
---------------
A preferred way report any problems or ask questions is the 
[issue tracker](https://github.com/fenderglass/hapdup/issues). 

In case you prefer personal communication, please contact Mikhail at fenderglass@gmail.com.
