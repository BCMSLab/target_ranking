# target_ranking
Code for the BMC Genomics paper (Integrating binding and expression data to predict transcription factors combined function).
The article describes the Bioconductor R package [target](http://github.com/MahShaaban/target)

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/cregart/) image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image that can be obtained and launched on any local 
machine running [docker](https://hub.docker.com/r/bcmslab/target/).

```bash
$ docker pull bcmslab/target:latest
$ docker run -it bcmslab/target:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in a form of research compendium. This includes the scripts to reproduce
the figures and tables in this manuscript. From within the container, [git](https://git-scm.com) can be used to clone the source code.

The following code clones the repository containing the source code.

```bash
$ git clone http://github.com/BCMSLab/target_ranking
```

## Generating figures and tables

The script to generate the figures and tables in the manuscirpt can be run through `Makefile`

```bash
$ cd target_ranking
$ make
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 3.6.2 (2019-12-12) on `x86\_64-pc-linux-gnu`.
