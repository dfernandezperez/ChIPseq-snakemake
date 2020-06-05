# Pasini's lab ChIP-seq pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.3-brightgreen.svg)](https://snakemake.bitbucket.io)

Snakemake-based ChIP-seq pipeline to be run in our PBS-based HPC using singularity containers. The singularity image that is used to run this pipeline is created from [this](https://github.com/dfernandezperez/Docker/blob/master/ChIP-seq/Dockerfile) docker container.

## Setup

The following files are located inside the folder `configuration`. In this folder you will find the files with raw data paths (`units.tsv`), sample metadata (`samples.tsv`), cluster configuration (`cluster.yaml`) and pipeline parameters -alignment, peak calling...- (`config.yaml`).

### Raw data

Paths to raw data are located in the file `units.tsv`. The file has the following structure:

| sample | lane | fq1 | fq2 |
|--------|------|-----|------|
| name_of_sample | name_of_lane_or_resequencing | path/to/forward.fastq | path/to/reverse.fastq |

* The first field correspond to the sample name. This field has to be the same as the sample name that is specified in the `samples.tsv` file (see below). It is recommended to NOT use underscores in the name of the samples, dashed are prefered. I still don't understand why sometimes I get errors if I use them so before fixing I strongly recommend to use dashes instead.

* The second field corresponds to `lane`. The idea of this field is to group fastq files corresponding to the same sample (or to samples that have to be merged). For example, if 1 sample arrived in 2 different lanes from a PE experiment, in total there will be 4 fastqs (2 forward and 2 reverse). In this case, one should enter the same sample 2 times, putting in the `lane` field the corresponding lanes (lane1 and lane2, for example). Actually one can write any word in this field, the idea is to group fastqs from the same sample. All the entries with the same name in the `sample` field with different `lane` will be merged in the same fastq. Here an example of how it would be with 1 sample that arrived in 2 lanes:

| sample | lane | fq1 | fq2 |
|--------|------|-----|------|
| foo | lane1 | path/to/forward_lane1.fastq | path/to/reverse_lane1.fastq |
| foo | lane2 | path/to/forward_lane2.fastq | path/to/reverse_lane2.fastq |

Here I am using lane1 and lane2 for consistency and making things more clear, but the following would also work:

| sample | lane | fq1 | fq2 |
|--------|------|-----|------|
| foo | potato | path/to/forward_lane1.fastq | path/to/reverse_lane1.fastq |
| foo | checazzo | path/to/forward_lane2.fastq | path/to/reverse_lane2.fastq |

* Finally the last 2 fields `fq1` and `fq2` correspond to the paths to the fastq files. `fq1` is the FORWARD read and  `fq2` the REVERSE. The order is very important because they will be sent in that order to the aligner.


### Sample metadata

All metadata and information regarding every sample is located in `samples.tsv`. The file has the following structure:

| NAME | INPUT | SPIKE | AB | USER | GENOME | RUN | IS_INPUT |
|------|-------|-------|----|------|--------|-----|----------|
| name_of_sample | input_to_use | If the sample contains spikein. true or false | antibody | user | versione of genome (i.e: mm10) | run of the sequencing | if the sample is an input |

* For every sample, the `NAME` field has to contain exactly the same name that was written in the `sample` column of the `units.tsv`.

* The `INPUT` field contains the name of the input corresponding to the given sample. It has to be the name of the input written in the fields `sample` and `NAME` from `units.tsv` and  `samples.tsv`.

* `SPIKE`: TRUE or FALSE based on if the sample contains spike-in or not.

* `AB`, `USER`, `RUN`: Metadata corresponding to each sample. If there's nothing to fill I usuallt write an X.

* `GENOME`: Version of the genome used for the alignment. It will be used for peak annotation with ChIPseeker. Right now the accepted values are mm9, mm10, hg19 and hg38.

* `IS_INPUT`: The options are TRUE or FALSE. If the sample is an input set it to TRUE. Also, in case the sample is an input sequenced just to calculate the ratio sample/spike-in that won't be used to call peaks, set it to TRUE.


### Configuration of pipeline parameters

In the root folder of this repository (I say this because there's in another folder a file with the same name) there is the file `config.yaml`. This files contains the configuration of the software and parameters used in the pipeline. Modify them as you wish. Check always that you are using the correct genome files corresponding to the version that you want to use. Also check the effective genome size that is used by deeptools to calculate the GC bias.


### Cluster configuration

`cluster.yaml` contains the per rule cluster parameters (ncpus, ram, walltime...). It can be modified as desired. In the future I want to remove this file in favour of the new [snakemake profiles](https://github.com/Snakemake-Profiles) system (see below), but I still need to understand a little bit better how it works and how to properly do the migration.


### Snakemake profiles

In Snakemake 4.1 [snakemake profiles](https://github.com/Snakemake-Profiles) were introduced. They are supposed to substitute the classic cluster.json file and make the execution of snakemake more simple. The parameters that will be passed to snakemake (i.e: --cluster, --use-singularity...) now are inside a yaml file (`config.yaml`) inside the profile folder (in the case of this repository is `snakemake_profile`). The `config.yaml` inside `snakemake_profile` contains the parameters passed to snakemake. So if you were executing snakemake as `snakemake --cluster qsub --use-singularity` the new `config.yaml` would be like this:

```yaml
cluster: qsub
use-singularity: true
```

## Execution of the pipeline

Once you have all the configuration files as desired, it's time to execute the pipeline. For that you have to execute the `execute_pipeline.sh` script, followed by the name of the rule that you want to execute. If any rule is given it will automatically execute the rule `all` (which would execute the standard pipeline). Examples:

```bash
./execute_pipeline.sh all
```

is equivalent to 

```bash
./execute_pipeline.sh
```

If you want to obtain also broad peaks...

```bash
./execute_pipeline.sh all_broad
```

At the end of the `Snakefile` you will find all the possible target rules and their corresponding output files.


## To Do's

* Migrate 100% to snakemake profiles and stop using the `cluster.yaml` configuration.
* Adding the input-normalization of spike-in samples in case the input is provided.
