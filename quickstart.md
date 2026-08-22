---
layout: default
title: "Quickstart"
---

# 1. Quickstart

## Understanding your dataset: read depth, informative reads, and resolution

# A brief summary of Micro-C processing and outputs

First, I will give a broad overview of how the data was processed, since it's easy to run a pipeline without understanding what went into it. In our pipeline (Resources:pipeline for Micro-C and RCMC), we use a standard sequencing aligner (bwa-mem2) with some extra flags to account for differences between typical paired-end libraries and Micro-C libraries. Then, we use pairtools to find read pairs where both ends aligned properly and where the read is unique. Finally, we store these uniquely mapped pairs in a .cool file, which allows us to store and access pairs efficiently. 

Since 3D genomics has low coverage compared to typical methods, we usually need to bin pairs from nearby genomic positions together to visualize and analyze data. People often refer to the size of the bin (how many genomic positions are pooled together into one bin) as "resolution," but it's hard to define resolution in the 3D genomics context. So, I will refer to this term as "binsize," and reserve "resolution" for metrics that actually measure the richness or quality of the dataset (more detail below). 

The main file format that we will be working with is called an .mcool file, which stores coolers binned at multiple binsizes together simultaneously. At each binsize, the cooler is "balanced" to account for systematic differences that could affect read counts, such as accessibility biases. Typically, our .mcools contain binsizes from 250bp to 10Mb. 

# Reporting the quality of your data

Often, the first thing you want to do after processing your data is see how good it is. The quickest and most informative way is to look at the `stats` output of `pairtools dedup`, which is the command that finds and removes duplicate pairs that were marked by `pairtools parse2` (see [pairtools docs](https://pairtools.readthedocs.io/en/latest/index.html) for more information. 

First, let's grab the quantities of interest from the stats file into a python dataframe:

```
def grab_metrics(statfile, *args):
    with open(statfile, 'r') as f:
        stats = [line.rsplit('\t', 1) for line in f]
        stats_dict = {
            key: int(value) / 1_000_000_000
            for key, value in stats
            if key in args
        }

    return os.path.basename(statfile).split('.')[0], stats_dict

metric_args = ['total', 'total_mapped','total_unmapped','total_nodups', 'total_dups', 'cis','trans', 'cis_1kb+', 'cis_20kb+']
metrics = [grab_metrics(statfile, *metric_args) for statfile in stat_files]

metric_df = pd.DataFrame.from_dict(dict(metrics), orient='index')

```

Why do we care about these quantities?
1. `total_dups` / `total_mapped` represents your duplicate fraction. Numbers >50% indicate that most of your reads were either sequencing (optical) duplicates, or PCR duplicates, which could be due to over-sequencing your library or reflect a low-complexity library.
2. `total_nodups` represents the number of uniquely mapped reads, which is the quantity that actually goes into the processed data and therefore represents the quantity of downstream analyzable reads.
3. the `cis` and `trans` ratios can be informative of your digestion and crosslinking, which is out of the scope of this tutorial to explain in detail, but in general we want a cis fraction of 60% at minimum.
4. `cis_1kb+` represents the number of "informative" reads, since very close pairs (<1kb) are not informative for 3D genomic analyses. If most of your reads are <1kb, the data may not be of analyzable quality. 

Lastly, I'll touch on the concept of "resolution."

## Visualization with Hi-Glass:

Hi-Glass is a web viewer that supports 3D genomics datasets in addition to 1D tracks. Please follow the "visualization" section of the [Micro-C/RCMC paper](https://www.nature.com/articles/s41596-026-01393-3) to set up a docker instance. (@Hansen Lab, we have a shared instance running. Please see the private docs for instructions.)

We often want to use 1D tracks like .bigwig and .bed files to view alongside our 3D genomics. To add 1D tracks, we use a package called [`clodius`](https://pypi.org/project/clodius/) to aggregate datasets in a zoom-able way, since they are not inherently multi-binsize like .mcools. There is detailed information on the [Hi-Glass page](https://pkerpedjiev.github.io/higlass-docs-beta/data_preparation.html) but here are some snippets:

```
## bedpe file, ex. loops:
clodius aggregate bedpe --chromsizes-filename /path/to/chromsizes loopfile.bedpe --output-file loopfile.multires
higlass-manage ingest loopfile.multires \
	--filetype bed2ddb \
	--datatype 2d-rectangle-domains \
  --project-name your_project_name

#bw
higlass-manage ingest your_bigwig --assembly your_assembly --project-name your_project_name

## bed:
clodius aggregate bedfile --chromsizes-filename path/to/chromsizes bedfile.bed
higlass-manage ingest --project-name project_name bedfile.bed.beddb

```

If you have Hi-Glass running on a server (@Hansen Lab, this applies to our shared instance; see private docs) then be sure to move the aggregated files to the media directory before you run `ingest`, and add the flag `--no-upload` to your `ingest` command. 

## Resources:

The cooler package is extremely well-documented (see Resources:cooler documentation page).

- [Mapping pipeline for Micro-C and RCMC](https://github.com/ahansenlab/MicroC_RCMC_analysis/) \
- [cooler documentation page](https://cooler.readthedocs.io) \
- [pairtools explanation of pair types](https://pairtools.readthedocs.io/en/latest/formats.html) \


## Citation

If you use this tutorial in your work, please cite:

> Author et al. *Tutorial Title*. Year.

---

**Next:** [P(s) curves and observed-over-expected](expected.md)
