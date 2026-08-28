---
title: "Overview of common genomics utilities"
---

Before we get into the analysis sections, I want to take a moment to go over some of the tools and packages you will be seeing consistently in (3D) genomics. I won't actually give tutorials on these tools since many exist, but mostly expose you to them so that it becomes clearer what you have to learn if the terms are strange to you. 

### Frames for your data

A very common way that we represent analyzed data is to tabulate it (make a table that both you and your code can read and interpret). You can dress up your tabulated data in different ways, but the core of it that you will have columns that indicate features of the dataset, and rows that indicate samples of the dataset. For example, columns could be chromosome, start of a genomic interval, and end of genomic interval, and rows could be 1 million genomic intervals. Rows could be loops, and columns could be the coordinates and strengths of those loops. Etc. (You could technically flip rows and columns, but don't do this. Almost all genomics wants your data tabulated in this way.

#### python: pandas, polars, bioframe

`pandas` provides the DataFrame, which is essentially a programmable spreadsheet. You will encounter pandas constantly for filtering, joining, grouping, reshaping, and summarizing tabular genomics data. Many of the `cooltools` and `cooler` functions are representing the data as dataframes (or viewframes; see below) under the hood. So, if you have familiarity with `pandas`, you can better understand what's going on in the code that runs those tools. 

`polars` is another DataFrame library that does many of the same things as `pandas`, but is designed to be faster and more memory-efficient for large datasets. The syntax is somewhat different, and I don't use it at all in this tutorial. But it's good to know that it's basically a faster and cleaner `pandas`. 

`bioframe` builds on the DataFrame concept specifically for genomic intervals. It provides functions for things you would otherwise have to implement yourself, such as finding overlaps between genomic regions. The base unit of `bioframe` is called a view/viewframe. It allows you to load specific genomics formats like .bed files, get genomic intervals, get all the chromosomal arms for your genome, etc. See below for a quick example, and see the [bioframe docs](https://bioframe.readthedocs.io/en/latest/) for more information. 

```
import bioframe

peaks = bioframe.read_table("peaks.bed", schema="bed3")
genes = bioframe.read_table("genes.bed", schema="bed3")

overlaps = bioframe.overlap(peaks, genes)
```

#### R: data.frame, GenomicRanges, InteractionSet/GenomicInteractions

The equivalent structure in R to a python dataframe is, well... `data.frame`, which is a base R package. You will also very frequently encounter `tibble`/`dplyr`, which are part of the [tidyverse ecosystem](https://tidyverse.org/).

In genomics specifically, you'll encounter Bioconductor data structures such as GRanges, which are designed specifically for genomic intervals:

```
library(GenomicRanges)

gtf <- import("genes.gtf")
genes <- gtf[gtf$type == "gene"]

peaks <- GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges(
        start = c(100, 500),
        end = c(200, 600)
    )
)

findOverlaps(peaks, genes)
```

Similar to `bioframe`, `GRanges` presumes that your data is structured in a specific order (in `GRanges`, sequence names and sequence ranges) allows you to interact with these more easily. The flip side is that if you have anything that deviates from this structure, `GRanges` can't handle it that easily (aside from putting it all in the 'metadata' of the object). This is particularly a problem for loops, which are basically two genomic ranges put together (one for each side of the loop). To address this, the `InteractionSet` package allows you to define two anchors as an interaction (`GenomicInteractions`) object. 

#### Command-line tools

While high-level programming languages like python have a lot of conveniences, it pays to become familiar with shell scripting to write and edit files quickly from the command line. We often want to look at, modify, or assess some characteristic of a tab-delimited file (like a BED file). Here are a few:

1. `awk`: for data parsing. For example, say you have a BEDPE file of loops (six columns) and want to get a BED file of anchors (three columns). You could do this with `awk 'BEGIN{OFS=FS="\t"} {print $1,$2,$3; print $4,$5,$6}' loops.bedpe > anchors.bed` (OFS and FS define the field separator, assuming your files are tab-delimited as they should be). See the `awk` man page at Resources:awk docs.
3. `grep`: for searching. For example: `grep 'chr1' peaks.bed > chr1_peaks.bed`
4. `sed` for modifications. For example, to remove the replicate suffix from all items in a file: `sed 's/_rep[0-9]*//' samples.txt > samples_clean.txt`
5. `cut` for snipping out things. For example, get the first three columns from a table: `cut -f1-3 table.tsv > output.bed`
6. `tail` and `head` for viewing portions of files (useful since genomics files are often too large to view with `cat` or load into `less`). Ex. to view the first 20 lines of a file: `head -n 20 peaks.bed`.

---

**Next:** [P(s) curves and observed/expected](expected.md)
