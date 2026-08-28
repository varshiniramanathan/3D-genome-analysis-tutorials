---
title: "Introduction"
---

# 1. Introduction

## Prerequisites

This tutorial is mostly comprised of basic python and bash, but doesn't require proficiency in either. It starts from the post-mapping stage and assumes that you used a pipeline similar to the [Hansen Lab's pipeline](https://github.com/ahansenlab/MicroC_RCMC_analysis/) to produce [cooler](https://cooler.readthedocs.io/en/latest/datamodel.html)-format datasets. I'll generally refer to 3D genomics datasets as Micro-C, although for the most part these analyses will also apply to Hi-C and RCMC (except, of course, whole-chromosome and whole-genome analyses do not apply to RCMC). 

While RCMC can usually be handled on a personal computer, whole-genome high-res 3D genomics datasets are large. Many of the code snippets make use of parallelization and assume computing and memory limits that will be difficult to achieve on a personal computer. 
---

**Next:** [Understanding your dataset, stats reporting, and Hi-Glass visualization](stats_n_vis.md)
