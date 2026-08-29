---
title: "P(s) curves and observed-over-expected"
---

### Concept overview: Expected genomic interaction probability


### Computing P(s) curves

Please see this [`cooltools` page](https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html) for a detailed overview. Below is a code snippet condensed from those docs, which will just plot a smoothed P(s) and its derivative. There are some convenience additions for naming and saving expected for multiple coolers that you are analyzing together (which is often the case, ex. for different cell types or conditions that you performed Micro-C on together).

It's always good practice to save your expecteds in a dedicated folder for Micro-C/Hi-C, as they can take a while to compute on-the-fly particularly at small binsizes. Then, if you use consistent file naming, you can always pull those expecteds when you need them. For RCMC, you can compute them on-the-fly no problem.

```

# gets multiple expecteds of many coolers in a dataset

import Pathlib
import matplotlib.pyplot as plt
import cooler 
import cooltools
from cooltools.lib import numutils, plotting
import numpy as np
import pandas as pd
import os
import bioframe

uc_expected_path = PATH_TO_SAVE_EXPECTED
uc_filenames = LIST_OF_FILENAMES
resolution = 2000 # as an example

# we want to avoid loading each cooler multiple times so let's load them all now.
# but if you are very memory-constrained or have many coolers, consider doing all the operations (loading P(s), plotting), together in a single loop and loading each cooler on each iteration
uc_coolers = [cooler.Cooler(f"{clr_path}::/resolutions/{resolution}" for clr_name in uc_filenames]

# just a sample way of saving files in an easily retrievable manner (matched properly to the cooler names)
uc_exp_fnames = [P_s_{Path(f).stem}_{resolution}bp.tsv" for f in uc_filenames]

# utility function for getting chrom arms 
def get_hg38_arms(clr=None, excl_ym=True):
    hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
    hg38_cens = bioframe.fetch_centromeres('hg38')
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

    if clr is not None:
        hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

    # we generally want to exclude the Y chrom and the mitochondrial DNA
    if excl_ym:
        hg38_arms = hg38_arms[~hg38_arms["chrom"].isin(['chrY','chrM'])]

    return hg38_arms

# get expected micro-c if it doesn't exist
cvds = [] 
for expfile, clr in zip(uc_exp_fnames, uc_coolers):

    # skipping the centromeres is a bit cleaner
    hg38_arms = get_hg38_arms(clr, excl_ym=True)

    if not (os.path.isfile(os.path.join(uc_expected_path, expfile))):
        # runs parallelized over 16 cores 
        cvd = cooltools.expected_cis(clr=clr, smooth=True, view_df=hg38_arms, aggregate_smoothed=True, nproc=16)
        cvd.to_csv(os.path.join(uc_expected_path, expfile), index=False, sep='\t')
        cvds.append(cvd)
    else:
        cvd = pd.read_csv(os.path.join(uc_expected_path, expfile), sep='\t')
        cvds.append(cvd)

```

That was most of the effort. in the line `cooltools.expected_cis`, we computed the average contact frequency at each genomic separation. This can be directly plotted as the P(s) curve. 

```
# same imports as above, and same variable definitions 

# plot the u-c together on one plot
#f, axs = plt.subplots(1, 2, figsize=(20, 8))

f, axs = plt.subplots(
    2, 1,
    figsize=(6, 9),
    gridspec_kw={'height_ratios': [3, 1]},
    sharex=True
)

cm = ["#33BBEE", "#CC3311", "#009988"] # sacrosanct hansen lab colors (ikykyk)

for i in range(len(coolers)):
    uc_ps = cvds[i]

    uc_ps.loc[uc_ps['dist_bp'] < 2*resolution] = np.nan # first few bins are unreliable
    uc_s = uc_ps['dist_bp']
    uc_ps = uc_ps['balanced.avg.smoothed.agg']

    # plot P(s) in log-log scale
    axs[0].loglog(uc_s, uc_ps/np.max(uc_ps), label= uc_exp_fnames[i]', color=cm[i], linewidth=2)

    der = np.gradient(np.log(uc_ps),
                  np.log(uc_s)

    axs[1].semilogx(uc_s, der, label= uc_exp_fnames[i], color=cm[i],linewidth=2)

    axs[0].set_ylabel('P(s) (max-normalized)')
    axs[1].set_ylabel('dP(s)/ds')
    axs[1].set_xlabel('s [bp]')

plt.legend(loc='lower right', frameon=False)
plt.show()

```

### Observed-over-expected contact maps 

Many quantitative analyses we want to do on Micro-C maps will require P(s) normalization in some form. This is because if we care about how unusual some sort of genomic interaction is (ex. if an interaction is enriched), we need to be comparing the typical interactions at that genomic distance. Otherwise all your signal will correspond to very close-by genomic elements that come into 3D proximity all the time. Many of the 3D genomics utilities, like those for loop calling and pileups, will do this for you. However, sometimes you may want to see an O/E map for a new analysis of yours, or simply just to visualize the difference.

Cooltools can divide a contact map matrix by its expected:

However, this requires storing the whole matrix in memory, which can quickly get memory-infeasible for whole chromosomes. If you want to quickly isolate and normalize regions of interest, use   `cooltools ObsExpSnipper()`. Here is a quick example that computes the O/E intensity of a list of loops. This is not the most computationally efficient way to run this code but it's easy-ish to read.

```
# compute_oe.py

import numpy as np
import os
import pandas as pd
import cooltools
from cooltools.api import snipping

def oe(PS, clr, loops, quantsize):
    loops_curr = []
    PS_df = pd.read_csv(PS, sep='\t')

    # by default, expected-normalization is done with the balanced.avg column (not smoothed or aggregated)
    # if you wish to change this, pass a different expected_value_col as I have done here 
    snipper = snipping.ObsExpSnipper(clr, PS_df, exected_value_col=balanced.avg.smoothed)

    # fetch matrices by chromosome
    for chrom in loops.iloc[:,0].unique():

        # select stores a sparse (memory-light) representation of a big region
        # by default, it's a chromosome in the cooler
        # if you wanted to select another region ex. part of a chrom, you need to set those regions as a view_df and pass it to ObsExpSnipper() above
        mat_chrom = snipper.select(chrom,
                                   chrom)
        loops_chrom = loops[loops['chr1']==chrom] 
     
        # for each loop, snip the loop center +/- quantsize
        for i, loop in loops_chrom.iterrows():

            left_coord = int(np.mean([loop['start1'], loop.iloc['end1']]))
            right_coord = innt(np.mean([loop['start2'], loop.iloc['end2']]))

            # represents a square of size 2*quantsize around the loop center 
            tup = (left_coord - int(quantsize), left_coord + int(quantsize),
                   right_coord - int(quantsize), right_coord + int(quantsize))

            # snip gets a dense (memory-heavy) representation of a small area (tup)
            snip = snipper.snip(mat_chrom, chrom, chrom, tup)
            score = np.average(np.average(snip))

            loops_curr.append(score)

    return loops_curr


if __name__=='__main__':
    # parameters
    resolution = some_integer
    quantsize = some_integer

    # lists
    microc_cooler_names = []
    PS_names = []

    # filepaths
    output_path = ''
    loopcall_path = ''

    microc_coolers = [utils.get_clr(clr_name, resolution) for clr_name in microc_cooler_names]

    for p, clr in zip(PS_names, microc_coolers):
        if not os.path.exists(p):
            print(f"making expected at {p}...")

            # note that this P(s) is not subset to chromosomal arms
            PS_df = cooltools.expected_cis(clr=clr, aggregate_smoothed=True, nproc=8)
            PS_df.to_csv(p, sep='\t', index=False)

    loops = pd.read_csv(loopcall_path, sep='\t', names=['chr1','start1','end1','chr2','start2','end2'])
    loops.sort_values(by=['chr1,'start1','start2'] inplace=True).reset_index(drop=True)

    loops_out = loops.copy()

    for i, cond in enumerate(microc_cooler_names):
        oe_scores_curr = oe(PS_names[i], microc_coolers[i], loops, quantsize)
        loops_out[f'Strength_{cond}'] = oe_scores_curr

```

Again, the comprehensive docs can be found [here for `cooltools.lib.numutils`](https://cooltools.readthedocs.io/en/latest/cooltools.lib.html#module-cooltools.lib.numutils) and [here for `cooltools.api.snipping`](https://cooltools.readthedocs.io/en/latest/cooltools.html#module-cooltools.api.snipping).

