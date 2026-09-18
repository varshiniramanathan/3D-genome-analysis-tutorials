---
title: "Calling, quantifying, and classifying loops"
---

### Concept overview: Loop calling


### Micro-C calling with Mustache

[Mustache](https://github.com/ay-lab/mustache) is our loop caller of choice for Micro-C. Mustache uses a stacked difference of Gaussians approach, which is a common way to detect edges in images. It's very flexible and has a lot of parameters, which can take some trial to optimize but will often yield very good results with the right parameter set. Please read the README on their repo to understand how to install and run it. Here, I will briefly discuss the main considerations for using Mustache for loop calling. 

1. Binsize (`-res` flag in Mustache)

Recalling some of the information from the [first section](stats_n_vis.md), binsize is how many genomic bp we aggregate together per "pixel" in the contact map. Smaller binsizes let us find smaller features, but if the map is too noisy Mustache will struggle to distinguish true loops from noise. The amount of data in a contact map decays with genomic distance. Assuming 2+ billion uniquely mapped reads:

- For CRE loops that are smaller but often shorter-range, you can use a 2kb or 1kb binsize. \
- For CTCF loops, you can use 2kb-5kb binsizes up to a few Mb in genomic separation. \
- For large interactions like Polycomb-mediated interactions, you may have some luck using 10kb binsizes up to many Mb of separation, but these interactions can be difficult to detect by typical loop calling methods.

2. Mustache parameters

The most useful Mustache parameters to tune (in my opinion) are the sparsity threshold (`-st`), distance (`-d`), and p-threshold (`-pt`). Increasing sparsity and decreasing the p-threshold will both reduce false positives, with the tradeoff of lowering loop detection. I prefer to increase the sparsity threshold from the default 0.88 to 0.92 or 0.95 for 1kb loop calling at high resolution because I only want to retain confident loops in high-signal areas (note that this is the opposite of Mustache's guidance, but will help avoid over-detecting in high-res datasets). I would set `-d` (distance) to be the maximum E-P loop distance you expect for 1kb/2kb bin size, and 2-5 Mb for larger binsizes depending on how deep your sequencing is and how long the loops you're interested in might be. A good rule of thumb is is to set your distance cutoff at the genomic separation where you cannot see interactions with your eye at your target binsize.

The problem is that CRE-CRE loops are small and faint. So to catch them, you need to be right on the line between noise and signal for your datasets, which takes a lot of trial and error. A good way to test this is to pick a small chromosome and run Mustache many times at a range of parameters, ex:

-res: 1kb, 2kb, 5kb \
-st: 0.7, 0.88, 0.92 \
-pt: 0.05, 0.1, 0.2 

A way you could do this in a for loop, while using different distances for different binsizes, would be:

```
#!/bin/bash

res_vals=(1kb 2kb 5kb) # mustache takes both kb and bp ex. 1000 would be interpreted the same as 1kb
st_vals=(0.7 0.88 0.92)
pt_vals=(0.1 0.2 0.4)

for res in "${res_vals[@]}"; do

     case "$res" in
        1kb) dist=200000 ;;
        2kb) dist=500000 ;;
        5kb) dist=2000000 ;;
    esac

    for st in ${st_vals[@]}"; do
        for pt in ${pt_vals[@]}"; do

            python3 ./mustache/mustache/mustache.py  \
                -f input.cool \
                -r "$res" \
                --st "$st" \
                --pt "$pt" \
                -d "$dist" \
                -o "mustache_${res}_${st}_${pt}"
        done
    done
done

```

This generates a lot of loop call files (3x3x3=27 per Micro-C dataset, which is a lot to manually check). If you want to start with a more restricted test, I would do something like: 

-res 1kb -st 0.92 -d 200000 (fine-scale loops <200 kb with stringent sparsity cutoff since 1kb has the smallest loops but also the most noise)
-res 2kb -st 0.88 -d 500000 (CRE (hopefully) loops <500 kb)
-res 5kb -st 0.88 or -st 0.7 -d 2000000 (probably mostly CTCF loops).

Notice that I didn't sweep `-pt` in this restricted test. This is because you can achieve roughly equivalent effects by filtering on the FDR later (see below), which saves some testing since you are only FDR-filtering on your favorite parameter set.

3. Filtering loopcalls

In your loop calling sweep, you want to get to a point where you are detecting most of the loops you want, plus some false positives. Then you can try to scale back the false positives. I tend to filter afterwards on q-value if I notice a lot of false positives. This can be done with the following bash script:

```
#bedpe_qval.sh:

#!/bin/bash
name="$1" # whatever mustache gave you, like loopcalls.tsv
qval="$2" # the threshold you want to use

# modify the print statement to drop any extra columns - you might want to drop column 7 for example 
tail -n +2 "$name" | awk -v q="$qval" '{OFS="\t"} ($7 < q) {print $1, $2, $3, $4, $5, $6, $7}' > "${name%.*}_q${qval}.bedpe"  
```

4. Merging loopcalls. 

To detect E-P loops at high confidence from 1kb/2kb and other loops at 5kb+, you want to run Mustache at different binsizes and filters and then merge those. Often, my 1kb/2kb loop calls are where I'm getting high-confidence CRE interactions, and 5kb is "everything else." Then, you want to merge these loopcalls in order of small -> big. 

Here is a short script that uses pairToPair (bedtools) to do this. You'll want to adjust the "slop" (distance between merge-able items) depending on what binsizes you called at. Then you'll have a final set. Be warned that these loops were called at different binsizes so the intervals will all be different sizes. If you don't want this, you can adjust the interval sizes using awk. 

```
p2p_clean.sh:

#!/bin/bash

f1="$1" # highest priority loops
f2="$2"
f3="$3" # lowest priority loops
odir="$4" # loop directory

# generate basenames for all 3 files 
bn1=$(basename "$f1" ".bedpe")
bn2=$(basename "$f2" ".bedpe")
bn3=$(basename "$f3" ".bedpe")

# compare f2 against f1 and f3 against f2/1
pairToPair -a "${odir}/$f2" -b "${odir}/$f1" -type "notboth" -slop 2500 > "${odir}/${bn2}_not_${bn1}.bedpe"
pairToPair -a "${odir}/$f3" -b "${odir}/$f1" -type "notboth" -slop 2500 > "${odir}/${bn3}_not_${bn1}.bedpe"
pairToPair -a "${odir}/$f3" -b "${odir}/$f2" -type "notboth" -slop 2500 > "${odir}/${bn3}_not_${bn2}.bedpe"
pairToPair -a "${odir}/${bn3}_not_${bn2}.bedpe" -b "${odir}/{bn3}_not_${bn1}.bedpe" -type "both" > "${odir}/${bn3}_only.bedpe"

cat "${odir}/$f1" "${odir}/${bn2}_not_${bn1}.bedpe" "${odir}/${bn3}_only.bedpe" > "${odir}/merged_loops.bedpe"
```

### RCMC loop calling with CHIRON

[CHIRON](https://github.com/ahansenlab/chiron) is the loop caller we made for high-resolution RCMC. Mustache can fail to pick up extremely small loops due to the nature of the difference-of-Gaussians method. To detect loops at many size scales, CHIRON uses a CNN that was trained on Micro-C and RCMC. Unlike Mustache, CHIRON only has the binsize parameter to change, but if your RCMC is high-res enough to use CHIRON over Mustache, you should probably stick with the default 1kb anyway. 

