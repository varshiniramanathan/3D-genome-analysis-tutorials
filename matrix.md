

---
title: "Contact maps are matrices"
---

This section will go over how to understand and manipulate matrices, especially in the context of 3D genomic contact maps. You will find yourself wrangling a matrix representation of a contact map if you do anything beyond visualization, so understanding how they work can be very helpful. 

### Understanding contact maps as matrices

In the introduction section, we discussed how Micro-C paired-end reads correspond to ligations of physically close DNA fragments, and how these pairs are binned by genomic position. Thus, each end of the pair represents one end of the interaction, and each bin represents the total interaction between a chunk of nearby genomic positions with another chunk of genomic positions. We represent these interactions in a matrix, where the rows and columns both represent all the genomic bins in the genome. Here is a neat animation by James Jusuf that shows how read pairs are binned and how these bins are represented as a matrix. Notice how the color fills in first close to the diagonal. That's because probabilistically, more reads come from ligations of nearby genomic positions, so as you add more reads you get more distal interactions. 

![A video of Micro-C reads distributing in a contact map matrix](https://github.com/user-attachments/assets/50051161-c8be-4d60-b425-31fb5d31b078)

<video src='ttps://github.com/user-attachments/assets/50051161-c8be-4d60-b425-31fb5d31b078' width=180/>


### Sparse and dense matrix representations

`cooltools` and `cooler` rely a lot on what's called a sparse matrix representation. Without going into details, the big problem with big 3D genomics is that it makes gigantic matrices For example, at 1kb binsize, the 250Mb human chr1 is a 250,000 x 250,000 matrix. However, very little of what we do needs to access this whole matrix at once. Also, as far as current tech goes, even the best contact map has a lot of empty bins. Such matrices, where a lot of the bins are empty, can be represented in a more efficient format that stores only the non-zero values. This is what we call a sparse matrix. On the other hand, a typical row, column grid-formatted matrix is considered a dense matrix. In general, `cooltools` uses the sparse matrix whenever possible and then when it actually comes time to plot something, extract some data, etc. on a smaller region, it switches to the dense matrix. 

---

**Next** -> [P(s) curves and observed/expected](expected.md)
