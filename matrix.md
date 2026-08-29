---
title: "Contact maps are matrices"
---

This section will go over how to understand and manipulate matrices, especially in the context of 3D genomic contact maps. You will find yourself wrangling a matrix representation of a contact map if you do anything beyond visualization, so understanding how they work can be very helpful. 

### Sparse and dense matrix representations

`cooltools` and `cooler` rely a lot on what's called a sparse matrix representation. Without going into details, the big problem with big 3D genomics is that it makes gigantic matrices For example, at 1kb binsize, the 250Mb human chr1 is a 250,000 x 250,000 matrix. However, very little of what we do needs to access this whole matrix at once. Also, as far as current tech goes, even the best contact map has a lot of empty bins. Such matrices, where a lot of the bins are empty, can be represented in a more efficient format that stores only the non-zero values. This is what we call a sparse matrix. On the other hand, a typical row, column grid-formatted matrix is considered a dense matrix. In general, `cooltools` uses the sparse matrix whenever possible and then when it actually comes time to plot something, extract some data, etc. on a smaller region, it switches to the dense matrix. 
