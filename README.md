# Animal Breeding Linear Least-Squares Problems

This is a copy of Markus Hegland's animal breeding linear least-squares
problem test set, originally posted on the [CERFACS ftp server](ftp://ftp.cerfacs.fr/pub/animal).

Each problem is an overdetermined linear least-squares problem with a rank
deficiency of 1.
The problems are named according to their relative size.
They all concern an animal breeding concept.
We refer to the documentation for more information, under `original/descr.ps`.

A few changes to the Fortran source files that generate the data were necessary
to get them to compile.
This repository supplies updated Fortran source files along with the data files.
The data is generated in [Harwell-Boeing format](http://math.nist.gov/MatrixMarket/collections/hb.html) by default.

## What is in this repository?

* Updated Fortran source files under `original/Conv`
* the problems in Harwell-Boeing format under `hb`
* the same problems in [Rutherford-Boeing format](https://www.cise.ufl.edu/research/sparse/matrices/DOC/rb.pdf) under `rb`
* solutions, also in Rutherford-Boeing format under `mls`.

I pre-generated the problems for you.
However, if you would like to generate them again yourself, you will need to
edit one of `original/Conv/conv.f` or `original/Conv/conv2.f`, uncomment the
section corresponding to the problem size you are interested in, and compile
the program.
Any recent version of `gfortran` should be sufficient, no external dependency
is required.

The tar archive containing the data files for the extreme-size problems appears
to be cut off and I was not able to un-archive it or repair it.
If you know of a way to repair it, or better yet, have done it, please consider
submitting a pull request!

## Solutions

Because each problem is rank deficient, the solution we elected to report
is the minimum least-squares solution, i.e., among all the least-squares
solutions, the one that has minimum Euclidean norm.
Each solution was generated by way of a full orthogonal decomposition of the
tall and thin problem matrix using Tim Davis' [`Factorize`](https://www.mathworks.com/matlabcentral/fileexchange/24119-don-t-let-that-inv-go-past-your-eyes--to-solve-that-system--factorize-) Matlab code.
The problems in Harwell-Boeing format are read using [hb_to_msm.m](https://people.sc.fsu.edu/~jburkardt/m_src/hb_to_msm/hb_to_msm.html).

Each minimum least-squares solution was generated using the Matlab script
`write_mls.m` found under `matlab`.
My computer runs out of memory when trying to generate the minimum least-squares
solution for `very2`.
I have 16Gb of RAM.
If you have more RAM and are able to generate the solution, please consider
submitting a pull request!

It is important to note that each solution corresponds to a *scaled* problem,
where each column of the matrix is scaled by its Euclidean norm if it is
nonzero.
For this reason, the solutions are named, e.g., `small_scaled_mls.rb`.

The Matlab code writes the solutions as simple vectors in a text file.
The matrices, right-hand sides and solutions are subsequently converted to Rutherford-Boeing format using the Julia package [`HarwellRutherfordBoeing.jl`](https://github.com/JuliaSparse/HarwellRutherfordBoeing.jl).

## References

Here are the original references for this test set.
The second one is included in this repository, under `original/descr.ps`.

```bibtex
@Inbook{hegland-1990,
  author = {Hegland, M.},
  editor = {Burkhart, H.},
  title = {On the computation of breeding values},
  bookTitle = {CONPAR 90---VAPP IV: Joint International Conference on Vector and Parallel Processing Zurich, Switzerland, September 10--13, 1990 Proceedings},
  year = {1990},
  publisher = {Springer Berlin Heidelberg},
  address = {Berlin, Heidelberg},
  pages = {232--242},
  isbn = {978-3-540-46597-3},
  doi = {10.1007/3-540-53065-7_103},
}

@TechReport{hegland-1993,
  author = {Hegland, M.},
  title = {Description and Use of Animal Breeding Data for Large Least Squares Problems},
  institution = {CERFACS},
  year = {1993},
  type = {Technical Report},
  number = {TR/PA/93/50},
  address = {Toulouse, France},
}
```
