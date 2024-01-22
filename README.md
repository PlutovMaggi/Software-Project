# Software-Project

In this project we implemented a version of the Normalized Spectral Clustering algorithm based on [1, 2].

Algorithm 1 - The Normalized Spectral Clustering Algorithm

Given a set of n points X = x1, x2, . . . , xN ∈ Rd the algorithm is:

1: Form the weighted adjacency matrix W from X

2: Compute the normalized graph Laplacian Lnorm

3: Determine k and obtain the largest k eigenvectors u1, . . . , uk of Lnorm

4: Let U ∈ Rn×k be the matrix containing the vectors u1, . . . , uk as columns

5: Form the matrix T ∈ Rn×k from U by renormalizing each of U ’s rows to have unit length, that is set tij = uij /(∑j (uij)^2 )1/2

6: Treating each row of T as a point in Rk, cluster them into k clusters via the K-means algorithm

The full description and calculations are explained in the following  PDF file:
[project.pdf](https://github.com/PlutovMaggi/Software-Project/files/14001252/project.pdf)


Including the following files:
1. spkmeans.py: Python interface. 
2. spkmeans.h: C header file.
3. spkmeans.c: C interface.
4. spkmeansmodule.c: Python C API wrapper.
5. setup.py: The setup file.
6. comp.sh: compilation script.


References:

[1] Andrew Ng, Michael Jordan, and Yair Weiss. On spectral clustering: Analysis and an algorithm.
Advances in neural information processing systems, 14:849–856, 2001.

[2] Ulrike Von Luxburg. A tutorial on spectral clustering. Statistics and computing, 17(4):395–416,
2007.
8
