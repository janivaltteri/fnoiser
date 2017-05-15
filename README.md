# fnoiser

A 1/f noise generator package for R

Contains currently three functions:

1) fn1(t,g) generates a single vector of coloured noise with 
given length t and colour g

2) fn2(t,g,r) generates a t-by-2 array of two coloured noises
with a correlation r

3) fnn(n,t,g,r) generates a t-by-n array of n coloured noises.
(correlation is not currently implemented)

This is currently being developed for educational purposes.

The package works by summation of sine waves, which is a rather 
slow method for generating 1/f noise.

