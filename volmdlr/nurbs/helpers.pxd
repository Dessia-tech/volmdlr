# helpers.pxd
from libcpp.vector cimport vector

cdef vector[double] linspace(double start, double stop, int num, int decimals)

cdef double binomial_coefficient(int k, int i)
