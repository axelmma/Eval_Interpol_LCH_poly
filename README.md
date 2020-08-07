# Eval_Interpol_LCH_poly
This is a in implementation in C using the FLINT library [HJP15] for the multipoint evaluation and interpolation of polynomials represented in the LCH-basis over finite fields of characteristic p as described in the paper by Axel Mathieu-Mahias and Michaël Quisquater : 

ISSAC '20: Proceedings of the 45th International Symposium on Symbolic and Algebraic ComputationJuly 2020 Pages 344–351 https://doi.org/10.1145/3373207.3404009.



# Algorithms for multipoint evaluations :

* With higher memory locality of reference
* With lower memory locality of reference (referred to as suboptimal in the program)
* Following the original approach as in [LAHC16].

# Algorithm for multipoint interpolation : 
* With higher memory locality of reference

# Subroutines :
Our algorithms are based on the subroutines for fast standard multipoint evaluation and interpolation of standard polynomials over finite fields of characteristic p.
They can be found in Ch.10 of [GG03].
Also, higher memory locality of reference is achieved by performing perfect shuffles as in [YEMR13]. 


# References : 
[LAHC16] Lin, S., Al-Naffouri, T. Y., Han, Y. S., and Chung, W. (2016). Novel polynomial basis with fast fourier transform and its
application to reed-solomon erasure codes. IEEE Trans. Information Theory, 62(11):6284–6299.

[GG03] Gathen, J. V. Z. and Gerhard, J. (2003). Modern Computer Algebra. Cambridge University Press, New York, NY,
USA, 2 edition.

[YEMR13] Yang, Q., Ellis, J., Mamakani, K., and Ruskey, F. (2013). In-place permuting and perfect shuffling using involutions.
Information Processing Letters, 113(10):386 – 391.

[HJP15] Hart, W., Johansson, F., and Pancratz, S. (2015). FLINT: Fast Library for Number Theory. Version 2.5.2,
http://flintlib.org.
