#ifndef DEF_LCH
#define DEF_LCH

#include <gmp.h>
#include "math.h"
#include "flint.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_nmod_poly.h"
#include "fmpz.h"
#include "fq_poly.h"
#include <nmod_poly.h>
#include <string.h>
#include <time.h> 
#include <sys/dir.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

/***************************************/
/********** MISCELLANEOUS   ************/
/***************************************/

/*Conversions*/
void 		  lu_to_fq(fq_nmod_ctx_t ctx, fq_nmod_t * field_element, unsigned long a, int j, unsigned long p);
void 		  lu_to_fmpz(fq_nmod_ctx_t ctx, fmpz_t * i_fmpz, unsigned long a, unsigned long p);
unsigned long fmpz_to_lu(fq_nmod_ctx_t ctx, fmpz_t * i_fmpz , int r, unsigned long p);

/*Shuffling*/
void shuffle(fq_nmod_ctx_t ctx, fq_nmod_t * coeffs_vec, unsigned long p, int r, int k, unsigned long block_size, int block_number);
void revSwap(fq_nmod_ctx_t ctx, fq_nmod_t * coeffs_vec, unsigned long p, int r, int k, unsigned long block_size, int block_number);

/*Printing*/
void print_vector_fq(fq_nmod_ctx_t ctx, fq_nmod_t * vector, unsigned long size);
void print_binary_tree(int * binary_tree, int depth);

/*Zero / Clear*/
void zero_subproducts_tree(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, unsigned long p, int depth);
void clear_subproducts_tree(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, unsigned long p, int depth);
void zero_p_evaluation_points(fq_nmod_ctx_t ctx, fq_nmod_t * p_evaluation_points, unsigned long p);
void zero_p_linear_modulus(fq_nmod_ctx_t ctx, fq_nmod_poly_t * p_linear_modulus, unsigned long p);
void zero_binary_tree(int * binary_tree, int depth, unsigned long p);

/***************************************/
/**********  p-point Eval   ************/
/***************************************/

//Build the perfect binary tree associated to p
void build_binary_tree(int * binary_tree, int depth, unsigned long p);

//Generate the linear modulus for the standard p-point evaluations
void generate_p_linear_modulus(fq_nmod_ctx_t ctx, fq_nmod_poly_t * p_linear_modulus, fq_nmod_t * p_evaluation_points, unsigned long p);

fq_nmod_poly_t * initialize_subproducts_tree(fq_nmod_ctx_t ctx, unsigned long p, int depth);

void 			 building_up_subproducts_tree(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, fq_nmod_poly_t * linear_moduli, 
											 					 int * binary_tree, int current_floor, int current_node, int depth, int * index);

//Evaluate the linearized polynomial L_{U_0^k} at the basis element v_{ind_basis}
fq_nmod_t * eval_L_U0k(fq_nmod_ctx_t ctx, int ind_basis, int k, fq_nmod_t * basis, unsigned long p);

//Compute the set Sxk of cardinality p for a given iteration k in {0, ... ,j-1} and for a given x in U_{k+1}
void Sxk(fq_nmod_ctx_t ctx, fq_nmod_t * points, fq_nmod_t* basis, fq_nmod_t elem_coset, unsigned long p, int i, int k, int r);

//Generate a univariate polynomial of F_{p^r}[X] of degree less than p by setting its coefficients 
void generate_polynomial_standard_multipoint_evaluation(fq_nmod_ctx_t ctx, fq_nmod_poly_t * f, fq_nmod_t * coefficients, unsigned long start, unsigned long p);

//For the original evaluation.
void generate_p_evaluation_points(fq_nmod_ctx_t ctx,fq_nmod_t * basis, fq_nmod_t * field_element, fq_nmod_t * p_evaluation_points, fq_nmod_poly_t * L, int k, unsigned long p);

void generate_linearized_polynomials(fq_nmod_ctx_t ctx, fq_nmod_poly_t * L, fq_nmod_t * basis, int i, int j, unsigned long p);

void standard_p_point_evaluation(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, int * binary_tree, fq_nmod_poly_t *linear_moduli,
									 				fq_nmod_poly_t poly, fq_nmod_t * coeffs_vec, int current_floor, int current_node, int depth, int * p_index, int * p_index_bis);

/*********************************************************/
/****  Multi-point evaluation of LCH polynomials  ********/
/*********************************************************/

//Evaluate the polynomial represented in the LCH-basis over F_p^r with higher memory locality of reference. This function corresponds to algorithm 1 of the article.
void multipoint_evaluation_LCH_polynomial_optimized(fq_nmod_ctx_t ctx, fq_nmod_t * basis, int *binary_tree, fq_nmod_t * coefficients, int j, unsigned long p);

//Evaluate the polynomial represented in the LCH-basis over F_p^r with lower memory locality of reference.
void multipoint_evaluation_LCH_polynomial_suboptimal(fq_nmod_ctx_t ctx, fq_nmod_t * basis, int * binary_tree, fq_nmod_t * coefficients, int j, unsigned long p);

//Evaluate the polynomial represented in the LCH-basis over over F_p^r following the article "Novel Polynomial Basis With Fast Fourier Transform and Its Application to Reed–Solomon Erasure Codes"
void original_evaluation_LCH_polynomial(fq_nmod_ctx_t ctx, fq_nmod_t* basis, fq_nmod_poly_t * L, fq_nmod_t * coefficients, int j, unsigned long p);

/**********************************************************************/
/********   Multi-point interpolation of LCH polynomials   ************/
/**********************************************************************/

void 			 p_point_interpolation(fq_nmod_ctx_t ctx, fq_nmod_t * basis, fq_nmod_poly_t * subproducts_tree, fq_nmod_poly_t * linear_moduli, int * binary_tree, 
																 fq_nmod_t * evaluations, unsigned long p, int k, int depth, int* p_index_evaluations);

fq_nmod_poly_t * linear_combination_for_linear_moduli(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, fq_nmod_poly_t * linear_moduli, int * binary_tree, 
													  fq_nmod_t * evaluations, fq_nmod_poly_t * r0, fq_nmod_poly_t *r1, int depth,  int current_floor, int current_node, int * p_index_bis);

void			 multipoint_interpolation_LCH_polynomial_optimized(fq_nmod_ctx_t ctx, fq_nmod_t * basis, int *binary_tree, fq_nmod_t * coefficients, int j, unsigned long p);

#endif