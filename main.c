#include <string.h>
#include "LCH.h"


int main (int argc, char ** argv){	

	if (argc != 3 ){
		printf("Usage : ./exe p r \n>p is a prime number, r is the degree of the extension field\n");
		return -1;
	}
	
	srand(time(NULL));
	fmpz_t 	 p;							//Characteristic of the field
	fmpz_init_set_ui(p,atoi(argv[1]));
	int r = atoi(argv[2]);				//Degree of the extension


	FILE * f;
	char file[128] = "";
	sprintf(file,"Multi-point_LCH_r_%d.dat", r);
	f = fopen (file,"a+");

	unsigned long seed1 = 1, 
				  seed2 = 1; //Fix the seeds			   			
	//seed1 = rand();
	//seed2 = rand();

	flint_rand_t   state;
	flint_randinit(state);
	flint_randseed(state, seed1, seed2); 

	//Sets the affine space U = U_i^j = <v_i, ... , v_{i+j-1}>.
	int i = 0,
		j = r,	 					
		tree_depth,
		*  perfect_binary_tree;	

	unsigned long p_ulong, nb_coeff;
	
	fq_nmod_ctx_t ctx;								
	fq_nmod_t 	  gen,
				  * coefficients, 					
				  basis[r]; 
							
	fq_nmod_ctx_init(ctx, p, r, "x");			
	fq_nmod_init(gen,ctx);
	fq_nmod_gen(gen, ctx);   						
				 					
	p_ulong 	= fmpz_get_ui(p);					
	tree_depth 	= floor(log2(p_ulong));				
	nb_coeff 	= pow(p_ulong,j);					//The LCH-polynomial has p^j coefficients

	
	//Generate basis of F_{p^r}
	fq_nmod_init(basis[0],ctx);
	fq_nmod_one(basis[0], ctx);
	for (int i = 1 ; i < r ; i ++){	
		fq_nmod_init(basis[i],ctx);
		fq_nmod_mul(basis[i], basis[i-1], gen, ctx);
	}
  	
  	//Generate the coefficients of the LCH-polynomial
	coefficients = malloc(sizeof(*coefficients) * nb_coeff);
	for (unsigned long a = 0; a < nb_coeff ; a ++){
		fq_nmod_init(coefficients[a], ctx);
		fq_nmod_randtest_not_zero(coefficients[a], state, ctx);
	}

	//Build the perfect binary tree associated to p for the standard p-point evaluations
	perfect_binary_tree = malloc(sizeof(*perfect_binary_tree) * ((1 << (tree_depth+1)) - 1) );
	build_binary_tree(perfect_binary_tree, tree_depth, p_ulong);
	
	//Print parameters.
	printf("\n\n> PARAMETERS : p = %lu | r = %d | U_i^j = <v_%d, ... , v_{%d}>\n", p_ulong, r, i, i+j-1);
	printf("***** BINARY TREE / Index : value *****\n"); print_binary_tree(perfect_binary_tree, tree_depth);
	//fq_nmod_print_pretty(basis[i],ctx); printf("   ");

	//Multi-point evaluations and interpolation.
	double t1,t2,dt;
	fprintf(f,"%lu		", fmpz_get_ui(p));

	t1 = clock();
	multipoint_evaluation_LCH_polynomial_optimized(ctx, basis, perfect_binary_tree, coefficients, j, p_ulong);			
	t2 = clock();
	dt = (double)((t2 - t1) / CLOCKS_PER_SEC);
	fprintf(f,"%f 	", dt);

	t1 = clock();
	multipoint_evaluation_LCH_polynomial_suboptimal(ctx, basis, perfect_binary_tree, coefficients,j, p_ulong);
	t2 = clock();
	dt = (double)((t2 - t1) / CLOCKS_PER_SEC);
	fprintf(f,"%f 	", dt);
			
	//Generate the linearized polynomials for the original evaluation
	fq_nmod_poly_t L[j];
	for (int a = 0 ; a < j ; a ++)
		fq_nmod_poly_init(L[a],ctx);
	generate_linearized_polynomials(ctx, L, basis, i,j, p_ulong);	

	t1 = clock();
	original_evaluation_LCH_polynomial(ctx, basis, L, coefficients, j, p_ulong);			
	t2 = clock();
	dt = (double)((t2 - t1) / CLOCKS_PER_SEC);
	fprintf(f,"%f 	\n", dt);
					
	multipoint_interpolation_LCH_polynomial_optimized(ctx, basis, perfect_binary_tree,coefficients, j, p_ulong);
			
	//********** Clear *********** 
	free(perfect_binary_tree);
	for (int a = 0 ; a < nb_coeff ; a ++)
		fq_nmod_clear(coefficients[a],ctx);
	free(coefficients);

	for (int i = 0 ; i < r ; i ++){	
		fq_nmod_poly_clear(L[i],ctx);
		fq_nmod_clear(basis[i],ctx);
	}

	fclose(f);
	flint_randclear(state);
	return EXIT_SUCCESS;		
}