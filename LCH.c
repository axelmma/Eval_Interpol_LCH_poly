#include "LCH.h"


/***************************************/
/**********  Conversions    ************/
/***************************************/

/*Convert an unsigned long into an fq_nmod_t*/
void lu_to_fq(fq_nmod_ctx_t ctx, fq_nmod_t * field_element, unsigned long a, int j, unsigned long p){
	int ind = 0, qi, ri;
	fq_nmod_t tmp;
	fq_nmod_init(tmp, ctx);

	if (a != 0){
		qi = a;
		ind = 0;			
		while (qi > 0){
			ri = qi%p;
			qi = qi/p;
			
			nmod_poly_set_coeff_ui(*field_element, ind++, ri);
		}
	}
}

/*Convert an unsigned long into an fmpz_t*/
void lu_to_fmpz(fq_nmod_ctx_t ctx, fmpz_t * i_fmpz, unsigned long a, unsigned long p){
	int ind = 0, qi, ri;

	if (a != 0){
		qi = a;
		ind = 0;			
		while (qi > 0){
			ri = qi%p;
			qi = qi/p;
			
			fmpz_set_ui(i_fmpz[ind], ri);
			ind ++;
		}
	}
}

/*Convert a fmpz_t into an unsigned long*/
unsigned long fmpz_to_lu(fq_nmod_ctx_t ctx, fmpz_t * i_fmpz , int r, unsigned long p){
	//int r = fq_nmod_ctx_degree(ctx);
	unsigned long tmp = 1;
	unsigned long i_lu = 0;

	//On multiplie les coefficients par la puissances de p correspondante pour obtenir la valeur de l'indice.
	for (int i = 0; i < r ; i ++){
		i_lu += (fmpz_get_ui(i_fmpz[i]) * tmp );
		tmp *=p;
	}	
	return i_lu;
}

/***************************************/
/**********  Shufflings	    ************/
/***************************************/

/*In-place perfect shuffle
The algorithm an be found here : https://arxiv.org/pdf/1204.1958.pdf (page 5) by Qiungxun Yang, John Ellis, Khalegh Mamakani and Frank Ruskey.*/
void shuffle(fq_nmod_ctx_t ctx, fq_nmod_t * coeffs_vec, unsigned long p, int r, int k, unsigned long block_size, int block_number){
	revSwap(ctx, coeffs_vec, p, r, k-1,block_size, block_number);
	revSwap(ctx, coeffs_vec, p, r, k, block_size, block_number);
}

void revSwap(fq_nmod_ctx_t ctx, fq_nmod_t * coeffs_vec, unsigned long p, int r, int k, unsigned long block_size, int block_number){
	unsigned long j = 0;
	fmpz_t * i_fmpz = malloc(sizeof(*i_fmpz)*r);

	for (int i = 0 ; i < r ; i++)
		fmpz_init(i_fmpz[i]);
		
	for (int i = 1 ; i < block_size ; i ++){
		lu_to_fmpz(ctx, i_fmpz, i, p);

		//Reverse i 
		for (int j = 0 ; j < k/2 ; j ++)
			fmpz_swap(i_fmpz[j], i_fmpz[(k-1)-j]);		

		j = fmpz_to_lu(ctx, i_fmpz, r, p);
		if (i < j)
			fq_nmod_swap(coeffs_vec[(block_number * block_size) + i],coeffs_vec[(block_number * block_size) + j],ctx);
		
		for (int a = 0 ; a < r ; a++)
			fmpz_zero(i_fmpz[a]);
	}
	//free(i_fmpz);
}	



/***************************************/
/********** 	Printing    ************/
/***************************************/

void print_vector_fq(fq_nmod_ctx_t ctx, fq_nmod_t * vector, unsigned long size){
	for (unsigned long i = 0 ; i < size ; i ++){
		fq_nmod_print_pretty(vector[i], ctx); printf("  ");
	}
	printf("\n\n");
}

void print_binary_tree(int * binary_tree, int depth){	
	for (int i = 0; i < (1 << (depth+1))- 1 ; i ++){
		printf("(%d) = %d\n",i, binary_tree[i]);	
	}
}

/***************************************/
/********  Zero and clear   ************/
/***************************************/

void zero_subproducts_tree(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, unsigned long p, int depth){
	int nb_nodes = (1 << (depth + 1)) + p - 1;							
 	for (int i = 0 ; i < nb_nodes ; i ++)
 		fq_nmod_poly_zero(subproducts_tree[i], ctx);		
}

void clear_subproducts_tree(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, unsigned long p, int depth){
	int nb_nodes = (1 << (depth + 1)) + p - 1;							
 	for (int i = 0 ; i < nb_nodes ; i ++)
 		fq_nmod_poly_clear(subproducts_tree[i], ctx);		
}

void zero_p_evaluation_points(fq_nmod_ctx_t ctx, fq_nmod_t * p_evaluation_points, unsigned long p){
	for (int a = 0 ; a < p ; a ++){
		fq_nmod_zero(p_evaluation_points[a], ctx);
	}
}

void zero_p_linear_modulus(fq_nmod_ctx_t ctx, fq_nmod_poly_t * p_linear_modulus, unsigned long p){
	for (int a = 0 ; a < p ; a ++){
		fq_nmod_poly_zero(p_linear_modulus[a], ctx);
	}
}

void zero_binary_tree(int * binary_tree, int depth, unsigned long p){

	for (int i = 0 ; i < ((1 << depth) - 1) ; i ++){
		binary_tree[(i<<1)+1] = 0;
		binary_tree[(i<<1)+2] = 0;
	}
}


/******************************************************************************
					p-point evaluations of standard polynomials
 ******************************************************************************/

void build_binary_tree(int * binary_tree, int depth, unsigned long p){
	binary_tree[0] = p;
	for (int i = 0 ; i < ((1 << depth) - 1) ; i ++){
		binary_tree[(i<<1)+1] = binary_tree[i]>>1;
		binary_tree[(i<<1)+2] = (binary_tree[i]>>1) + (binary_tree[i] & 1);
	}
}

/*For each of the p evaluation points as input, this fonction computes its corresponding linear modulo (m_i = X - X_i)*/
void generate_p_linear_modulus(fq_nmod_ctx_t ctx, fq_nmod_poly_t * p_linear_modulus, fq_nmod_t * p_evaluation_points, unsigned long p){
	fq_nmod_poly_t X;
	fq_nmod_poly_t tmp;

	fq_nmod_poly_init(tmp,ctx);
	fq_nmod_poly_init(X,ctx);
	fq_nmod_poly_gen(X,ctx);

	for (unsigned long i = 0 ; i < p ; i ++){
		fq_nmod_poly_set_coeff(tmp, 0, p_evaluation_points[i],ctx);
		fq_nmod_poly_sub(p_linear_modulus[i],X, tmp, ctx);
	}

	fq_nmod_poly_clear(X,ctx);
	fq_nmod_poly_clear(tmp,ctx);
}


/*Initialize the nods of the subproduct tree as polynomials*/
fq_nmod_poly_t * initialize_subproducts_tree(fq_nmod_ctx_t ctx, unsigned long p, int depth){
 	int nb_nodes = (1 << (depth + 1)) + p - 1;
 	fq_nmod_poly_t * subproducts_tree = malloc(sizeof(*subproducts_tree) * nb_nodes); 								
 	for (int i = 0 ; i < nb_nodes ; i ++)						
 			fq_nmod_poly_init(subproducts_tree[i],ctx);		 						  
 		
 	return subproducts_tree;
}

/*Recursively build up the subproducts tree for p-point evaluation/interpolation of standard polynomials*/
void building_up_subproducts_tree(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, fq_nmod_poly_t * linear_moduli, int * binary_tree, int current_floor, int current_node, int depth, int * index){
		if (current_floor == depth){
			if (binary_tree[current_node] == 1){
				fq_nmod_poly_set(subproducts_tree[current_node], linear_moduli[*index], ctx);
				(*index)++;
				return;
			}
			else{
				fq_nmod_poly_mul_KS(subproducts_tree[current_node], linear_moduli[*index] , linear_moduli[*index + 1], ctx);	
 				(*index) += 2;
 				return;
			} 
		}
		else{			
			building_up_subproducts_tree(ctx, subproducts_tree, linear_moduli, binary_tree, current_floor + 1, (current_node << 1) + 1, depth, index);
			building_up_subproducts_tree(ctx, subproducts_tree, linear_moduli, binary_tree, current_floor + 1, (current_node << 1) + 2, depth, index);

			fq_nmod_poly_mul_KS(subproducts_tree[current_node], subproducts_tree[(current_node << 1) + 1 ], subproducts_tree[(current_node << 1) + 2 ], ctx);
			return;
		}
}

/*Recursively compute the evaluation of the linearized polynomial L_{U_0^k} at v_{ind_basis} where v_{ind_basis} is an element of the basis of F_{p^r}
  - L_{U_0^k}(v_{ind_basis}) is computed from the polynomial L_{U_0^0} for which  L_{U_0^0}(v) = v for any v in F_{p^r}. Then L_{U_0^{k+1}} is obtained
  from the recursive definition of the linearized polynomials.
*/
fq_nmod_t * eval_L_U0k(fq_nmod_ctx_t ctx, int ind_basis, int k, fq_nmod_t * basis, unsigned long p){
	
	//For the case we actually just want to return L_U00(v) = v
	if (k == 0){
		fq_nmod_t * eval = malloc(sizeof(*eval));
		fq_nmod_init(*eval,ctx);
		fq_nmod_set(*eval, basis[ind_basis],ctx);
		return eval;
	}
	if (k == 1){
		fq_nmod_t tmp, * C;
		C = malloc(sizeof(*C));
		
		fq_nmod_init(tmp,ctx);
		fq_nmod_init(*C,ctx);
		
		fq_nmod_pow_ui(*C, basis[ind_basis], p, ctx);
		fq_nmod_pow_ui(tmp, basis[0], p-1, ctx);
		
		fq_nmod_mul(tmp,tmp,basis[ind_basis],ctx);

		fq_nmod_sub(*C,*C,tmp,ctx);
		//fq_nmod_clear(tmp,ctx);
		return C;
	}


	fq_nmod_t * A, * B, *C;
	A = malloc(sizeof(*A));
	B = malloc(sizeof(*B));
	C = malloc(sizeof(*C));
	
	fq_nmod_init(*A,ctx);
	fq_nmod_init(*B,ctx);
	fq_nmod_init(*C,ctx);
	
	A = eval_L_U0k(ctx, ind_basis,k-1,basis,p);
	B = eval_L_U0k(ctx, k-1,k-1,basis,p);

	fq_nmod_pow_ui(*B, *B, p-1, ctx);
	fq_nmod_mul(*B,*B, *A , ctx);
	fq_nmod_pow_ui(*A, *A, p, ctx);
	fq_nmod_sub(*C,*A,*B,ctx);
	
	fq_nmod_clear(*B,ctx);
	fq_nmod_clear(*A,ctx);
	free(B);
	free(A);
	return C;
}


/*Compute the set of p points S_{x,k} for the p-point evaluation/interpolation of standard polynomials.
  This follows formula (24) of the article.
*/
void Sxk(fq_nmod_ctx_t ctx, fq_nmod_t * points, fq_nmod_t* basis, fq_nmod_t elem_coset, unsigned long p, int i, int k, int r){

	fq_nmod_t *A,*B;
	fq_nmod_t sum,tmp;
	fq_nmod_init(sum,ctx);
	fq_nmod_init(tmp,ctx);
	A = malloc(sizeof(*A));
	B = malloc(sizeof(*B));
	
	fq_nmod_init(*A,ctx);
	fq_nmod_init(*B,ctx);
	unsigned long coeff;


	for (int c = 0 ; c < p ; c ++){

		if ( k == 0 ){
			fq_nmod_set(*A, basis[0],ctx);
		}
		else{
			A = eval_L_U0k(ctx, i+k, k, basis, p);
		}
		
		//Not executed if k == r-1
		for(int j = i+k+1 ; j < r; j ++ ){						//without loss of generality we set mu = 0 so j = i+k+1
			
			coeff = nmod_poly_get_coeff_ui(elem_coset, j);

			if ( k == 0 )
				fq_nmod_set(*B, basis[j],ctx);
			else
				B = eval_L_U0k(ctx, j, k, basis, p);
			
			fq_nmod_mul_ui(tmp,*B,coeff,ctx);
			fq_nmod_add(sum,sum,tmp,ctx); 
		}
		if (k != r-1){
			fq_nmod_mul_ui(tmp,*A,c,ctx);
			fq_nmod_add(points[c],sum, tmp,ctx);
			fq_nmod_zero(sum,ctx);
		}
		else
			fq_nmod_mul_ui(points[c],*A,c,ctx);	
	}

	fq_nmod_clear(*A,ctx);
	fq_nmod_clear(*B,ctx);
	free(A);
	free(B);
} 

/*Generate the polynomial of degree less than p for the p-point evaluation/interpolation by setting its p first coefficients*/
void generate_standard_polynomial(fq_nmod_ctx_t ctx, fq_nmod_poly_t * f, fq_nmod_t * coefficients, unsigned long start, unsigned long p){
	for (int d = 0; d < p ; d++){
		fq_nmod_poly_set_coeff(*f, d, coefficients[start + d],ctx);
	}
}

/*For a given x \in U_{k+1} compute the set {L_k(x) + c \cdot L_k{v_{i+k}} | c \in Fp}
-> (For original evaluation of LCH polynomials) */
void generate_p_evaluation_points(fq_nmod_ctx_t ctx,fq_nmod_t * basis, fq_nmod_t * field_element, fq_nmod_t * p_evaluation_points, fq_nmod_poly_t * L, int k, unsigned long p){
	fq_nmod_t eval, tmp;
	fq_nmod_init(eval,ctx);
	fq_nmod_init(tmp,ctx);

	fq_nmod_poly_evaluate_fq_nmod(eval, L[k], basis[k], ctx);
	fq_nmod_poly_evaluate_fq_nmod(tmp,  L[k], *field_element,ctx);
	
	for (int c = 0 ; c < p ; c++){
		fq_nmod_mul_ui(p_evaluation_points[c], eval, c, ctx);
		fq_nmod_add(p_evaluation_points[c], p_evaluation_points[c], tmp, ctx);
	}		
	fq_nmod_clear(eval,ctx);
	fq_nmod_clear(tmp, ctx);
}

/*
  (For original evaluation of LCH polynomials). 
  Used to generate the p points for standard p-point evaluations
  The linearized polynomials are defined according to the recursive definition given in relation (1) of the article.
*/
void generate_linearized_polynomials(fq_nmod_ctx_t ctx, fq_nmod_poly_t * L, fq_nmod_t * basis, int i, int j, unsigned long p){
	fq_nmod_t tmp;		
	fq_nmod_poly_t poly_tmp;
	fq_nmod_poly_init(poly_tmp,ctx);
	fq_nmod_init(tmp,ctx);

	fq_nmod_poly_gen(L[0],ctx);														
	for (int k = 0; k < j-1 ; k ++){
		fq_nmod_poly_pow(L[k+1],L[k], p, ctx);										
		fq_nmod_poly_evaluate_fq_nmod(tmp,L[k],basis[i+k], ctx);	    			
		fq_nmod_pow_ui(tmp,tmp, p-1, ctx );											
		fq_nmod_poly_scalar_mul_fq_nmod(poly_tmp, L[k], tmp, ctx);  	      		
		fq_nmod_poly_sub(L[k+1], L[k+1], poly_tmp, ctx);			     		
	}
	 fq_nmod_clear(tmp,ctx);	 
	 fq_nmod_poly_clear(poly_tmp,ctx);
}

//Compute the p-point evaluation of standard polynomials as in Ch.10 of Modern Computer Algebra by Gathen and Gerhard.
void standard_p_point_evaluation(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, int * binary_tree, fq_nmod_poly_t *linear_moduli, fq_nmod_poly_t poly, fq_nmod_t * coeffs_vec, int current_floor, int current_node, int depth, int * p_index, int * p_index_bis){
	fq_nmod_poly_t r0;
	fq_nmod_poly_t r1;
	fq_nmod_poly_init(r0,ctx);
	fq_nmod_poly_init(r1,ctx);

	if (current_floor == depth)
	{
		if (binary_tree[current_node] == 1){
			fq_nmod_poly_get_coeff(coeffs_vec[*p_index], poly, 0, ctx);
			fq_nmod_poly_clear(poly, ctx);
			(*p_index)++;
			(*p_index_bis)++;
			return;
		}
		else{
			fq_nmod_poly_rem(r0, poly, linear_moduli[*p_index_bis], ctx);   
			fq_nmod_poly_get_coeff(coeffs_vec[(*p_index)], r0, 0, ctx);
			fq_nmod_poly_rem(r1, poly, linear_moduli[(*p_index_bis)+1],ctx);
			fq_nmod_poly_get_coeff(coeffs_vec[(*p_index)+1], r1, 0, ctx);
			fq_nmod_poly_clear(r0, ctx);
			fq_nmod_poly_clear(r1, ctx);
			(*p_index)+=2;
			(*p_index_bis)+=2;
			return;	

		}
	}
	fq_nmod_poly_rem(r0, poly, subproducts_tree[(current_node<<1)+1], ctx);	
	fq_nmod_poly_rem(r1, poly, subproducts_tree[(current_node<<1)+2], ctx);
			 
	standard_p_point_evaluation(ctx, subproducts_tree, binary_tree, linear_moduli, r0, coeffs_vec, current_floor + 1, (current_node<<1)+1, depth, p_index,p_index_bis);
	standard_p_point_evaluation(ctx, subproducts_tree, binary_tree, linear_moduli, r1, coeffs_vec, current_floor + 1, (current_node<<1)+2, depth, p_index,p_index_bis);
	return;
}

/******************************************************************************
	Multi-point evaluations of polynomials represented in the LCH-basis
 ******************************************************************************/

/*This function executes algorithm 1 of the article for the fast multipoint evaluation of a polynomial represented in the LCH-basis over F_p^r*/
void multipoint_evaluation_LCH_polynomial_optimized(fq_nmod_ctx_t ctx, fq_nmod_t * basis, int *binary_tree, fq_nmod_t * coefficients, int j, unsigned long p){
	int index, index_bis,
		* p_index = &index,
		* p_index_bis = &index_bis,
		depth = floor(log2(p));
	
	unsigned long block_size,				//For perfect shuffles
				  block_number = 0,			//
				  index_bloc,				//
				  nb_block,
				  step,
				  pj = pow(p,j),
				  nb_eval_points = pj,
				  block_size_vand = pj;
	
	fq_nmod_poly_t * p_linear_moduli = malloc(sizeof(*p_linear_moduli)*p),			  //The linear polynomials for building the subproduct tree
				   * poly_monomial_basis = malloc(sizeof(*poly_monomial_basis)),	  //The standard polynomial evaluated at p points
				   * subproducts_tree = initialize_subproducts_tree(ctx,p, depth); 	  //The subproduct tree for the p-point evaluations
	
	fq_nmod_t 	   * p_evaluation_points = malloc(sizeof(*p_evaluation_points)*p), 	
			  	   * coset_element = malloc(sizeof(*coset_element));

	fq_nmod_init(*coset_element, ctx);
	fq_nmod_poly_init(*poly_monomial_basis, ctx);

	for (int a = 0 ; a < p ; a ++){
		fq_nmod_init(p_evaluation_points[a], ctx);
		fq_nmod_poly_init(p_linear_moduli[a], ctx);
	}

	
	printf("		*** MULTIPOINT EVALUATION WITH HIGHER (!) MEMORY LOCALITY OF REFERENCE of polynomials represented in the LCH-basis ***\n\n");
	printf("***** INPUT POLYNOMIAL COEFFICIENTS *****\n"); print_vector_fq(ctx, coefficients, pj);
	//Step 2.
	for (int k = j-1 ; k >= 0 ; k --){
		block_size_vand  /= p;															//There are p^k identical Vandermonde matrices within each block
		step = block_size_vand * p;														//Next element of U_{k+1} in lexicographic order (p^{k+1})
		index_bloc = 0;																	//Index of blocks of Vandermonde matrices.
		index = 0;																		//For standard p-point evaluation. 
		index_bis = 0;
	
		//Step 3.
		if(k == j-1){																	//B_{j-1} = P_{1,j}.	
			shuffle(ctx, coefficients,p,j,j, pj, 0);									//Only one block one size p^j
		}
		else{																			//For any k = j-2 ... 0												
			block_number = 0;			
			block_size = pow(p,(k+2));													
			nb_block = pj/block_size;

			//Shuffling of blocks of size p^{k+2}
			for (unsigned long a = 0 ; a < nb_block ; a++){ 							//Shuffling of the same block k+1 times to the left to perform a right shift.
				for (int b = 0 ; b < k+1 ; b ++)				
					shuffle(ctx, coefficients,p, k+2, k+2, block_size, block_number);
				block_number ++;
			}

			block_number = 0;
			block_size /= p;		
			nb_block *= p; 		
			for (unsigned long a = 0 ; a < nb_block; a ++)								//Shuffling of blocks of size p^{k+1}
					shuffle(ctx, coefficients,p,k+1,k+1, block_size, block_number++);
		}																				
		
		//Step 4. Loop on all elements of U_{k+1} with respect to the total order defined over F_{p^r}
		for (unsigned long x = 0 ; x < nb_eval_points ; x += step){		
			*p_index = 0;
			lu_to_fq(ctx, coset_element, x, j, p);
			//Step 5.		
			Sxk(ctx, p_evaluation_points, basis, *coset_element, p, 0, k, j);			//Generate the p points for the standard p-point evaluations as in relation (24)
			
			generate_p_linear_modulus(ctx, p_linear_moduli, p_evaluation_points, p);	//Generate the p linear modulus to build up the subproducts tree for the standard p-point evaluations
			building_up_subproducts_tree(ctx, subproducts_tree, p_linear_moduli, binary_tree, 0, 0, depth, p_index);		

			//Step 6. Loop on all index of Vandermonde matrices for each block
			for (unsigned long l = 0; l < block_size_vand ; l ++){
				*p_index = 0;															//Reset for next iteration. Pointer on the indexes of the linear moduli for p-point eval
				//Step 7.
				generate_standard_polynomial(ctx, poly_monomial_basis, coefficients, (index_bloc * (block_size_vand * p)) + (l * p) , p);	//Generate the standard polynomial as in relation (23) from the coefficients of the LCH polynomial 
				//Steps 8 and 9.
				standard_p_point_evaluation(ctx, subproducts_tree, binary_tree, p_linear_moduli, *poly_monomial_basis, coefficients, 0, 0, depth, p_index_bis, p_index);
				fq_nmod_poly_zero(*poly_monomial_basis,ctx);
			}
			index_bloc ++;	//Increment for the next sequence of Vandermonde matrices (for the next x)

			zero_p_evaluation_points(ctx, p_evaluation_points, p);
			zero_p_linear_modulus(ctx, p_linear_moduli, p);
			zero_subproducts_tree(ctx, subproducts_tree,p, depth);
			fq_nmod_zero(*coset_element,ctx);
		}
	}

	printf("***** OUTPUT POLYNOMIAL EVALUATIONS *****\n");
	print_vector_fq(ctx, coefficients, nb_eval_points);
	
	//Clear
	clear_subproducts_tree(ctx, subproducts_tree, p, depth);

	for (int a = 0 ; a < p ; a ++){
		fq_nmod_clear(p_evaluation_points[a], ctx);
		fq_nmod_poly_clear(p_linear_moduli[a], ctx);
	}
	fq_nmod_clear(*coset_element, ctx);
	fq_nmod_poly_clear(*poly_monomial_basis, ctx);
}



void multipoint_evaluation_LCH_polynomial_suboptimal(fq_nmod_ctx_t ctx, fq_nmod_t * basis, int * binary_tree, fq_nmod_t * coefficients, int j, unsigned long p){
	int index, index_bis,
		* p_index = &index,
		* p_index_bis = &index_bis,
		depth = floor(log2(p));

	unsigned long card_coset = 1,  
				  nb_blocks_vand = 0,
				  index_vand,
				  step,
				  nb_eval_points = pow(p, j);
	
	fq_nmod_poly_t * subproducts_tree = initialize_subproducts_tree(ctx,p, depth),
				   * p_linear_moduli = malloc(sizeof(*p_linear_moduli)*p),
				   * poly_monomial_basis  = malloc(sizeof(*poly_monomial_basis)); 
	
	fq_nmod_t * p_evaluation_points = malloc(sizeof(*p_evaluation_points)*p), 
			  * coset_element = malloc(sizeof(*coset_element));

	fq_nmod_init(*coset_element,ctx);
	fq_nmod_poly_init(*poly_monomial_basis,ctx);

	for (int a = 0 ; a < p ; a ++){
		fq_nmod_init(p_evaluation_points[a], ctx);
		fq_nmod_poly_init(p_linear_moduli[a], ctx);
	}

	/*Evaluation*/
	printf("\n\n		*** MULTIPOINT EVALUATION WITH LOWER (!) MEMORY LOCALITY OF REFERENCE of polynomials represented in the LCH-basis ***\n\n");
	printf("***** INPUT POLYNOMIAL COEFFICIENTS *****\n"); print_vector_fq(ctx, coefficients, nb_eval_points);
	for (int k = j-1 ; k >= 0 ; k --){
		step 			= pow(p, k+1);					//Next elem of U_{k+1} in lexicographic order
		nb_blocks_vand  = step/p;	    				//Number of blocks of Vandermonde matrices : There are p^k at iteration k
		if (k != j-1)  
			card_coset 	*= p ; 							//Number of distinct Vandermonde matrices within each block : there are #U_{k+1} at iteration k
		index_vand 		= 0;							//Index of Vandermonde matrices within a block
		index 			= 0;
		index_bis		= 0;							//For standard p point evaluation
		
	
		for (int i = 0 ; i < j-k ; i ++)
			shuffle(ctx, coefficients, p, j, j, nb_eval_points, 0);


		for (unsigned long l = 0; l < nb_blocks_vand ; l ++){
			for (unsigned long x = 0 ; x < nb_eval_points ; x += step){
				index = 0;
				lu_to_fq(ctx, coset_element, x, j, p);		
				Sxk(ctx, p_evaluation_points, basis, *coset_element, p, 0, k, j);
				
				generate_p_linear_modulus(ctx, p_linear_moduli, p_evaluation_points, p);
				building_up_subproducts_tree(ctx, subproducts_tree, p_linear_moduli, binary_tree, 0, 0, depth, p_index);

				generate_standard_polynomial(ctx, poly_monomial_basis, coefficients, (l * card_coset * p) + (index_vand * p) , p);

				index = 0;
				standard_p_point_evaluation(ctx, subproducts_tree, binary_tree, p_linear_moduli, *poly_monomial_basis, coefficients, 0, 0, depth, p_index_bis, p_index);

				zero_p_evaluation_points(ctx, p_evaluation_points, p);
				zero_p_linear_modulus(ctx, p_linear_moduli, p);
				zero_subproducts_tree(ctx, subproducts_tree,p, depth);
				fq_nmod_zero(*coset_element,ctx);
				fq_nmod_poly_zero(*poly_monomial_basis,ctx);
				index_vand ++;
			}
			index_vand = 0;	
		}

		for (int a = 0; a < k ; a ++)
			shuffle(ctx, coefficients, p, j, j, nb_eval_points, 0);
	}
	
	printf("***** OUTPUT POLYNOMIAL EVALUATIONS *****\n");
	print_vector_fq(ctx, coefficients, nb_eval_points);

	/*Clear*/
	clear_subproducts_tree(ctx, subproducts_tree, p, depth);
	for (int a = 0 ; a < p ; a ++){
		fq_nmod_clear(p_evaluation_points[a], ctx);
		fq_nmod_poly_clear(p_linear_moduli[a], ctx);
	}
	fq_nmod_clear(*coset_element,ctx);
	fq_nmod_poly_clear(*poly_monomial_basis,ctx);
}


void original_evaluation_LCH_polynomial(fq_nmod_ctx_t ctx, fq_nmod_t* basis, fq_nmod_poly_t * L, fq_nmod_t * coefficients, int j, unsigned long p){
	unsigned long card_coset = 0,  
				  nb_blocks = 0,
				  index_sub_block,
				  step, 
				  nb_eval_points = pow(p, j); 

	fq_nmod_t * p_evaluation_points = malloc(sizeof(*p_evaluation_points)*p), 
			  * coset_element 		= malloc(sizeof(*coset_element));

	fq_nmod_poly_t  * poly_monomial_basis  = malloc(sizeof(*poly_monomial_basis)); 

	fq_nmod_poly_init(*poly_monomial_basis, ctx);
	fq_nmod_init(*coset_element,ctx);


	for (int a = 0 ; a < p ; a ++)
		fq_nmod_init(p_evaluation_points[a], ctx);
	
	printf("\n\n 				*** ORIGINAL EVALUATION of polynomials represented in the LCH-basis ***\n\n");
	printf("***** INPUT POLYNOMIAL COEFFICIENTS *****\n"); print_vector_fq(ctx, coefficients, nb_eval_points);
	for (int k = j-1 ; k >= 0 ; k --){
		nb_blocks  = pow(p, k);	
		step = nb_blocks * p;    		
		card_coset = pow(p, (j-k-1)); 		
	 	index_sub_block = 0;

	 	for (int i = 0 ; i < j-k ; i ++)
			shuffle(ctx, coefficients, p, j, j, nb_eval_points, 0);
			
		for (unsigned long l = 0; l < nb_blocks ; l ++){
			for (unsigned long x = 0 ; x < nb_eval_points ; x += step){
				lu_to_fq(ctx, coset_element, x, j, p);		
				generate_p_evaluation_points(ctx, basis, coset_element, p_evaluation_points, L, k, p);
				generate_standard_polynomial(ctx, poly_monomial_basis, coefficients, (l * card_coset * p) + (index_sub_block * p) , p);
				
				//Horner 
				for (int d = 0 ; d < p ; d ++)
					fq_nmod_poly_evaluate_fq_nmod(coefficients[(l * card_coset * p) + (index_sub_block * p) + d], *poly_monomial_basis, p_evaluation_points[d], ctx);
				
				zero_p_evaluation_points( ctx, p_evaluation_points, p);			
				fq_nmod_poly_zero(*poly_monomial_basis, ctx);
				fq_nmod_zero(*coset_element, ctx);
				index_sub_block ++;
			}
			index_sub_block = 0;	
		}

		for (int a = 0; a < k ; a ++)
		shuffle(ctx, coefficients, p, j, j, nb_eval_points, 0);

	}
	printf("***** OUTPUT POLYNOMIAL EVALUATIONS *****\n");
	print_vector_fq(ctx, coefficients, nb_eval_points);

	for (int a = 0 ; a < p ; a ++){
		fq_nmod_clear(p_evaluation_points[a], ctx);
	}
}


/************************************************************************************************/
/********   Multi-point interpolation of polynomials represented in the LCH-basis.   ************/
/************************************************************************************************/

void p_point_interpolation(fq_nmod_ctx_t ctx, fq_nmod_t * basis, fq_nmod_poly_t * subproducts_tree, fq_nmod_poly_t * linear_moduli, int * binary_tree, fq_nmod_t * evaluations, unsigned long p, int k, int depth, int* p_index_evaluations){

	fq_nmod_t * evals_m_prime = malloc(sizeof(*evals_m_prime));
	fq_nmod_init(*evals_m_prime, ctx);

	evals_m_prime = eval_L_U0k(ctx, k, k, basis, p);					//i = 0 || m'(x) = -L_i(v_i)^{p-1}

	fq_nmod_pow_ui(*evals_m_prime, *evals_m_prime, p-1,ctx);						
	fq_nmod_neg(*evals_m_prime,*evals_m_prime,ctx);
	fq_nmod_inv(*evals_m_prime,*evals_m_prime,ctx);								
	
	//c_i = v_i. Constant. We multiply it after. 
	fq_nmod_poly_t * poly = malloc(sizeof(*poly));
	fq_nmod_poly_init(*poly,ctx);
	
	int index = 0;
	int * p_index = &index;

	poly = linear_combination_for_linear_moduli(ctx, subproducts_tree, linear_moduli, binary_tree, evaluations, NULL, NULL, depth, 0,0, p_index);
	fq_nmod_poly_scalar_mul_fq_nmod(*poly, *poly, *evals_m_prime, ctx);		//s_i.P(x)
	
	//Update the vecteur of evaluations for next iteration
	for (int i = 0; i < p ; i++){
		fq_nmod_zero(evaluations[i],ctx);
		fq_nmod_poly_get_coeff(evaluations[i], *poly,i,ctx); ;
	}

	fq_nmod_clear(*evals_m_prime,ctx);	
	fq_nmod_poly_clear(*poly,ctx);

}

/**/
fq_nmod_poly_t * linear_combination_for_linear_moduli(fq_nmod_ctx_t ctx, fq_nmod_poly_t * subproducts_tree, fq_nmod_poly_t * linear_moduli,
														int * binary_tree, fq_nmod_t * evaluations, fq_nmod_poly_t * r0, fq_nmod_poly_t *r1, int depth,  int current_floor, int current_node, int * p_index){

	if (current_floor == depth){
		fq_nmod_poly_t *poly = malloc(sizeof(*poly));
		fq_nmod_poly_init(*poly,ctx);

		if (binary_tree[current_node] == 1){
			fq_nmod_poly_set_fq_nmod(*poly, evaluations[*p_index], ctx);
			(*p_index)++;
			return poly;
		}
		else if (binary_tree[current_node] == 2){
			fq_nmod_poly_t tmp;
			fq_nmod_poly_init(tmp, ctx);

			fq_nmod_poly_scalar_mul_fq_nmod(tmp, linear_moduli[(*p_index)+1], evaluations[(*p_index)],ctx);	
			fq_nmod_poly_scalar_mul_fq_nmod(*poly, linear_moduli[(*p_index)], evaluations[(*p_index)+1],ctx);

			fq_nmod_poly_add(*poly, *poly, tmp, ctx);
			(*p_index)+=2;
						return poly;
		}
	}

	fq_nmod_poly_t *poly = malloc(sizeof(*poly));
	fq_nmod_poly_init(*poly, ctx);
	fq_nmod_poly_t * A = malloc(sizeof(*A));
	fq_nmod_poly_t * B = malloc(sizeof(*B));
	fq_nmod_poly_init(*A, ctx);
	fq_nmod_poly_init(*B, ctx);

	fq_nmod_poly_set(*A, *linear_combination_for_linear_moduli(ctx, subproducts_tree, linear_moduli, binary_tree, evaluations, A, B,  depth, current_floor+1, (current_node<<1)+1, p_index ), ctx);
	fq_nmod_poly_set(*B, *linear_combination_for_linear_moduli(ctx, subproducts_tree, linear_moduli, binary_tree, evaluations, A, B,  depth, current_floor+1, (current_node<<1)+2, p_index ), ctx);	


	fq_nmod_poly_t tmp;
	fq_nmod_poly_init(tmp, ctx);
	fq_nmod_poly_mul_KS(tmp, subproducts_tree[(current_node<<1)+2], *A, ctx);
	fq_nmod_poly_mul_KS(*poly, subproducts_tree[(current_node<<1)+1], *B, ctx);
	fq_nmod_poly_add(*poly, *poly, tmp, ctx);
	return poly;
}



void multipoint_interpolation_LCH_polynomial_optimized(fq_nmod_ctx_t ctx, fq_nmod_t * basis, int *binary_tree, fq_nmod_t * coefficients, int j, unsigned long p){
	int index, index_bis,
		* p_index = &index,
		* p_index_bis = &index_bis,
		depth = floor(log2(p));
	
	unsigned long block_size,				//For perfect shuffles
				  block_number = 0,			//
				  index_bloc,				//
				  nb_block,
				  step,
				  pj = pow(p,j),
				  nb_eval_points = pj,
				  block_size_vand = 1;
	
	fq_nmod_poly_t * p_linear_moduli = malloc(sizeof(*p_linear_moduli)*p),			  //The linear polynomials for building the subproduct tree
				   * poly_monomial_basis = malloc(sizeof(*poly_monomial_basis)),	  //The standard polynomial evaluated at p points
				   * subproducts_tree = initialize_subproducts_tree(ctx,p, depth); 	  //The subproduct tree for the p-point evaluations
	
	fq_nmod_t 	   * p_evaluation_points = malloc(sizeof(*p_evaluation_points)*p), 	
			  	   * coset_element = malloc(sizeof(*coset_element)),
			  	   * evals = malloc(sizeof(*evals)*p);
	
	fq_nmod_init(*coset_element, ctx);
	fq_nmod_poly_init(*poly_monomial_basis, ctx);

	for (int a = 0 ; a < p ; a ++){
		fq_nmod_init(p_evaluation_points[a], ctx);
		fq_nmod_poly_init(p_linear_moduli[a], ctx);
		fq_nmod_init(evals[a],ctx);
	}

	//Evaluation of the factorized product
	printf("\n\n 		*** MULTIPOINT INTERPOLATION WITH HIGHER MEMORY LOCALITY OF REFERENCE of polynomials represented in the LCH-basis ***\n\n");
	printf("***** INPUT POLYNOMIAL EVALUATIONS *****\n"); print_vector_fq(ctx, coefficients, nb_eval_points);
	for (int k = j-1 ; k >= 0 ; k --){
		if (k != j-1)
			block_size_vand  *= p;														//There are p^k identical Vandermonde matrices within each block
		
		step = block_size_vand * p;														//Next element of U_{k+1} in lexicographic order (p^{k+1})
		index_bloc = 0;																	//Index of blocks of Vandermonde matrices.
		index = 0;																		//For standard p-point evaluation. 
		index_bis = 0;
	
		//Loop on all elements of U_{k+1} in lexicographic order
		for (unsigned long x = 0 ; x < nb_eval_points ; x += step){		
			*p_index = 0;
			lu_to_fq(ctx, coset_element, x, j, p);		
			Sxk(ctx, p_evaluation_points, basis, *coset_element, p, 0, j-k-1, j);
		
			generate_p_linear_modulus(ctx, p_linear_moduli, p_evaluation_points, p);
			building_up_subproducts_tree(ctx, subproducts_tree, p_linear_moduli, binary_tree, 0, 0, depth, p_index);
			
			//Loop on all index of Vandermonde matrices for each block
			for (unsigned long l = 0; l < block_size_vand ; l ++){

				for (int i = 0 ; i < p ; i ++)
					fq_nmod_set(evals[i], coefficients[*p_index_bis + i], ctx);

				p_point_interpolation(ctx, basis, subproducts_tree, p_linear_moduli, binary_tree, evals, p,j-k-1, depth, p_index_bis);
			
			for (int i = 0 ; i < p ; i ++)
					fq_nmod_set(coefficients[*p_index_bis + i], evals[i], ctx);
			(*p_index_bis)+=p;

				zero_p_evaluation_points(ctx, evals, p);
				fq_nmod_poly_zero(*poly_monomial_basis,ctx);
			}
			index_bloc ++;

			zero_p_evaluation_points(ctx, p_evaluation_points, p);
			zero_p_linear_modulus(ctx, p_linear_moduli, p);
			zero_subproducts_tree(ctx, subproducts_tree,p, depth);
			fq_nmod_zero(*coset_element,ctx);
		}


		if(j-k-1 == j-1){																	//B_{j-1} = P_{1,j}.
			for (int b = 0 ; b < j-1 ; b ++)	
				shuffle(ctx, coefficients,p,j,j, pj, 0);									
		}
		else{			
			block_number = 0;
			block_size = pow(p,(j-k));
			nb_block = pj/block_size;

					
			for (unsigned long a = 0 ; a < nb_block; a ++){								
				for (int b = 0 ; b < j-k-1 ; b ++)
					shuffle(ctx, coefficients,p, j-k, j-k, block_size, block_number);
				block_number++;
			}
			////////////////////////////////																						
			block_number = 0;			
			block_size *= p;															
			nb_block /= p;

			
			for (unsigned long a = 0 ; a < nb_block ; a++){ 									
				shuffle(ctx, coefficients,p, (j-k)+1, (j-k)+1, block_size, block_number++);
			}
		}			
	}

	printf("***** OUTPUT POLYNOMIAL COEFFICIENTS *****\n");
	print_vector_fq(ctx, coefficients, nb_eval_points);
	
	//Clear
	clear_subproducts_tree(ctx, subproducts_tree, p, depth);

	for (int a = 0 ; a < p ; a ++){
		fq_nmod_clear(p_evaluation_points[a], ctx);
		fq_nmod_poly_clear(p_linear_moduli[a], ctx);
	}
	fq_nmod_clear(*coset_element, ctx);
	fq_nmod_poly_clear(*poly_monomial_basis, ctx);
}


