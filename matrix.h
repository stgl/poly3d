/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/****************************************************************************
* FILE: matrix.h
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* Header file for matrix.c
*****************************************************************************/


/***************************** Includes/Defines *****************************/
#define		TRUE				1			/* true flag					*/
#define		FALSE				0			/* false flag					*/
#define		INVERSE_ROT			TRUE		/* inverse coord rotation flag	*/
#define		FORWARD_ROT			FALSE		/* forward coord rotation flag	*/
#define		INVERSE_TRANS		TRUE		/* inverse coord trans flag		*/
#define		FORWARD_TRANS		FALSE		/* forward coord trans flag		*/


/****************************** Type Definitions ****************************/
typedef double		component_t;
typedef component_t	vector_t[3];
typedef component_t	matrix_t[3][3];


/************************* ANSI function declarations ***********************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void		add_matrices(matrix_t m1, matrix_t m2, matrix_t m3);
void		add_vectors(vector_t v1, vector_t v2, vector_t v3);
void		cauchy(matrix_t stress, vector_t normal_v, vector_t traction_v);
void		copy_vector(vector_t v, vector_t v_copy);
void		copy_matrix(matrix_t m, matrix_t m_copy);
void		cross_product(vector_t v1, vector_t v2, vector_t v3);
component_t	dot_product(vector_t v1, vector_t v2);
void		initialize_matrix(matrix_t m, component_t value);
void		initialize_vector(vector_t v, component_t value);
void		normalize_vector(vector_t v);
void		matrix_mult(matrix_t a, matrix_t b, matrix_t c);
void		matrix_vector_mult(matrix_t a, vector_t b, vector_t c);
void		principal(matrix_t stress, vector_t prin, matrix_t traj);
void		rotate_tensor(int inverse_rot, matrix_t a, matrix_t s);
void		rotate_vector(int inverse_rot, matrix_t a, vector_t x);
void		scalar_matrix_mult(component_t s, matrix_t m1, matrix_t m2);
void		scalar_vector_mult(component_t s, vector_t v1, vector_t v2);
void		subtract_matrices(matrix_t m1, matrix_t m2, matrix_t m3);
void		subtract_vectors(vector_t v1, vector_t v2, vector_t v3);
void		transpose_matrix(matrix_t a);
void		transform_position_vector(int inverse_trans, vector_t origin,
			matrix_t a, vector_t x);
component_t	vector_magnitude(vector_t v);

float		distance(vector_t v, vector_t w);


/************************** K&R function declarations ***********************/
#else
void		add_matrices();
void		add_vectors();
void		cauchy();
void		copy_vector();
void		copy_matrix();
void		cross_product();
component_t	dot_product();
void		initialize_matrix();
void		initialize_vector();
void		normalize_vector();
void		matrix_mult();
void		matrix_vector_mult();
void		principal();
void		rotate_tensor();
void		rotate_vector();
void		scalar_matrix_mult();
void		scalar_vector_mult();
void		subtract_vectors();
void		subtract_matrices();
void		transpose_matrix();
void		transform_position_vector();
component_t	vector_magnitude();

float		distance();
#endif
