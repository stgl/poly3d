//#define DEBUG 
/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: poly3d.c, version 0.0 (beta)
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* A 3-D polygonal boundary element code based on superposition of angular
* dislocations (Yoffe, 1956; Comninou & Dunders, 1966).  Translated and
* expanded from a triangular element FORTRAN code written by M. Jeyakumaran
* at Northwestern University (Jeyakumaran et al., 1992).
*
* See the "Poly3D Users Manual" for more information.
*****************************************************************************/


/**************************** Revision History *******************************
*
* If you change this code in ANY way, describe and initial those changes
* below.  Put a comment line containing the word CHANGE (all uppercase)
* where the changes are made, so they will be easy to find in the source code.
*
* DATE		FILE			NAME					INTLS
* --------  --------------  ----------------------  -----
* Jun-93	poly3d.c		Andrew L. Thomas		ALT
*		Version 0.0 (beta) of Poly3D completed
* Nov-97        poly3d.c                Yann Lagalaye
*		Newpoly3d (execut) Correction of "shadow effect"
*****************************************************************************/


/******************************** Includes **********************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include "matrix.h"
#include "safetan.h"
#include "getoptPoly3D.h"
#include "getwords.h"
#include "infcoeff.h"
#include "elastic.h"
#include "pi.h"
#include "nrutil.h"
#include "nr.h"


/******************************** Constants *********************************/
/*-------------
 MISCELLANEOUS
--------------*/
#define MAXFILE	 			256			/* Max length of file names			*/
#define MAX_ERROR_MSG		256			/* Max length of error messages		*/
#define MAXLINE				256			/* Max line length for getwords()	*/
#define MAXWORDS			50			/* Max # words on getwords() line	*/
#define GLOBAL_NAME			"global"	/* Global coord system name			*/
#define ELT_CSYS_NAME		"elocal"	/* Element coord system name		*/
#define END_STMT			"end"		/* End flag for input file sections	*/
#define COMMENT_CHAR		'*'			/* Input file comment character		*/
#define CONTINUE_CHAR		'\\'		/* Input file line continuation char*/
#define ERROR				-1			/* Return value	for function errors	*/
#define FALSE				0			/* False flag						*/
#define TRUE				1			/* True flag						*/
#define BVECTOR_BC			0			/* Burgers vector BC component flag	*/
#define TRACTION_BC			1			/* Traction BC component flag		*/

/*-----------
 PROGRAM INFO
-------------*/
#define PROGRAM			"poly3d.c"		
#define VERSION			"Beta-Release"
#ifdef __DATE__						
#define COMPILE_DATE	__DATE__		
#else
#define COMPILE_DATE "Date Unavailable"
#endif

/*--------------------------------------
 NUMERICAL LIMITS USED WHEN CHECKING...
---------------------------------------*/
#define SWAP_TINY			1.0e-10		/* ...if vertices must be swapped	*/
#define TINY_ANGLE			0.5*PI/180.	/* ...if elt coord sys can be calc	*/
#define BVERT_TINY			1.0e-10		/* ...if point lies below a vertex	*/
#define COPLANAR_LIMIT		30.			/* ...if elt vertices are co-planar	*/

/*-------------------------------------
 PRINT OPTION ARRAY SIZE AND POSITIONS
--------------------------------------*/
#define NUM_PR_OPTS			5
#define DISPL				0
#define	STRAIN				1
#define PSTRAIN				2
#define STRESS				3
#define PSTRESS				4

/*--------------------------------------------------------------
 CHARS USED IN INPUT FILE PRINT STRINGS TO ENABLE PRINT OPTIONS
---------------------------------------------------------------*/
#define DISPL_CHAR			'd'
#define STRESS_CHAR			's'
#define STRAIN_CHAR			'e'
#define TRACTION_CHAR		't'
#define BVECTOR_CHAR		'b'
#define PRINCIPAL_CHAR		'p'

/*----------------------------------------
 INPUT FILE FORMAT FOR DEFINING CONSTANTS
-----------------------------------------*/
#define CONST_NAME_POS		0
#define CONST_VALUE_POS		2
#define CONST_NUM_PARAMS	3

/*------------------------------------------------------
 INPUT FILE FORMAT FOR DEFINING USER COORDINATE SYSTEMS
-------------------------------------------------------*/
#define CS_NAME_POS			0
#define CS_PARENT_POS		1
#define CS_ORIGIN_POS		2
#define CS_ROT_POS			5
#define CS_ROT_ORDER_POS	8
#define CS_NUM_PARAMS		9

/*------------------------------------------------
 INPUT FILE FORMAT FOR DEFINING OBSERVATION GRIDS
-------------------------------------------------*/
#define OG_NAME_POS			0
#define OG_DIMEN_POS		1
#define OG_PRINT_OPS_POS	2
#define OG_INPUT_CSYS_POS	3
#define OG_OBSPT_CSYS_POS	4
#define OG_DATA_CSYS_POS	5
#define OG_BEGIN_POS		6
#define OG_END_POS			9
#define OG_NUMPTS_POS		12
#define OG_MIN_NUM_PARAMS	9

/*---------------------------------------
 INPUT FILE FORMAT FOR DEFINING VERTICES
----------------------------------------*/
#define V_CHAR				'v'
#define V_CHAR_POS			0
#define V_NAME_POS			1
#define V_CSYS_POS			2
#define V_X_POS				3
#define V_NUM_PARAMS		6

/*--------------------------------------
 INPUT FILE FORMAT FOR DEFINING OBJECTS
---------------------------------------*/
#define OBJ_CHAR			'o'
#define OBJ_CHAR_POS		0
#define OBJ_NAME_POS		1
#define	OBJ_PRINT_OPS_POS	2
#define OBJ_POS_CSYS_POS	3
#define OBJ_MIN_NUM_PARAMS	2

/*---------------------------------------
 INPUT FILE FORMAT FOR DEFINING ELEMENTS
----------------------------------------*/
#define E_CHAR				'e'
#define E_CHAR_POS			0
#define E_NUM_VERT_POS		1
#define E_BC_CSYS_POS		2
#define E_BC_TYPE_POS		3
#define E_BC_POS			4
#define E_VERTEX_POS		7
#define E_MIN_NUM_PARAMS	10

/*------------------------------------
 FORMAT FOR PRINTING ELEMENT GEOMETRY
-------------------------------------*/
#define ELT_GEOM_LABELS \
	" ELT Vertex Name                  X1         X2         X3\n"
#define ELT_GEOM_UNDLNS \
	"---- -------------------- ---------- ---------- ----------\n"
#define ELT_GEOM_FMT \
	"%4d %-20s %10.4f %10.4f %10.4f\n"


/*-----------------------------------------
 FORMAT FOR PRINTING OBSERVATION GRID DATA
------------------------------------------*/
#define OG_LOC_LABELS \
	"      X1       X2       X3 "
#define OG_LOC_UNDLNS \
	"-------- -------- -------- "
#define OG_LOC_FMT \
	"%8.3f %8.3f %8.3f "
#define OG_DISPL_TITLE \
	"\nDISPLACEMENTS:\n\n"
#define OG_DISPL_LABELS \
	"        U1         U2         U3\n"
#define OG_DISPL_UNDLNS \
	"---------- ---------- ----------\n"
#define OG_DISPL_FMT \
	"%10.3e %10.3e %10.3e\n"
#define OG_STRAIN_TITLE \
	"\nSTRAINS:\n\n"
#define OG_STRAIN_LABELS \
	"       E11        E22        E33        E12        E23        E13\n"
#define OG_STRAIN_UNDLNS \
	"---------- ---------- ---------- ---------- ---------- ----------\n"
#define OG_STRAIN_FMT \
	"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"
#define OG_STRESS_TITLE \
	"\nSTRESSES:\n\n"
#define OG_STRESS_LABELS \
	"     SIG11      SIG22      SIG33      SIG12      SIG23      SIG13\n"
#define OG_STRESS_UNDLNS \
	"---------- ---------- ---------- ---------- ---------- ----------\n"
#define OG_STRESS_FMT \
	"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"
#define OG_PSTRAIN_TITLE \
	"\nPRINCIPAL STRAINS:\n\n"
#define OG_PSTRAIN_LABELS \
	"    N1     N2     N3         E1     N1     N2     N3         E2     N1     N2     N3         E3\n"
#define OG_PSTRAIN_UNDLNS \
	"------ ------ ------ ---------- ------ ------ ------ ---------- ------ ------ ------ ----------\n"
#define OG_PSTRAIN_FMT \
	"%6.3f %6.3f %6.3f %10.3e %6.3f %6.3f %6.3f %10.3e %6.3f %6.3f %6.3f %10.3e\n"
#define OG_PSTRESS_TITLE \
	"\nPRINCIPAL STRESSES:\n\n"
#define OG_PSTRESS_LABELS \
	"    N1     N2     N3       SIG1     N1     N2     N3       SIG2     N1     N2     N3       SIG3\n"
#define OG_PSTRESS_UNDLNS \
	"------ ------ ------ ---------- ------ ------ ------ ---------- ------ ------ ------ ----------\n"
#define OG_PSTRESS_FMT \
	"%6.3f %6.3f %6.3f %10.3e %6.3f %6.3f %6.3f %10.3e %6.3f %6.3f %6.3f %10.3e\n"


/*-------------------------------
 FORMAT FOR PRINTING OBJECT DATA
--------------------------------*/
#define OBJ_LOC_LABELS \
	" ELT      X1C      X2C      X3C "
#define OBJ_LOC_UNDLNS \
	"---- -------- -------- -------- "
#define OBJ_LOC_FMT \
	"%4d %8.3f %8.3f %8.3f "
#define OBJ_DISPL_TITLE \
	"\nDISPLACEMENTS:\n\n"
#define OBJ_DISPL_LABELS \
	"        B1      U1(+)      U1(-)         B2      U2(+)      U2(-)         B3      U3(+)      U3(-) " 
#define OBJ_DISPL_UNDLNS \
	"---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- "
#define OBJ_DISPL_FMT \
	"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e "
#define OBJ_STRESS_TITLE \
	"\nSTRESSES (TRACTIONS):\n\n"
#define OBJ_STRESS_LABELS \
	"        T1         T2         T3 "
#define OBJ_STRESS_UNDLNS \
	"---------- ---------- ---------- "
#define OBJ_STRESS_FMT \
	"%10.3e %10.3e %10.3e "
#define OBJ_BC_CSYS_LABELS \
	"Coord Sys\n"
#define OBJ_BC_CSYS_UNDLNS \
	"---------\n"
#define OBJ_BC_CSYS_FMT \
	"%s\n"


/********************************** Macros **********************************/
#define RADIANS(A) ((A)*PI/180.)		/* Convert degrees to radians		*/
#define MAX(A,B) (((A) > (B)) ? (A):(B))/* MAX macro						*/


/******************************** Structures ********************************/
struct csys_s {							/* -- COORDINATE SYSTEM STRUCTURE -	*/
	char			*name;				/* Coordinate system name			*/
	double			origin[3];			/* Coord sys origin	(global)		*/
	double			local_rot[3][3];	/* (To) global rotation matrix		*/
	struct csys_s	*next;				/* Ptr to next c.s. in linked list	*/
};
typedef struct csys_s csys_t;

struct obs_grid_s {						/* -- OBSERVATION GRID STRUCTURE --	*/
	char			*name;				/* Observation grid name			*/
	int				dimension;			/* Dimension of observation grid	*/
	double			begin[3];			/* Obs grid beginning coords		*/
	double			end[3];				/* Obs grid ending coords			*/
	int				numpts[3];			/* No of obs pts along x1,x2,x3		*/
	int				print[NUM_PR_OPTS];	/* Print options array				*/
	csys_t			*endpt_csys;		/* Input coordinate system			*/
	csys_t			*obspt_csys;		/* Observation grid coord system	*/
	csys_t			*outp_csys;			/* Output coord sys for obs grid	*/
	struct obs_grid_s *next;			/* Ptr to next o.l. in linked list	*/
};
typedef struct obs_grid_s obs_grid_t;

struct vert_s {							/* ------- VERTEX STRUCTURE -------	*/
	char			*name;				/* Vertex name						*/
	double			x[3];				/* Vertex coordinates (global)		*/
	csys_t			*csys;				/* Coordinate system for vertex		*/
	struct vert_s	*next;				/* Ptr to next vertex in linked list*/
};
typedef struct vert_s vert_t;

struct disloc_seg_s {					/* --- DISLOC SEGMENT STRUCTURE ---	*/
	double 			elt_b[3][3];		/* Proj of element b to segment b	*/
	double			trend;				/* Strike of plunging leg of d.s.	*/
	double			plunge;				/* Plunge of plunging leg of d.s.	*/
	double			local_rot[3][3];	/* Local-to-global rotation matrix	*/
	vert_t			*vert[2];			/* Dislocation segment vertices		*/
};
typedef struct disloc_seg_s disloc_seg_t;

struct elt_s {							/* ------ ELEMENT STRUCTURE -------	*/
	int				num_vertices;		/* Number of vertices				*/
	int				bc_type[3];			/* Boundary condition type array	*/
	double			bc[3];				/* Boundary condition magnitudes	*/
	csys_t			elt_csys;			/* Element-local coordinate system	*/
	double			*b[3];				/* Burgers vector array				*/
	disloc_seg_t	 *disloc_seg;		/* Dislocation segment array		*/
	csys_t			*bc_csys;			/* Ptr to coord sys for element BCs	*/
	struct elt_s	*next;				/* Ptr to next elt in linked list	*/
};
typedef struct elt_s elt_t;

struct obj_s {							/* ------ OBJECT STRUCTURE -------- */
	char			*name;				/* Object type						*/
	int				print[NUM_PR_OPTS];	/* Print options					*/
	csys_t			*pos_csys;			/* Position coordinate system		*/
	elt_t			*first_elt;			/* Pointer to first element			*/
	elt_t			*last_elt;			/* Pointer to last element			*/
	struct obj_s	*next;				/* Ptr to next obj in linked list	*/
};
typedef struct obj_s obj_t;


/*************************** External Variables *****************************/
char		*title1_E = NULL;			/* Problem title					*/
char		*title2_E = NULL;			/* Problem subtitle					*/
int			half_space_E = TRUE;		/* Half/whole space flag			*/
int			check_cond_num_E = TRUE;	/* Check matrix condition num flag	*/
double		cond_num_E = -1.0;			/* Matrix condition number			*/
char		infile_E[MAXFILE];			/* Input file name					*/
char		outfile_E[MAXFILE];			/* Output file name					*/
int			linenum_E = 0;				/* Current line # in input file		*/
int			num_elts_E = 0;				/* Number of elements				*/
int			below_vertex_E = FALSE;		/* Below vertex flag				*/
double		null_value_E = -999.0;		/* Null output value				*/

/*********************************************************************************************************/
/************************* NEW: added 98-12-09 to reflect the problem of observation point near vertices */
double		coef_exclu_E   = 0.0;		/* coef exclusion value */
int			near_vertex_E = FALSE;		/* nearnest vertex flag for observation point */
/*       EXPLANATIONS: 
    Let's obs_pt be an observation point.
    Put flag near_vertex_E to FALSE.
    For each vertex in current project:
    Let's d be the mean length of segments containing v
    Let's l be the distance from v to a choosen observation point
    Then if l<d*coef_exclu => near_vertex_E = TRUE
    
    Then, before printing computed values for this observation point:
      if near_vertex_E=TRUE => print null_value_E
      otherwise print computed values
      
    So coef_exclu_E = 1.0, means 100% of the mean length of segments containing v.
*/
/*********************************************************************************************************/

FILE		*ifp_E	= stdin;			/* Input file ptr (default stdin)	*/
FILE		*ofp_E	= stdout;			/* Output file ptr (default stdout) */
FILE		*tempfp_E[NUM_PR_OPTS];		/* Temporary file ptrs				*/

double		shear_mod_E		= -1.0;		/* Shear modulus					*/
double		psn_ratio_E		= -1.0;		/* Poisson's ratio					*/
double		youngs_mod_E	= -1.0;		/* Young's modulus					*/
double		bulk_mod_E		= -1.0;		/* Bulk modulus						*/
double		lame_lambda_E	= -1.0;		/* Lame's lambda					*/

int			print_elt_geom_E = FALSE;	/* Print element geometry flag		*/
char		*elt_geom_csys_name_E =NULL;/* Element geometry coord sys name	*/

int			rem_stress_bc_E = TRUE;		/* Remote stress vs strain bc flag	*/
double		rem_stress_E[3][3];			/* Remote stress tensor				*/
double		rem_strain_E[3][3];			/* Remote strain tensor				*/

double		*b_vector_E;				/* Burger's vector array			*/
double		**ic_matrix_E;				/* Influence coeff matrix			*/

csys_t		*first_csys_E     = NULL;	/* 1st memb of csys linked list		*/
obs_grid_t	*first_obs_grid_E = NULL;	/* 1st memb of obs grid linked list	*/
obj_t		*first_obj_E      = NULL;	/* 1st memb of obj linked list		*/
elt_t		*first_elt_E      = NULL;	/* 1st memb of elt linked list		*/
vert_t		*first_vert_E     = NULL;	/* 1st memb of vert linked list	*/


/************************* ANSI Function Declarations ***********************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
double	array_max_norm(double **a, int start_row, int end_row, int start_col,
		int	end_col);
int		calc_elt_parameters(elt_t *current_elt);
void	close_temp_files(void);
void	copy_temp_files(void);
void	determine_burgers_vectors(void);
void	displ_strain(int calc_displ, int calc_strain, double x[3],
		double displ[3], double strain[3][3], elt_t *omit_elt);
void	displ_strain_poly_elt(int calc_displ, int calc_strain, 
		elt_t *current_elt, double x[3], double displ[3],
		double strain[3][3], elt_t *omit_elt, int under);
void	print_obj_data(void);
csys_t	*find_csys(char *name);
vert_t	*find_vert(char *name);
int		get_double_var(double *var, char *var_name, char *word[],
		int numwords);
int		get_boolean_var(int *var, char *var_name, char *true_string,
		char *false_string, char *word[], int numwords, char *line);
int		get_text_var(char **var, char *var_name, char *word[], int numwords,
		char *line);
void	get_elt_info(elt_t **current_elt, obj_t *current_obj, int numwords,
		char *word[], char *line);
void	get_obj_info(obj_t **current_obj, int numwords, char *word[],
		char *line);
void	get_program_args(void);
void	get_vert_info(vert_t **current_vert, int numwords, char *word[],
		char *line);
void	displ_strain_ics_poly_elt(int calc_displ, int calc_strain, 
		elt_t *current_elt, double x[3],
		double displ_ic[3][3], double strain_ic[3][3][3],
		elt_t *omit_elt);
void	print_obs_grid_data(void);
int		open_files();
void	open_temp_files(int print[]);
int		parse_command_line_args(int argc, char *argv[]);
void	p_error(char *error_msg, char *line);
void	display_msg(char *_msg);
void	print_elt_data(obj_t *current_obj,
		elt_t *current_elt, int elt_num, double displ[3],
		double  stress[3][3]);
void	print_elt_geometry(void);
void	print_obj_titles(obj_t *current_obj);
void	print_obs_grid_titles(obs_grid_t *current_obs_grid);
void	print_obs_pt_data(obs_grid_t *current_obs_grid, double x[3],
		double displ[3], double strain[3][3]);
void	print_problem_info(void);
int		read_csystems();
int		read_objs_elts_verts();
void	read_infile(void);
int		read_line(char *line, char *word[]);
int		read_observation_grids();
int		read_constants();
void	setup_global_coords(void);


/************************** K&R Function Declarations ***********************/
#else
double	array_max_norm();
int		calc_elt_parameters();
void	close_temp_files();
void	copy_temp_files();
void	determine_burgers_vectors();
void	displ_strain();
void	displ_strain_poly_elt();
void	print_obj_data();
csys_t	*find_csys();
vert_t	*find_vert();
int		get_double_var();
int		get_boolean_var();
int		get_text_var();
void	get_elt_info();
void	get_obj_info();
void	get_program_args();
void	get_vert_info();
void	displ_strain_ics_poly_elt();
void	print_obs_grid_data();
int		open_files();
void	open_temp_files();
int		parse_command_line_args();
void	p_error();
void	display_msg();
void	print_elt_data();
void	print_elt_geometry();
void	print_obj_titles();
void	print_obs_grid_titles();
void	print_obs_pt_data();
void	print_problem_info();
int		read_csystems();
int		read_objs_elts_verts();
void	read_infile();
int		read_line();
int		read_observation_grids();
int		read_constants();
void	setup_global_coords();
#endif


/******************************* Function: main ******************************
* In:	argc	- number of command line arguments
*		argv	- array of command line arguments
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
         main(int argc, char *argv[])
#else
main(argc, argv)
int	 argc;
char *argv[];
#endif

{
	time_t   start, t_reading, t_matrix, t_problem, t_object, t_obs,finish; 
	double   elapsed_time;
	/* Use get_program_args() or parse_command_line_args() to get
	   file names and program options.
	-------------------------------------------------------------*/
#if FPROMPT
	get_program_args();
#else
	parse_command_line_args(argc,argv);
#endif
 
	/* Start counter */ 
	time( &start ); 

	/* Open the input and output files
	----------------------------------*/
	open_files();

	/* Read input file and set up the problem
	-----------------------------------------*/
	read_infile(); 
 
	time( &t_reading ); 
	elapsed_time = difftime( t_reading, start ); 
	printf( "\nReading file       : %6.0f seconds.", elapsed_time );

	/* Solve for burger's vector for each element
	---------------------------------------------*/ 
	determine_burgers_vectors(); 
 
	time( &t_matrix ); 
	elapsed_time = difftime( t_matrix,t_reading ); 
	printf( "\nEquations/matrix   : %6.0f seconds.", elapsed_time ); 

	/* Print the problem info
	-------------------------*/
	print_problem_info(); 
 
	time( &t_problem ); 
	elapsed_time = difftime( t_problem,t_matrix ); 
	printf( "\nPrint pb info      : %6.0f seconds.", elapsed_time );

	/* Calculate displs and tractions on elements
	---------------------------------------------*/ 
	print_obj_data(); 
 
	time( &t_object ); 
	elapsed_time = difftime( t_object,t_problem ); 
	printf( "\nD and T on elements: %6.0f seconds.", elapsed_time ); 

	/* Calculate displacements and stresses along observation grids
	---------------------------------------------------------------*/ 
	print_obs_grid_data(); 
 
	time( &t_obs ); 
	elapsed_time = difftime( t_obs,t_object ); 
	printf( "\nD and S on obs     : %6.0f seconds.", elapsed_time ); 
 
	/* End counter, and print elapsed time */ 
	time( &finish ); 
	elapsed_time = difftime( finish, start ); 
	printf( "\nProgram takes %6.0f seconds.\n", elapsed_time ); 

	return(0);
}


/*********************** Function: print_problem_info ***********************
* Prints general problem information to the output file.
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_problem_info(void)
#else
void print_problem_info()
#endif

{
	/* Print program name, version, and date
	----------------------------------------*/
	fprintf(ofp_E,"OUTPUT FROM: %s, version %s\n",PROGRAM, VERSION);
	fprintf(ofp_E,"   COMPILED: %s\n",COMPILE_DATE);

	/* Print input file name and problem titles
	-------------------------------------------*/
	fprintf(ofp_E,"\n INPUT FILE: %s\n",infile_E);
	fprintf(ofp_E,  "     TITLE1: %s\n",title1_E);
	fprintf(ofp_E,  "     TITLE2: %s\n",title2_E);

	/* Print elastic constant values
	--------------------------------*/
	fprintf(ofp_E,"\nELASTIC CONSTANTS:\n");
	fprintf(ofp_E,  "    Shear Modulus   = %f\n",shear_mod_E);
	fprintf(ofp_E,  "    Poisson's Ratio = %f\n",psn_ratio_E);
	fprintf(ofp_E,  "    Young's Modulus = %f\n",youngs_mod_E);
	fprintf(ofp_E,  "    Bulk Modulus    = %f\n",bulk_mod_E);
	fprintf(ofp_E,  "    Lame's Lambda   = %f\n",lame_lambda_E);

	/* Print the null output value
	------------------------------*/
	fprintf(ofp_E,"\nNULL OUPUT VALUE = %f\n",null_value_E);

   /* Print the coef exclusion value
	------------------------------*/
	fprintf(ofp_E,"\nCOEF EXCLUSION VALUE = %f\n",coef_exclu_E);

	/* Print condition number of the influence coefficient matrix
	-------------------------------------------------------------*/
	fprintf(ofp_E,"\nCONDITION NUMBER = ");
	if (cond_num_E < 0.0) {
		fprintf(ofp_E,"(no traction bc's -> no matrix needed)\n");
	} else {
		if (check_cond_num_E)
			fprintf(ofp_E,"%f\n",cond_num_E);
		else
			fprintf(ofp_E,"(not requested)\n");
	}

	/* Print element geometries (if requested)
	------------------------------------------*/
	if (print_elt_geom_E)
		print_elt_geometry();
}


/************************** Function: print_obj_data ************************
* Calculates and prints object data to the output file.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_obj_data(void)
#else
void print_obj_data()
#endif

{

	int			elt_num;
	int			calc_displ;
	int			calc_strain;
	obj_t		*current_obj;
	double		displ[3];
	double		strain[3][3];
	double		stress[3][3];
	elt_t		*current_elt;

	current_obj = first_obj_E;
	current_elt = first_elt_E;
	while (current_elt != NULL) {

		/* If this element starts new object......
		-----------------------------------------*/
		if (current_elt == current_obj->first_elt) {

			/* Reset element number
			----------------------*/
			elt_num = 1;

			if (current_obj->print[DISPL] || current_obj->print[STRESS]) {

				/* Open temp files
				------------------*/
				open_temp_files(current_obj->print);
	
				/* Print object titles
				----------------------*/
				print_obj_titles(current_obj);

			} else {

				/* Skip this object if no output requested
				------------------------------------------*/
				current_elt = current_obj->last_elt;
			}
		}

		if (current_obj->print[DISPL] || current_obj->print[STRESS]) {

			calc_displ  = current_obj->print[DISPL];
			calc_strain = current_obj->print[STRESS];
			displ_strain(calc_displ,calc_strain,
				current_elt->elt_csys.origin,displ,strain,
				current_elt);
			strain_to_stress(strain,shear_mod_E,lame_lambda_E,
				stress);

			// Print the data for this elt
			print_elt_data(current_obj, current_elt, elt_num,
				displ, stress);

			elt_num++;
		}

		if (current_elt == current_obj->last_elt) {
			if (current_obj->print[DISPL] || current_obj->print[STRESS]) {
				
				/* Copy temp files to main output file, then close them.
				--------------------------------------------------------*/
				copy_temp_files();
				close_temp_files();
			}
			current_obj = current_obj->next;
		}

		current_elt = current_elt->next;
	}
}


/********************** Function: print_obs_grid_data ******************
* Loops through linked list of observation grids, caculating and printing
* the requested displacement, strain, and stress data for each to the
* output file.
************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_obs_grid_data(void)
#else
void print_obs_grid_data()
#endif
{
	obs_grid_t	*current_obs_grid;
	int		i, j, k;
	double	x[3];
	double	dx[3];
	double	displ[3];
	double	strain[3][3];
	double	*begin;
	double	*end;
	int		*numpts;
	char	   error_msg[MAX_ERROR_MSG];
	int		*print;
	int		calc_displ;
	int		calc_strain;

   /* Declaration for the bug correction */
   double   x_copy[3];


	/* Loop over each observation grid
	----------------------------------*/
	current_obs_grid = first_obs_grid_E;
	while (current_obs_grid != NULL) {

		/* Determine data required by obs grid print options
		----------------------------------------------------*/
		print = current_obs_grid->print;
		calc_displ = print[DISPL];
		calc_strain = (print[STRAIN] || print[PSTRAIN] || print[STRESS]
			|| print[PSTRESS]);

      /* Remenber that begin_pt & end_pt are in GLOBAL coordinate system.*/
		begin  = current_obs_grid->begin;
		end    = current_obs_grid->end;
		numpts = current_obs_grid->numpts;
		subtract_vectors(end,begin,dx);

		/* Open temp files
		------------------*/
		open_temp_files(current_obs_grid->print);

		/* Print observation grid titles
		--------------------------------*/
		print_obs_grid_titles(current_obs_grid);

		/* Process the observation grid
		-------------------------------*/
		switch (current_obs_grid->dimension) {

			case 0:
				copy_vector(begin,x);
            /* Transform current observation point into global CSys
            ----------------------------------------------------*/
            transform_position_vector(INVERSE_TRANS,
                                      current_obs_grid->endpt_csys->origin,
                                      current_obs_grid->endpt_csys->local_rot,
                                      x);
				displ_strain(calc_displ,calc_strain,x,displ,strain,
					NULL);
				print_obs_pt_data(current_obs_grid,x,displ,strain);
				break;

			case 1:
				dx[2] /= (numpts[0]-1);
				dx[1] /= (numpts[0]-1);
				dx[0] /= (numpts[0]-1);
				for (i=0; i < numpts[0]; i++) {
					x[2] = begin[2] + dx[2]*i;
					x[1] = begin[1] + dx[1]*i;
					x[0] = begin[0] + dx[0]*i;
               /* Transform current observation point into global CSys
               ----------------------------------------------------*/
               transform_position_vector(INVERSE_TRANS,
                                         current_obs_grid->endpt_csys->origin,
                                         current_obs_grid->endpt_csys->local_rot,
                                         x);
					displ_strain(calc_displ,calc_strain,x,displ,
						strain,NULL);
					print_obs_pt_data(current_obs_grid,x,displ,strain);
				}
				break;

			case 2:
			case 3:
				for (i=0; i < 3; i++) {
					dx[i] /= (numpts[i] == 1) ? 1 : (numpts[i]-1);
				}
				for (i=0; i < numpts[2]; i++) {
					x[2] = begin[2] + dx[2]*i;
					for (j=0; j < numpts[1]; j++) {
						x[1] = begin[1] + dx[1]*j;
						for (k=0; k < numpts[0]; k++) {
							x[0] = begin[0] + dx[0]*k;
                     copy_vector(x,x_copy);
                     /* Transform current observation point into global CSys
                     ----------------------------------------------------*/
                     transform_position_vector(INVERSE_TRANS,
                                               current_obs_grid->endpt_csys->origin,
                                               current_obs_grid->endpt_csys->local_rot,
                                               x_copy);
							displ_strain(calc_displ,calc_strain,/****/x_copy/****/,displ,strain,NULL);
							print_obs_pt_data(current_obs_grid, /****/x_copy/****/,displ,strain);
						}
					}
				}
				break;

			default:
				sprintf(error_msg,
					"Invalid dimension (%d) for observation grid",
					current_obs_grid->dimension);
				p_error(error_msg,NULL);

		}

		/* Copy temp files to main output file, then close them.
		--------------------------------------------------------*/
		copy_temp_files();
		close_temp_files();

		current_obs_grid = current_obs_grid->next;

	} /*while*/
}


/************************ Function: displ_strain ****************************
* Calculates the total displacement and/or strain at a point due to ALL
* elements. 
*
* In:	   calc_displ	- calculate displacements flag
*		   calc_strain - calculate strains flag
*		   x			   - coords (global) of pt at which to calc displ & strain
*		   omit_elt	   - element to omit when calculating displs (NULL = none)
*
* Out:	displ  		- displacement vector (global coords)
*		   strain 		- strain tensor (global coords)
*
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     displ_strain(int calc_displ, int calc_strain, double x[3],double displ[3], double strain[3][3], elt_t *omit_elt)
#else
void displ_strain(calc_displ, calc_strain, x, displ, strain, omit_elt)
int				calc_displ;
int				calc_strain;
double			x[3];
double			displ[3];
double			strain[3][3];
elt_t		      *omit_elt;
#endif

{
	elt_t	*current_elt;
	double	elt_displ[3];
	double	elt_strain[3][3];
/*Declarations for correction of the "shadow effect" */
         int i;
         double orient;
      	double under_plane;
      	int under;
      	double normal [3];
      	double data1 [3];
      	double data [3];
      	vert_t *verta;
      	vert_t *vertb;
      	vert_t *vertc;
      	double seg1 [3];
      	double seg2 [3];
      	double inside;
      	double inside_vector [3];
      	int inside_test;
      	int inside_test_fin;
      	double x3_global [3];
      	disloc_seg_t *disloc_seg;
/* end of the declaration */

/*Declarations for coef_exclu */
         double dist;
/* end of the declaration */

	initialize_vector(displ,0.0);
	initialize_matrix(strain,0.0);

	/* Loop over each element
	-------------------------*/
	current_elt = first_elt_E;
   near_vertex_E = FALSE;
	while (current_elt != NULL)
   {
		/* Test UNDER to determine if the data point is inside the "shadow zone"
		---------------------------------------------------------------------*/
		under = 0;
		inside_test_fin = 1;
		
		/* Determine if the element has a positive side up or down (orient)
		   and determine if the data is under the plane 
			defined by the element (under_plane) 
		-----------------------------------------------------------------*/
		disloc_seg = current_elt->disloc_seg;
		x3_global [0] = 0;
		x3_global [1] = 0;
		x3_global [2] = 1;
		verta = disloc_seg [0].vert[0];
		vertb = disloc_seg [0].vert[1];
		vertc = disloc_seg [1].vert[1]; 
		subtract_vectors (vertb->x, verta->x, seg1);
		subtract_vectors (vertc->x, verta->x, seg2);
		cross_product (seg1, seg2, normal);
		orient = dot_product (normal, x3_global);
		subtract_vectors (x, verta->x, data1);
		under_plane = dot_product (normal, data1);

      /* Determine if obs_pt is nearnest vertex verta, vertb or vertc
         1) compute d = mean length of the 3 disloc_seg
         2) compute distance from obs_pt to verta and see if this diance is < d*coef_exclu
         3) compute distance from obs_pt to vertb and see if this diance is < d*coef_exclu
         4) compute distance from obs_pt to vertc and see if this diance is < d*coef_exclu
         if condition is TRUE, break the loop over elements and set the flag near_vertex_E=TRUE
            and return.
            
         We use a new function "distance" define in matrix.h & matrix.c
      */
      dist = 0.0;
      dist  = distance(verta->x,vertb->x);
      dist += distance(vertb->x,vertc->x);
      dist += distance(vertc->x,verta->x);
      dist *= coef_exclu_E/3.0;
      
      if (distance(x,verta->x)<=dist || distance(x,vertb->x)<=dist || distance(x,vertc->x)<=dist)
      {
         near_vertex_E = TRUE;
         break;
      }
		
		/* Determine if the data point is in the "rigid body" 
			(cf. explanation of the bug)
		---------------------------------------------------*/

		for (i = 0; i < current_elt->num_vertices; i++)
		{
			inside = 0;
			verta = disloc_seg [i].vert[0];
			vertb = disloc_seg [i].vert[1];
			subtract_vectors (vertb->x, verta->x, seg1);
		   subtract_vectors (x, verta->x, data);
			cross_product (seg1, data, inside_vector);
			inside = dot_product (x3_global, inside_vector);
			
			if (orient > 0) 
			{
               if (inside > 0 && under_plane < 0)
                  inside_test = 1;
				   else
                  inside_test = 0;
			}
			if (orient < 0)
			{
            if (inside < 0 && under_plane >  0)
               inside_test = 1;
				else
               inside_test = 0;
			}
			if (orient == 0)
				under = 0;

			inside_test_fin *= inside_test;  
		}

		/*Gives a value 1 to under for data points under the element
		and positive side up
		Gives a value 2  to under for data points under the element
		and positive side down
		Gives a value 0  to under for data points not under the element
		---------------------------------------------------------------*/

		if (inside_test_fin > 0)
      {
			if (orient > 0)
            under = 1;
		   if (orient < 0)
            under = 2;
		}
      else
         under = 0;

		/* Calculate displacement and strain due to this element
		--------------------------------------------------------*/
		displ_strain_poly_elt(calc_displ,calc_strain,current_elt,
			x,elt_displ,elt_strain,omit_elt,under);

		/* Add the contribution of this element to the total displ & strain
		-----------------------------------------------------------------*/
		add_vectors(displ,elt_displ,displ);
		add_matrices(strain,elt_strain,strain);
		current_elt = current_elt->next;
	}

	/* Adjust strain for remote strain
	----------------------------------*/
   if (near_vertex_E!=TRUE)
	   add_matrices(strain,rem_strain_E,strain);
}


/****************** Function: displ_strain_poly_elt ********************
* Calculates the displacement and/or strain at a point due to a single
* polygonal element.
*
* In:	calc_displ	- calculate displacements flag
*		calc_strain - calculate strains flag
*		current_elt	- element being considered
*		x			- coords (global) of pt at which to calc displ & strain
*		omit_elt	- element to omit when calculating displs (NULL = none)
* Out:	elt_displ  	- displacement vector (global coords) due to this elt
*		elt_strain	- strain tensor (global coords) due to this elt
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     displ_strain_poly_elt(int calc_displ, int calc_strain,elt_t *current_elt, double x[3], double elt_displ[3],double elt_strain[3][3], elt_t *omit_elt, int under)
#else
void displ_strain_poly_elt(calc_displ, calc_strain, current_elt, x,
	elt_displ, elt_strain, omit_elt,under)
int		calc_displ;
int		calc_strain;
elt_t	*current_elt;
//elt_t   **current_elt;
double	x[3][3];
double	elt_displ[3];
double	elt_strain[3][3];
elt_t	*omit_elt;
int under;
#endif

{
	double	displ_ic[3][3];
	double	strain_ic[3][3][3];
	int		i, j, k;
	/*Declarations for shadow effect correction*/
	double bglobal[3];
	double e[3][3];
	double eg[3][3];
	double bg[3][3];

	double value=0.0;
	/* Initialize displacement vector and strain tensor
	---------------------------------------------------*/
	initialize_vector(elt_displ,0.0);
	initialize_matrix(elt_strain,0.0);

	/* Calculate the displacement and strain influence coefficients
	---------------------------------------------------------------*/
	displ_strain_ics_poly_elt(calc_displ,calc_strain,current_elt,x,
		displ_ic,strain_ic,omit_elt);

	/* Superpose the contribution from each burger's vector component
	-----------------------------------------------------------------*/
	for (i=0; i < 3; i++)
   {
		for (j=0; j < 3; j++)
      {
			elt_displ[i] += displ_ic[j][i] * (*current_elt->b[j]); 	
         for (k=0; k < 3; k++)
				elt_strain[i][j] += strain_ic[k][i][j] * (*current_elt->b[k]);
		}
	}

 /* Shadow effect correction in case of data point under the element
---------------------------------------------------------------------------------*/
		
	if (under > 0)
   {
	   bglobal[0]=0;
      bglobal[1]=0;
      bglobal[2]=0;

		e[0][0]=1;e[0][1]=0;e[0][2]=0;
		e[1][0]=0;e[1][1]=1;e[1][2]=0;
		e[2][0]=0;e[2][1]=0;e[2][2]=1;

		eg[0][0]=1;eg[0][1]=0;eg[0][2]=0;
		eg[1][0]=0;eg[1][1]=1;eg[1][2]=0;
		eg[2][0]=0;eg[2][1]=0;eg[2][2]=1;

		for (i=0; i < 3; i++)
      {
			/* Transform the burger vector's component in the global coordinate system
			--------------------------------------------------------------------------- */
			rotate_vector(INVERSE_ROT,current_elt->elt_csys.local_rot,e[i]);
			scalar_vector_mult (*current_elt->b[i], e[i],bg[i]);
		}
		for (i=0; i < 3; i++)
      {
		   for (j=0; j<3; j++)
            bglobal[i]+=dot_product(bg[j], eg[i]);

			/*Corrects the displacement by the corresponding burger's vector (global coordsys)
			--------------------------------------------------------------------------------*/ 
			if (under == 1)
			   elt_displ[i] -= bglobal[i];
         else
			if (under == 2)
				elt_displ[i] += bglobal[i];
		}

	}/* end of if (under>0)*/
}

/***************** Function: displ_strain_ics_poly_elt **********************
* Calculates the displacement and/or strain influence coefficients at a point 
* due to a polygonal element.
*
* In:	calc_displ	- calculate displacements flag
*		calc_strain - calculate strains flag
*		current_elt	- element being considered
*		x			- coords (global) of pt at which to calc displ & strain
*		omit_elt	- element to omit when calculating displs (NULL = none)
*
* Out:	displ_ic[i][j]		- the jth component of displ (global coords)
*                        	  due to a unit ith Burgers vector component 
*                        	  (bc_coord_sys coords)
*		strain_ic[i][j][k]	- the jk component of strain (global coords)
*                        	  due to a unit ith Burgers vector component 
*                        	  (bc_coord_sys coords)
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     displ_strain_ics_poly_elt(int calc_displ, int calc_strain,elt_t *current_elt, double x[3], double displ_ic[3][3],double strain_ic[3][3][3], elt_t *omit_elt)
#else
void displ_strain_ics_poly_elt(calc_displ, calc_strain, current_elt, x, 
	displ_ic, strain_ic, omit_elt)
int		calc_displ;
int		calc_strain;
elt_t	*current_elt;
double	x[3];
double	displ_ic[3][3];
double	strain_ic[3][3][3];
elt_t	*omit_elt;
#endif

{
	int		i, j, k, l;
	int		seg;
	int		swap;
	vert_t	 *vert1;
	vert_t	 *vert2;
	double	depth1;
	double	depth2;
	double	beta;
	double	temp_double;
	double	temp_vector[3];
	double	r;
	double	z3;
	double	y1[3];
	double	y2[3];
	double	displ_ic1[3][3];
	double	displ_ic2[3][3];
	double	displ_ic3[3][3];
	double	displ_ic4[3][3];
	double	strain_ic1[3][3][3];
	double	strain_ic2[3][3][3];
	double	strain_ic3[3][3][3];
	double	strain_ic4[3][3][3];
	disloc_seg_t	 *disloc_seg;
	
	/* Initialize the influence coefficients
	----------------------------------------*/
	for (i=0; i < 3; i++) {
		initialize_vector(displ_ic[i],0.0);
		initialize_matrix(strain_ic[i],0.0);
	}

	/* Loop over the element's dislocation segments
	-----------------------------------------------*/
	disloc_seg = current_elt->disloc_seg;
	for (seg = 0; seg < current_elt->num_vertices; seg++) {

		/* Determine the segment vertices
		---------------------------------*/
		vert1 = disloc_seg[seg].vert[0];
		vert2 = disloc_seg[seg].vert[1];

		/* Compute the vectors from the segments vertices to the
                   data point
        Matrix<double> INFLUENCE(num_eqns, num_eqns);
        Vector<int> PIVOT;
        Vector<double> DISP;
        Vector<double> B_VECTOR(num_eqns);
        Vector<double> X;
        
        for (i=0; i<=num_eqns-1; i++){
		for (j=0; j<=num_eqns-1; j++){
                    INFLUENCE(j,i) = ic_matrix[j+1][i+1];
                }
        }
        printf("\nConverted influence matrix");
        
        for (i=0; i<=num_eqns-1; i++){
                B_VECTOR(i) = b_vector[i+1];
            
        }
        printf("\nConverted displacement vector");
        
        // Old code (nrutils)
        // d_ludcmp(ic_matrix, num_eqns, pivot, &d);
	// d_lubksb(ic_matrix, num_eqns, pivot, b_vector);

        printf("\nStarted LU decomposition");
        GetLU(INFLUENCE, PIVOT);
        printf("\nGot LU decomposition");
        X = B_VECTOR;
        SolveLU(INFLUENCE, PIVOT, X);
        printf("\nSolved system of equations");
        B_VECTOR = X;
        
        // Load results back into poly3D variables
         for (i=0; i<=num_eqns-1; i++){
		for (j=0; j<=num_eqns-1; j++){
                    ic_matrix[j+1][i+1] = INFLUENCE(j,i); 
                }
        }

        for (i=0; i<=num_eqns-1; i++){
                b_vector[i+1] = B_VECTOR(i);
            
        }
		---------------------------------------------------------*/
		subtract_vectors(x,vert1->x,y1);
		subtract_vectors(x,vert2->x,y2);

		/* Rotate vert-to-obs_point vectors to segment-local coords
		-----------------------------------------------------------*/
		rotate_vector(FORWARD_ROT,disloc_seg[seg].local_rot,y1);
		rotate_vector(FORWARD_ROT,disloc_seg[seg].local_rot,y2);

		depth1 = -vert1->x[2];
		depth2 = -vert2->x[2];
		beta   = PI/2.0 - disloc_seg[seg].plunge;

		if (((sqrt(y1[0]*y1[0]+y1[1]*y1[1]) < BVERT_TINY) && y1[2] >= 0.0) ||
			((sqrt(y2[0]*y2[0]+y2[1]*y2[1]) < BVERT_TINY) && y2[2] >= 0.0)) {
			below_vertex_E = TRUE;
			return;
		}

		/* If x lies along dipping leg of angular dislocations, swap
		   the vertex order, so singularity will be avoided
		------------------------------------------------------------*/
		swap = FALSE;
		z3 = y1[0]*sin(beta) + y1[2]*cos(beta);
		r  = vector_magnitude(y1);
		if ((r - z3) < SWAP_TINY) {
			swap = TRUE;
			copy_vector(y2,temp_vector);
			copy_vector(y1,y2);
			copy_vector(temp_vector,y1);
			for (i=0; i < 2; i++) {
				y1[i] *= -1.0;
				y2[i] *= -1.0;
			}
			temp_double = depth2;
			depth2 = depth1;
			depth1 = temp_double;
			beta = PI - beta;
		} 
 
		/* Be carefull for the calculation of vertical elements 
		------------------------------------------------------*/ 
		if (beta==0.0) 
			beta = 1.0e-14;
			

		/* Calculate displacement influence coeffs
		------------------------------------------*/
		if (calc_displ) {

			/* Avoid displacement discontinuity when calculating displ
			   inf coeff of an element on itself
			----------------------------------------------------------*/
			if  (current_elt != omit_elt) {

				comninou_displ_ics(y1,depth1,beta,psn_ratio_E,
					half_space_E,displ_ic1);
				comninou_displ_ics(y2,depth2,beta,psn_ratio_E,
					half_space_E,displ_ic2);

				/* Superpose the angular dislocation influence coeffs into
		   	   	a dislocation segment influence coeff
				----------------------------------------------------------*/
				subtract_matrices(displ_ic1,displ_ic2,
					displ_ic3);

				/* Swap the vertices back to proper order (if necessary)
				--------------------------------------------------------*/
				if (swap) {
					scalar_vector_mult(-1.0,displ_ic3[2],
						displ_ic3[2]);
					for (i=0; i < 3; i++) {
						displ_ic3[i][0] *= -1.0;
						displ_ic3[i][1] *= -1.0;
					}
				}

				/* Transform from disloc segment to element influence coeffs
				------------------------------------------------------------*/
				for (i=0; i < 3; i++) {
					initialize_vector(displ_ic4[i],0.0);
					for (j=0; j < 3; j++) {
						for (k=0; k < 3; k++) {
							displ_ic4[i][j] +=
								disloc_seg[seg].elt_b[i][k] *
								displ_ic3[k][j];
						}
					}
				}
		
				/* Rotate from C&D to global coordinates
				----------------------------------------*/
				for (i=0; i < 3; i++) {
					rotate_vector(INVERSE_ROT,disloc_seg[seg].local_rot,
						displ_ic4[i]);
				}

				/* Superpose the contribution of this dislocation segment
				---------------------------------------------------------*/
				for (i=0; i < 3; i++) {
					add_vectors(displ_ic[i],displ_ic4[i],
						displ_ic[i]);
				}
	
			} 
		}

		/* Calculate strain influence coeffs
		------------------------------------*/
		if (calc_strain) {

			comninou_strain_ics(y1,depth1,beta,psn_ratio_E,
				half_space_E,strain_ic1);
			comninou_strain_ics(y2,depth2,beta,psn_ratio_E,
				half_space_E,strain_ic2);

			/* Superpose the angular dislocation influence coeffs into
		   	   a dislocation segment influence coeff
			----------------------------------------------------------*/
			for (i=0; i < 3; i++) {
				subtract_matrices(strain_ic1[i],strain_ic2[i],
					strain_ic3[i]);
			}

			/* Swap the vertices back to proper order (if necessary)
			--------------------------------------------------------*/
			if (swap) {
				scalar_matrix_mult(-1.0,strain_ic3[2],
					strain_ic3[2]);
				for (i=0; i < 3; i++) {
					strain_ic3[i][0][2] *= -1;
					strain_ic3[i][2][0] *= -1;
					strain_ic3[i][1][2] *= -1;
					strain_ic3[i][2][1] *= -1;
				}
			}


			for (i=0; i < 3; i++) {
				initialize_matrix(strain_ic4[i],0.0);
				for (j=0; j < 3; j++) {
					for (k=0; k < 3; k++) {
						for (l=0; l < 3; l++) {
							strain_ic4[i][j][k] +=
								disloc_seg[seg].elt_b[i][l] *
								strain_ic3[l][j][k];
						}
					}
				}
			}
		
			/* Rotate strain inf coeffs to global coords
			--------------------------------------------*/
			for (i=0; i < 3; i++) {
				rotate_tensor(INVERSE_ROT,disloc_seg[seg].local_rot, strain_ic4[i]);
			}

			/* Superpose the contribution of this dislocation segment
			---------------------------------------------------------*/
			for (i=0; i < 3; i++) {
				add_matrices(strain_ic[i],strain_ic4[i],
					strain_ic[i]);
			}

		}
	
	} /* loop over dislocation segments */

}


/**************************** Function: read_infile **************************
* Reads the input file and sets up the problem to be solved.
******************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     read_infile(void)
#else
void read_infile()
#endif

{

	/* Read problem constants
	-------------------------*/
	read_constants();

	/* Read coordinate systems
	--------------------------*/
	read_csystems();

	/* Read observation grids
	-------------------------*/
	read_observation_grids();

	/* Read elements and vertices
	-----------------------------*/
	read_objs_elts_verts();

}


/************************** Function: find_csys *************************
* Returns a pointer to the coordinate system named by name, or NULL if no
* such coordinate system exists.
*
* In:	name	- name of the coordinate system to find
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
csys_t*  find_csys(char *name)
#else
csys_t *find_csys(name)
char *name;
#endif

{
	csys_t *current_csys;

	current_csys = first_csys_E;

	while (current_csys != NULL) {
		if (!strcmp(name,current_csys->name))
			break;
		current_csys = current_csys->next;
	}

	return(current_csys);
}


/***************************** Function: find_vert ***************************
* Return a pointer to the vertex named by name, or NULL if no such vertex
* exists.
*
* In:	name	- name of the vertex to find
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
vert_t*  find_vert(char *name)
#else
vert_t *find_vert(name)
char *name;
#endif

{
	vert_t	*current_vert;

	current_vert = first_vert_E;

	while (current_vert != NULL) {
		if (!strcmp(name,current_vert->name))
			break;
		current_vert = current_vert->next;
	}

	return(current_vert);
}


/***************************** Function: p_error *****************************
* Prints an error message to stderr and calls exit().  If line != NULL, the
* line and line number (from the input file) on which the error occured are
* printed as well.
*
* In:	error_msg	- error message to print
* 		line		- input file line number at which the error occurred
*					  (NULL = N/A)
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     p_error(char *error_msg, char *line)
#else
void p_error(error_msg, line)
char	*error_msg;
char	*line;
#endif

{
	fprintf(stderr,"\nerror: %s",error_msg);
	if (line == NULL) {
		fprintf(stderr,"\n");
	} else {
		fprintf(stderr," (%s, line %d)\n",infile_E,linenum_E);
		fprintf(stderr,"       %s\n",line);
	}
	exit(1);

}

/***************************** Function: display_msg *****************************
* Prints an message to stderr.
*
* In:	error_msg	- error message to print
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     display_msg(char *_msg)
#else
void display_msg(_msg)
char	*_msg;
#endif

{
	fprintf(stderr,"%s\n",_msg);
}


/************************ Function: setup_global_coords **********************
* Defines the global coordinate system, making it the first member in the
* linked list of coordinate systems (first_csys_E).
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     setup_global_coords(void)
#else
void setup_global_coords()
#endif

{

	int		i;

	/* Set first coord system to global coordinates
	-----------------------------------------------*/
	first_csys_E = (csys_t *) calloc((size_t) 1,
		sizeof(csys_t));
	if (!first_csys_E)
		p_error("Cannot allocate memory (calloc) for global coord sys",
			NULL);

	first_csys_E->name = (char *) malloc((size_t) 
		strlen((GLOBAL_NAME)+1));
	if (!first_csys_E->name)
		p_error( "Cannot allocate memory for global coord system name",
			NULL);
	strcpy(first_csys_E->name,GLOBAL_NAME);

	for (i=0; i < 3; i++) {
		first_csys_E->origin[i] = 0;
	}
	initialize_matrix(first_csys_E->local_rot,0.0);
	for (i=0; i < 3; i++) {
		first_csys_E->local_rot[i][i] = 1.0;
	}

}


/*************************** Function: open_files ****************************
* Opens the input and output files named by the external variables infile_E
* and outfile_E
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      open_files(void)
#else
int open_files()
#endif

{
	char error_msg[MAX_ERROR_MSG];

	if (infile_E[0] != '\0') {
		if ((ifp_E = fopen(infile_E,"r")) == NULL) {
			sprintf(error_msg,"Cannot open the input file %s",infile_E);
			p_error(error_msg,NULL);
		} /*if*/
	}

	if (outfile_E[0] != '\0') {
		if ((ofp_E = fopen(outfile_E,"w")) == NULL) {
			sprintf(error_msg,"Cannot open the output file %s",outfile_E);
			p_error(error_msg,NULL);
		} /*if*/
	} /*if*/

	return(0);
}


/********************** Function: parse_command_line_args ********************
* Parses the command line argument using getoptPoly3D().
*
* In:	argc	- number of command line arguments
*		argv	- the command line arguments
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      parse_command_line_args(int argc, char *argv[])
#else
int parse_command_line_args(argc, argv)
int		argc;
char	*argv[];
#endif

{
	int		optarg;						/* command line option/argument		*/
	int		exit = FALSE;

	infile_E[0] = outfile_E[0] = '\0';
	while ((optarg = getoptPoly3D("i:o:",argc,argv)) != NO_MORE_ARGS && !exit) {
		switch (optarg) {
			case NO_SUCH_ARG:
			case FILE_ARG:
				exit = TRUE;
				break;
			/* -i <filename> names the input file */
			case 'i':
				if (strlen(getopt_arg_E) > MAXFILE-1) {
					p_error("Input file name too long",NULL);
				}
				strcpy(infile_E,getopt_arg_E);
				break;
			/* -o <filename> names the output file */
			case 'o':
				if (strlen(getopt_arg_E) > MAXFILE-1) {
					p_error("Output file name too long",NULL);
				}
				strcpy(outfile_E,getopt_arg_E);
				break;
		} /*switch*/
	} /*while*/

	if (exit) {
		fprintf(stderr,"\nUsage: poly [-i infile] [-o outfile]");
		fprintf(stderr,"\n       (arguments may occur in any order)\n\n");
		return(ERROR);
	} /*if*/

	return(0);
}


/*************************** Function: read_constants ************************
* Reads problem constants from the input file, skipping blank and comment
* lines.  Stops reading when a line beginning with END_STMT is reached.
* Sets up remote stress and strain tensors and calls calc_elas_consts() to
* calculate undefined elastic constants.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      read_constants(void)
#else
int read_constants()
#endif

{
	int		numwords;
	char	*word[MAXWORDS];
	char	line[MAXLINE];
	char	error_msg[MAX_ERROR_MSG];
	int		temp;
	double	s11r, s22r, s33r, s12r, s13r, s23r;

	s11r = s22r = s33r = s12r = s13r = s23r = 0.0;

	/* read list of constants
	-------------------------*/ 
	for (;;) {

		/* read line from input file, increment linenum
		-----------------------------------------------*/
		numwords = read_line(line,word);
		linenum_E++;

		/* skip blank and comment lines
		-------------------------------*/
		if (numwords == 0) 
			continue;

		/* exit loop when end of list reached
		-------------------------------------*/ 
		if (!strcmp(word[0],END_STMT))
			break;

		/* parse constants
		------------------*/

		if (get_text_var(&title1_E,"title1",word,numwords,line))
			continue;

		if (get_text_var(&title2_E,"title2",word,numwords,line))
			continue;

		if (get_double_var(&shear_mod_E,"shear_mod",word,numwords))
			continue;

		if (get_double_var(&psn_ratio_E,"psn_ratio",word,numwords))
			continue;

		if (get_double_var(&youngs_mod_E,"youngs_mod",word,numwords))
			continue;

		if (get_double_var(&bulk_mod_E,"bulk_mod",word,numwords))
			continue;

		if (get_double_var(&lame_lambda_E,"lame_lambda",word,numwords))
			continue;

		if (get_double_var(&null_value_E,"null_value",word,numwords))
			continue;
/*********************************************** ADDED 98-12-09 */
      if (get_double_var(&coef_exclu_E,"coef_exclu",word,numwords))
			continue;
/***************************************************************/

		if (get_boolean_var(&rem_stress_bc_E,"rem_bc_type",
			"stress","strain",word,numwords,line))
			continue;

		if (get_double_var(&s11r,"s11r",word,numwords))
			continue;

		if (get_double_var(&s22r,"s22r",word,numwords))
			continue;

		if (get_double_var(&s33r,"s33r",word,numwords))
			continue;

		if (get_double_var(&s12r,"s12r",word,numwords))
			continue;

		if (get_double_var(&s13r,"s13r",word,numwords))
			continue;

		if (get_double_var(&s23r,"s23r",word,numwords))
			continue;

		if (get_boolean_var(&half_space_E,"half_space",
			"yes","no",word,numwords,line))
			continue;

		if (get_boolean_var(&check_cond_num_E,"check_cond_num",
			"yes","no",word,numwords,line))
			continue;

		if (get_boolean_var(&print_elt_geom_E,"print_elt_geom",
			"yes","no",word,numwords,line))
			continue;

		if (get_text_var(&elt_geom_csys_name_E,"elt_geom_csys"
			,word,numwords,line))
			continue;

		/* Print error message if line has incorrect format
		---------------------------------------------------*/
		else {
			p_error("Unknown constant, or incorrect format",
				line);
		}
	}

	/* Calculate elastic constants
	------------------------------*/
	if ((temp = calc_elas_consts(&shear_mod_E, &psn_ratio_E, &youngs_mod_E,
		&bulk_mod_E, &lame_lambda_E)) != 0) {
		sprintf(error_msg,"Too %s elastic constants defined (define two)",
			(temp == EC_TOO_FEW) ? "few" : "many");
		p_error(error_msg,NULL);
	} 

	/* Set up remote boundary condition stress and strain tensors
	-------------------------------------------------------------*/
	rem_stress_E[0][0] = s11r;
	rem_stress_E[1][1] = s22r;
	rem_stress_E[2][2] = s33r;
	rem_stress_E[0][1] = rem_stress_E[1][0] = s12r;
	rem_stress_E[0][2] = rem_stress_E[2][0] = s13r;
	rem_stress_E[1][2] = rem_stress_E[1][2] = s23r;
	copy_matrix(rem_stress_E,rem_strain_E);
	if (rem_stress_bc_E) {
		stress_to_strain(rem_strain_E,youngs_mod_E,psn_ratio_E,
			rem_strain_E);
	} else {
		strain_to_stress(rem_stress_E,shear_mod_E,lame_lambda_E,
			rem_stress_E);
	}

	return(0);
}


/************************ Function: read_csystems ***********************
* Reads user coordinate systems from the input file, skipping blank and
* comment lines.  Stops reading when a line beginning with END_STMT is
* reached.  Adds coordinate systems to the linked list whose first member
* (global coords) is given by first_csys_E.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      read_csystems(void)
#else
int read_csystems()
#endif

{
	int		numwords;
	char	*word[MAXWORDS];
	char	error_msg[MAX_ERROR_MSG];
	char	line[MAXLINE];
	char	temp_char;
	int		i,j,k;
	double	rot_matrix[3][3][3];
	double	local_rot[3][3];
	int		rot_order[3];
	double	rot[3];
	csys_t	*parent;
	csys_t	*current_csys;
	

	/* Initialize coordinate rotation matrices
	------------------------------------------*/
	initialize_matrix(local_rot,0.0);
	for (i=0; i < 3; i++) {
		initialize_matrix(rot_matrix[i],0.0);
	}

	/* Set first coordinate system to global coordinates
	----------------------------------------------------*/
	setup_global_coords();
	current_csys = first_csys_E;
			
	/* Read in coordinate systems
	-----------------------------*/
	for (;;) {

		/* read line from input file, increment linenum
		-----------------------------------------------*/
		numwords = read_line(line,word);
		linenum_E++;

		/* Skip blank and comment lines
		-------------------------------*/
		if (numwords == 0)
			continue;

		/* Exit loop when end of list reached
		-------------------------------------*/ 
		if (!strcmp(word[0],END_STMT))
			break;

		/* Check for proper number of parameters
		----------------------------------------*/
		else if (numwords != CS_NUM_PARAMS) {
			sprintf(error_msg,
				"Too %s parameters specified to define coord system",
				(numwords < CS_NUM_PARAMS) ? "few" : "many");
			p_error(error_msg,line);
		}

		/* Check if coordinate system name already taken
		------------------------------------------------*/
		if (find_csys(word[CS_NAME_POS]) != NULL) {
			p_error("Coordinate system already exits",line);
		}

		/* Allocate memory in linked list for local coordinate system
		-------------------------------------------------------------*/
		current_csys->next = (csys_t *)
			calloc((size_t) 1,sizeof(csys_t));
		if (!current_csys->next)
			p_error("Cannot allocate memory (calloc) for local coord sys",
				line);
		current_csys = current_csys->next;

		/* Assign coordinate system name
		--------------------------------*/
		current_csys->name = (char *) malloc((size_t) 
			strlen(word[CS_NAME_POS])+1);
		if (!current_csys->name)
			p_error("Cannot allocate memory for coord system name",
			 line);
		strcpy(current_csys->name,word[CS_NAME_POS]);
		
		/* Get coordinate system parent
		-------------------------------*/
		if ((parent = 
			find_csys(word[CS_PARENT_POS])) == NULL) {
			p_error("Undefined coordinate system",line);
		}

		/* Get coordinate system origin and convert to global coords
		------------------------------------------------------------*/
		for (i=0; i < 3; i++) {
			current_csys->origin[i] = atof(word[CS_ORIGIN_POS+i]);
		}
		transform_position_vector(INVERSE_TRANS,
			parent->origin,parent->local_rot,current_csys->origin);

		/* Get rots about x1,x2,x3 axes of parent and convert to radians
		----------------------------------------------------------------*/
		for (i=0; i < 3; i++) {
			rot[i] = RADIANS(atof(word[CS_ROT_POS+i]));
		}

		/* Get the rotation order
		-------------------------*/
		for (i=0; i < 3; i++) {
			temp_char = word[CS_ROT_ORDER_POS][i];
			if (temp_char < '1' || temp_char > '3') {
				p_error("Invalid axis for rotation order",
					line);
			}
			rot_order[i] = temp_char - '1';
		}

		/* Set up the rotation matrices
		-------------------------------*/
		for (i=0; i < 3; i++) {
			j = i+1;
			k = i+2;
			if (j > 2) j -= 3;
			if (k > 2) k -= 3;
			rot_matrix[i][i][i] = 1.0;
			rot_matrix[i][j][j] =
				rot_matrix[i][k][k] = cos(rot[i]);
			rot_matrix[i][k][j] = 
				-(rot_matrix[i][j][k] =  sin(rot[i]));
		}

		/* Calculate the global-to-local coordinate rotation matrix
		-----------------------------------------------------------*/
		matrix_mult(rot_matrix[rot_order[0]],parent->local_rot,local_rot);
		matrix_mult(rot_matrix[rot_order[1]],local_rot,local_rot);
		matrix_mult(rot_matrix[rot_order[2]],local_rot,local_rot);
		copy_matrix(local_rot,current_csys->local_rot);

	}

	return(0);
}


/********************** Function: read_observation_grids ********************
* Reads observation grids from the input file, skipping blank and
* comment lines.  Stops reading when a line beginning with END_STMT is
* reached.  Adds observation grids to the linked list whose first member
* (global coords) is given by first_obs_grid_E.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      read_observation_grids(void)
#else
int read_observation_grids()
#endif

{
	int		numwords;
	char	*word[MAXWORDS];
	char	line[MAXLINE];
	int		i;
	int		dimension;
	char	error_msg[MAX_ERROR_MSG];
	int		correct_num_params;
	obs_grid_t	 *current_obs_grid;
	int		num_ones;
	int		numpts;
	char	temp_char;

	/* Read in observation grids
	----------------------------*/
	for (;;) {

		/* read line from input file, increment linenum
		-----------------------------------------------*/
		numwords = read_line(line,word);
		linenum_E++;

		/* skip blank and comment lines
		-------------------------------*/
		if (numwords == 0)
			continue;

		/* exit loop when end of list reached
		-------------------------------------*/ 
		if (!strcmp(word[0],END_STMT))
			break;

		/* Get grid dimension
		---------------------*/
		dimension = atoi(word[OG_DIMEN_POS]);

		/* Check for proper number of parameters
		----------------------------------------*/
		switch (dimension) {
			case 0:
				correct_num_params = (numwords == OG_MIN_NUM_PARAMS);
				break;
			case 1:
				correct_num_params = (numwords == (OG_MIN_NUM_PARAMS + 4));
				break;
			case 2:
			case 3:
				correct_num_params = (numwords == (OG_MIN_NUM_PARAMS + 6));
				break;
			default:
				p_error("Invalid dimension for observation grid",
					line);
		}

		if (!correct_num_params) {
			sprintf(error_msg,
				"Incorrect number of parameters for %d-D grid",dimension);
			p_error(error_msg,line);
		}

		/* Allocate memory for observation grid
		---------------------------------------*/
		if (first_obs_grid_E == NULL) {
			first_obs_grid_E = (obs_grid_t *)
				calloc((size_t) 1,sizeof(obs_grid_t));
			if (!first_obs_grid_E)
				p_error("Cannot allocate memory (calloc) for obs grid",
					line);
			current_obs_grid = first_obs_grid_E;
		} else {
			current_obs_grid->next = (obs_grid_t *)
				calloc((size_t) 1,sizeof(obs_grid_t));
			if (!current_obs_grid->next)
				p_error("Cannot allocate memory (calloc) for obs grid",
					line);
			current_obs_grid = current_obs_grid->next;
		}

		/* Set the observation grid dimension
		-------------------------------------*/
		current_obs_grid->dimension = dimension;

		/* Get observation grid name
		----------------------------*/
		current_obs_grid->name = (char *) malloc((size_t) 
			strlen(word[OG_NAME_POS])+1);
		if (!current_obs_grid->name)
			p_error("Cannot allocate memory for observation grid name",
			 line);
		strcpy(current_obs_grid->name,word[OG_NAME_POS]);

		/* Get the print options
		------------------------*/
		i = 0;
		current_obs_grid->print[DISPL]   = FALSE;
		current_obs_grid->print[STRAIN]  = FALSE;
		current_obs_grid->print[STRESS]  = FALSE;
		current_obs_grid->print[PSTRAIN] = FALSE;
		current_obs_grid->print[PSTRESS] = FALSE;
		while ((temp_char = word[OG_PRINT_OPS_POS][i]) != '\0') {
			switch (temp_char) {
				case DISPL_CHAR:
					current_obs_grid->print[DISPL] = TRUE;
					break;
				case STRAIN_CHAR:
					current_obs_grid->print[STRAIN] = TRUE;
					break;
				case STRESS_CHAR:
					current_obs_grid->print[STRESS] = TRUE;
					break;
				case PRINCIPAL_CHAR:
					i++;
					switch (word[OG_PRINT_OPS_POS][i]) {
						case STRAIN_CHAR:
							current_obs_grid->print[PSTRAIN] = TRUE;
							break;
						case STRESS_CHAR:
							current_obs_grid->print[PSTRESS] = TRUE;
							break;
						default:
							p_error("Invalid observation grid print option",
								line);
					}
					break;
				default:
					p_error("Invalid observation grid print option",line);
			}
			i++;
		}
			

		/* Get the input coordinate system
		----------------------------------*/
		if ((current_obs_grid->endpt_csys = 
			find_csys(word[OG_INPUT_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",line);
		}

		/* Get the observation point coordinate system
		----------------------------------------------*/
		if ((current_obs_grid->obspt_csys = 
			find_csys(word[OG_OBSPT_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",line);
		}

		/* Get the output coodinate system
		-----------------------------------*/
		if ((current_obs_grid->outp_csys = 
			find_csys(word[OG_DATA_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",line);
		}

		/* Get the beginning & ending coordinates
		-----------------------------------------*/
		for (i=0; i < 3; i++) {
			current_obs_grid->begin[i] = atof(word[OG_BEGIN_POS+i]);
			if (dimension != 0) {
				current_obs_grid->end[i] = atof(word[OG_END_POS+i]);
			}
		}

		/* Get number of points along each coord axis in grid
		----------------------------------------------------*/
		num_ones = 0;
		for (i=0; i < ((dimension == 2) ? 3 : dimension); i++) {
			current_obs_grid->numpts[i] = numpts =
				atoi(word[OG_NUMPTS_POS+i]);
			if (numpts < 2) {
				if (dimension == 2 && numpts == 1) {
						num_ones++;
				} else {
					sprintf(error_msg,
						"%dD observation grid axes require %d or more points",
						dimension,((dimension == 2) ? 1 : 2));
					p_error(error_msg,line);
				}
			}
		}
		if (dimension == 2 && num_ones != 1) {
			p_error(
				"1 (& only 1) 2D observation grid axis requires 1 point",
				line);
		}

      /*             98-12-09
		   Do not convert the begin & end pts to global coords !!!!!
         So, no transformation in global coordinate system is required
		*/

      /*
		transform_position_vector(INVERSE_TRANS,
			current_obs_grid->endpt_csys->origin,
			current_obs_grid->endpt_csys->local_rot,
			current_obs_grid->begin);
		transform_position_vector(INVERSE_TRANS,
			current_obs_grid->endpt_csys->origin,
			current_obs_grid->endpt_csys->local_rot,
			current_obs_grid->end);
      */

	}

	return(0);
}


/********************* Function: read_objs_elts_verts *********************
* Reads objects, elements, and vertices from the input file, skipping blank
* and comment lines.  Stops reading when a line beginning with END_STMT is
* reached.  Objects, element and vertices are stored in seperate linked-lists.
* Pointers between these lists are set to indicate which elements belong to
* an object and which vertices belong to an element.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      read_objs_elts_verts(void)
#else
int read_objs_elts_verts()
#endif

{
	int		numwords;
	char	*word[MAXWORDS];
	char	line[MAXLINE];
	vert_t *current_vert;
	obj_t *current_obj;
	elt_t *current_elt;

	/* Read in vertices and elements
	--------------------------------*/
	for (;;) {

		/* read line from input file, increment linenum
		-----------------------------------------------*/
		numwords = read_line(line,word);
		linenum_E++;

		/* skip blank and comment lines
		-------------------------------*/
		if (numwords == 0)
			continue;

		/* exit loop when end of list reached
		-------------------------------------*/ 
		if (!strcmp(word[0],END_STMT))
			break;

		/* read vertex info
		-------------------*/
		if (word[0][0] == V_CHAR) {
			get_vert_info(&current_vert,numwords,word,line);
		}

		/* Read object info
		-------------------*/
		else if (word[0][0] == OBJ_CHAR) {
			get_obj_info(&current_obj,numwords,word,line);
		}

		/* Read element info
		--------------------*/
		else if (word[0][0] == E_CHAR) {
			get_elt_info(&current_elt,current_obj,numwords,
				word,line);
		}

	}

	/* Make sure last object contains element(s)
	--------------------------------------------*/
	if (current_obj->first_elt == NULL)
		p_error("Last object contains no elements",NULL);

	return(0);
}


/*********************** Function: print_obj_titles *************************
* Prints the object name and column titles to the output (& temporary output)
* files.
*
* In:	current_obj	- the object for which titles are to be printed
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_obj_titles(obj_t *current_obj)
#else
void print_obj_titles(current_obj)
obj_t *current_obj;
#endif

{
	fprintf(ofp_E,"\n\n====================================================\n");
	fprintf(ofp_E,    "              OBJECT: %s\n",current_obj->name);
	fprintf(ofp_E,    "ELT CENTER COORD SYS: %s\n",
		current_obj->pos_csys->name);
	fprintf(ofp_E,    "====================================================\n");

	/* Print titles
	---------------*/
	if (current_obj->print[DISPL]) {
		fprintf(tempfp_E[DISPL],OBJ_DISPL_TITLE);
		fprintf(tempfp_E[DISPL],OBJ_LOC_LABELS);
		fprintf(tempfp_E[DISPL],OBJ_DISPL_LABELS);
		fprintf(tempfp_E[DISPL],OBJ_BC_CSYS_LABELS);
		fprintf(tempfp_E[DISPL],OBJ_LOC_UNDLNS);
		fprintf(tempfp_E[DISPL],OBJ_DISPL_UNDLNS);
		fprintf(tempfp_E[DISPL],OBJ_BC_CSYS_UNDLNS);
	}
	if (current_obj->print[STRESS]) {
		fprintf(tempfp_E[STRESS],OBJ_STRESS_TITLE);
		fprintf(tempfp_E[STRESS],OBJ_LOC_LABELS);
		fprintf(tempfp_E[STRESS],OBJ_STRESS_LABELS);
		fprintf(tempfp_E[STRESS],OBJ_BC_CSYS_LABELS);
		fprintf(tempfp_E[STRESS],OBJ_LOC_UNDLNS);
		fprintf(tempfp_E[STRESS],OBJ_STRESS_UNDLNS);
		fprintf(tempfp_E[STRESS],OBJ_BC_CSYS_UNDLNS);
	}
}


/********************** Function: print_obs_grid_titles **********************
* Prints the observation grid name and column titles to the output
* (& temporary output) files.
*
* In:	current_obs_grid	- the obs grid for which titles are to be printed
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_obs_grid_titles(obs_grid_t *current_obs_grid)
#else
void print_obs_grid_titles(current_obs_grid)
obs_grid_t *current_obs_grid;
#endif

{
	fprintf(ofp_E,"\n\n====================================================\n");
	fprintf(ofp_E,   "%d-D OBSERVATION GRID: %s\n",
		current_obs_grid->dimension, current_obs_grid->name);
	fprintf(ofp_E,    " OBS POINT COORD SYS: %s\n",
		current_obs_grid->obspt_csys->name);
	fprintf(ofp_E,    "    OUTPUT COORD SYS: %s\n",
		current_obs_grid->outp_csys->name);
	fprintf(ofp_E,    "====================================================\n");

	/* Print titles to temporary files
	----------------------------------*/
	if (current_obs_grid->print[DISPL]) {
		fprintf(tempfp_E[DISPL],OG_DISPL_TITLE);
		fprintf(tempfp_E[DISPL],OG_LOC_LABELS);
		fprintf(tempfp_E[DISPL],OG_DISPL_LABELS);
		fprintf(tempfp_E[DISPL],OG_LOC_UNDLNS);
		fprintf(tempfp_E[DISPL],OG_DISPL_UNDLNS);
	}
	if (current_obs_grid->print[STRAIN]) {
		fprintf(tempfp_E[STRAIN],OG_STRAIN_TITLE);
		fprintf(tempfp_E[STRAIN],OG_LOC_LABELS);
		fprintf(tempfp_E[STRAIN],OG_STRAIN_LABELS);
		fprintf(tempfp_E[STRAIN],OG_LOC_UNDLNS);
		fprintf(tempfp_E[STRAIN],OG_STRAIN_UNDLNS);
	}
	if (current_obs_grid->print[PSTRAIN]) {
		fprintf(tempfp_E[PSTRAIN],OG_PSTRAIN_TITLE);
		fprintf(tempfp_E[PSTRAIN],OG_LOC_LABELS);
		fprintf(tempfp_E[PSTRAIN],OG_PSTRAIN_LABELS);
		fprintf(tempfp_E[PSTRAIN],OG_LOC_UNDLNS);
		fprintf(tempfp_E[PSTRAIN],OG_PSTRAIN_UNDLNS);
	}
	if (current_obs_grid->print[STRESS]) {
		fprintf(tempfp_E[STRESS],OG_STRESS_TITLE);
		fprintf(tempfp_E[STRESS],OG_LOC_LABELS);
		fprintf(tempfp_E[STRESS],OG_STRESS_LABELS);
		fprintf(tempfp_E[STRESS],OG_LOC_UNDLNS);
		fprintf(tempfp_E[STRESS],OG_STRESS_UNDLNS);
	}
	if (current_obs_grid->print[PSTRESS]) {
		fprintf(tempfp_E[PSTRESS],OG_PSTRESS_TITLE);
		fprintf(tempfp_E[PSTRESS],OG_LOC_LABELS);
		fprintf(tempfp_E[PSTRESS],OG_PSTRESS_LABELS);
		fprintf(tempfp_E[PSTRESS],OG_LOC_UNDLNS);
		fprintf(tempfp_E[PSTRESS],OG_PSTRESS_UNDLNS);
	}

}


/************************** Function: print_elt_data ***********************
* Calculates and prints displacement and traction data for an element, given:
*
* In:	current_obj	- the object to which this element belongs
*		current_elt	- the element for which data is to be printed
*		elt_num		- the element number within the object
* Out:	stress		- the stress tensor at the element center due to ALL elts
*		displ		- the displacement of the element center due to ALL elts
*					  EXCEPT this element
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_elt_data(obj_t *current_obj, elt_t *current_elt, int elt_num,double displ[3], double stress[3][3])
#else
void print_elt_data(current_obj, current_elt, elt_num, 
	displ, stress)
obj_t	*current_obj;
elt_t	*current_elt;
int		elt_num;
double	displ[3];
double	stress[3][3];
#endif

{
	double	half_b[3];
	double	displ_pos[3];
	double	displ_neg[3];
	double	normal_vector[3];
	double	traction[3];
	double	b[3];
	double	center[3];
	int		i;


	copy_vector(current_elt->elt_csys.origin,center);
	transform_position_vector(FORWARD_TRANS,
		current_obj->pos_csys->origin,
		current_obj->pos_csys->local_rot,
		center);

	/* Print element location info
	------------------------------*/
	for (i=0; i < NUM_PR_OPTS; i++) {
		if (current_obj->print[i]) {
			fprintf(tempfp_E[i],OBJ_LOC_FMT,elt_num,center[0],
				center[1],center[2]);
		}
	}

	/* Print displacement data if requested
	---------------------------------------*/
	if (current_obj->print[DISPL]) {

		/* Rotate displacement vector to the bc coordinate sys
		------------------------------------------------------*/
		rotate_vector(FORWARD_ROT,current_elt->bc_csys->local_rot,displ);

		/* Cempute the absolute displacements of the pos and neg sides
                   of the element by adding the displacement discontinuity to
	   	   the calculated displacment
		---------------------------------------------------------------*/
		for (i=0; i < 3; i++) {
			b[i] = *current_elt->b[i];
		}
		scalar_vector_mult(0.5,b,half_b);
		add_vectors(displ,half_b,displ_pos);
		subtract_vectors(displ,half_b,displ_neg);

		fprintf(tempfp_E[DISPL],OBJ_DISPL_FMT,
			b[0],displ_pos[0],displ_neg[0],
			b[1],displ_pos[1],displ_neg[1],
			b[2],displ_pos[2],displ_neg[2]);
	}

	/* Print stress (traction) data if requested
	--------------------------------------------*/
	if (current_obj->print[STRESS]) {

		/* Calculate the traction vector on the element plane and
	   	rotate to bc coordinates
		---------------------------------------------------------*/
		normal_vector[0] = normal_vector[1] = 0.0;
		normal_vector[2] = -1.0;
		rotate_vector(INVERSE_ROT,current_elt->elt_csys.local_rot,
			normal_vector);
		cauchy(stress,normal_vector,traction);
		rotate_vector(FORWARD_ROT,
			current_elt->bc_csys->local_rot,traction);


		fprintf(tempfp_E[STRESS],OBJ_STRESS_FMT,traction[0],
			traction[1],traction[2]);
	}

	/* Print BC coord sys name
	--------------------------*/
	for (i=0; i < NUM_PR_OPTS; i++) {
		if (current_obj->print[i]) {
			fprintf(tempfp_E[i],OBJ_BC_CSYS_FMT,
				current_elt->bc_csys->name);
		}
	}

}


/************************ Function: print_obs_pt_data ***********************
* Prints stress, strain, and displacement data for an observation point.
*
* In:	current_obs_grid	- the obs grid to which the obs point belongs
*		x					- coordinates (global) of the observation point
*		displ				- the displacement (global coords) at the obs pt
*		strain				- the strain (global coords) at the obs point
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_obs_pt_data(obs_grid_t *current_obs_grid, double x[3],double displ[3], double strain[3][3])
#else
void print_obs_pt_data(current_obs_grid, x, displ, strain)
obs_grid_t	 *current_obs_grid;
double	x[3];
double	displ[3];
double	strain[3][3];
#endif

{
	double	stress[3][3];
	double	prin[3];
	double	traj[3][3];
	double x_copy[3];
	csys_t *obspt_csys;
	csys_t *outp_csys;
	int		i;

	obspt_csys = current_obs_grid->obspt_csys;
	outp_csys = current_obs_grid->outp_csys;

	copy_vector(x,x_copy);
	transform_position_vector(FORWARD_TRANS,obspt_csys->origin,
		obspt_csys->local_rot,x_copy);

	/* Print observation point location to temp files
	-------------------------------------------------*/
	for (i=0; i < NUM_PR_OPTS; i++) {
		if (current_obs_grid->print[i]) {
			fprintf(tempfp_E[i],OG_LOC_FMT,x_copy[0],x_copy[1],x_copy[2]);
		}
	}

	/* Print displacement data to temp file
	---------------------------------------*/
	if (current_obs_grid->print[DISPL]) {
		rotate_vector(FORWARD_ROT,outp_csys->local_rot,displ);
		if (below_vertex_E || near_vertex_E)
			initialize_vector(displ,null_value_E);
		fprintf(tempfp_E[DISPL],OG_DISPL_FMT,displ[0],displ[1],displ[2]);
	}

	/* Print stress data to temp file
	---------------------------------*/
	if (current_obs_grid->print[STRESS] || current_obs_grid->print[PSTRESS]) {

		strain_to_stress(strain, shear_mod_E, lame_lambda_E, stress);
		rotate_tensor(FORWARD_ROT,outp_csys->local_rot,stress);

		if (current_obs_grid->print[STRESS]) {
			if (below_vertex_E || near_vertex_E)
				initialize_matrix(stress,null_value_E);
			fprintf(tempfp_E[STRESS],OG_STRESS_FMT,stress[0][0],stress[1][1],
				stress[2][2], stress[0][1], stress[1][2], stress[0][2]);
		}

		if (current_obs_grid->print[PSTRESS]) {
			/* Calculate principal stresses
			-------------------------------*/
			if (below_vertex_E || near_vertex_E) {
				initialize_vector(prin,null_value_E);
				initialize_matrix(traj,null_value_E);
			} else {
				principal(stress, prin, traj);
			}
			fprintf(tempfp_E[PSTRESS],OG_PSTRESS_FMT,
				traj[0][0], traj[0][1], traj[0][2], prin[0],
				traj[1][0], traj[1][1], traj[1][2], prin[1],
				traj[2][0], traj[2][1], traj[2][2], prin[2]);
		}
	}


	/* Print strain data to temp file
	---------------------------------*/
	if (current_obs_grid->print[STRAIN] || current_obs_grid->print[PSTRAIN]) {

		rotate_tensor(FORWARD_ROT,outp_csys->local_rot,strain);

		if (current_obs_grid->print[STRAIN]) {
			if (below_vertex_E || near_vertex_E)
				initialize_matrix(strain,null_value_E);
			fprintf(tempfp_E[STRAIN],OG_STRAIN_FMT,strain[0][0],strain[1][1],
				strain[2][2], strain[0][1], strain[1][2], strain[0][2]);
		}

		if (current_obs_grid->print[PSTRAIN]) {
			/* Calculate principal strains
			------------------------------*/
			if (below_vertex_E || near_vertex_E) {
				initialize_vector(prin,null_value_E);
				initialize_matrix(traj,null_value_E);
			} else {
				principal(strain, prin, traj);
			}
			fprintf(tempfp_E[PSTRAIN],OG_PSTRAIN_FMT,
				traj[0][0], traj[0][1], traj[0][2], prin[0],
				traj[1][0], traj[1][1], traj[1][2], prin[1],
				traj[2][0], traj[2][1], traj[2][2], prin[2]);
		}
	}

	below_vertex_E = FALSE;

   /* NEW: 98-12-09 */
   near_vertex_E  = FALSE;
}


/******************** Function: determine_burgers_vectors ********************
* Sets up a system of linear equations to solve for unknown Burgers vector
* components.  Each traction boundary condition component given in the input
* file leads to one equation and one unknown.  Solves the system of equations
* using the functions d_ludcmp() and d_lubksb() adapted from the book
* "Numerical Recipes in C, 2nd Ed." (Press et al., 1992).
*
* Modifies:
*	b_vector[]          - static array holding (previously) unknown Burgers
*                         vector components
*   "current_elt"->b[i] - pointer to ith component of Burgers vector.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     determine_burgers_vectors()
#else
void determine_burgers_vectors()
#endif

{
	elt_t	*elt1;
	elt_t	*elt2;
	int		num_eqns = 0;
	int		eqn_num = 0;
	double **ic_matrix;
	double **ic_matrix_inv;
	double *b_vector;
	double *column;
	double	norm1;
	double	norm2;
	int		*pivot;
	int		i, j;
	int		row;
	int		rplus;
	int		col;
	int		cplus;
	int		index;
	double	stress_ic[3][3][3];
	double	traction_ic[3][3];
	double	displ_ic[3][3];
	double	strain_ic[3][3][3];
	double	d;
	double	normal_vector[3];
 
#ifdef DEBUG 
	FILE *outMat, *outVec; 
#endif

	/* Determine number of simultaneous equations
	---------------------------------------------*/
	elt1 = first_elt_E;
	while (elt1 != NULL) {
		for (i=0; i < 3; i++) {
			if (elt1->bc_type[i] == TRACTION_BC)
				num_eqns++;
		}
		elt1 = elt1->next;
	}

	/* Allocate memory for matrix and vectors. Note that indicies range
	   from [1..num_eqns] rather than [0..num_eqns-1] for compatibility
	   with "Numerical Recipes" functions used below
	-------------------------------------------------------------------*/
	ic_matrix		= dmatrix(1,num_eqns,1,num_eqns);
	ic_matrix_inv	= dmatrix(1,num_eqns,1,num_eqns);
	column			= dvector(1,num_eqns);
	b_vector		= dvector(1,num_eqns);
	pivot			= ivector(1,num_eqns); 
 
	printf( "\n\n  -->nbr eq= %i", num_eqns );

	/* Set up vectors
	-----------------*/
	eqn_num = 1;
	elt1 = first_elt_E;
	while (elt1 != NULL) {
		for (i=0; i < 3; i++) {
			switch (elt1->bc_type[i]) {
				case BVECTOR_BC:
					elt1->b[i] = &elt1->bc[i];
					break;
				case TRACTION_BC:
					elt1->b[i] = &b_vector[eqn_num];
					b_vector[eqn_num] = elt1->bc[i];
					eqn_num++;
					break;
			}
		}
		elt1 = elt1->next;
	}
 
 
#ifdef DEBUG 
	outMat = fopen("out_matrix.txt","w"); 
	outVec = fopen("out_vector.txt","w"); 
	for (i=1; i<=num_eqns; i++) 
		fprintf(outVec,"%f\n",b_vector[i]); 
	fprintf(outVec,"\n"); 
#endif 
 

	row		= 1;
	index	= 0;
	elt1 = first_elt_E;
	while (elt1 != NULL) {
		index++;
		if (elt1->bc_type[0] == TRACTION_BC || 
			elt1->bc_type[1] == TRACTION_BC ||
			elt1->bc_type[2] == TRACTION_BC) {

			/* Determine the element normal
			-------------------------------*/
			normal_vector[0] = normal_vector[1] = 0.0;
			normal_vector[2] = -1.0;
			rotate_vector(INVERSE_ROT,elt1->elt_csys.local_rot, normal_vector); 

			/* Go to first column of influence coefficient matrix
			-----------------------------------------------------*/
			col = 1;

			/* Loop over elements
			---------------------*/
			elt2 = first_elt_E;
			while (elt2 != NULL) {

				/* Calculate the stress inf coeffs due to this element
				------------------------------------------------------*/
				displ_strain_ics_poly_elt(FALSE,TRUE,elt2,elt1->elt_csys.origin, 
					displ_ic, strain_ic, elt2);
				if (below_vertex_E)
					p_error(
						"One elt's center lies directly under another's vertex",
						NULL);
				for (i=0; i < 3; i++) {
					strain_to_stress(strain_ic[i],shear_mod_E,lame_lambda_E,
						stress_ic[i]);
				}

				/* Resolve them into traction inf coeffs
				----------------------------------------*/
				for (i=0; i < 3; i++) {
					cauchy(stress_ic[i],normal_vector,traction_ic[i]);
					rotate_vector(FORWARD_ROT,elt1->bc_csys->local_rot,
						traction_ic[i]);
				}

				/* Set up inf coeff matrix terms due to this element
				----------------------------------------------------*/
				rplus = 0;
				for (i=0; i < 3; i++) {
					if (elt1->bc_type[i] == TRACTION_BC) {
						cplus = 0;
						for (j=0; j < 3; j++) {
							switch (elt2->bc_type[j]) {
								case BVECTOR_BC:
									b_vector[row+rplus] -= traction_ic[j][i]*
										(*elt2->b[j]);
									break;
								case TRACTION_BC:
									ic_matrix[row+rplus][col+cplus] =
										traction_ic[j][i];
									cplus++;
									break;
							} /*switch*/
						} /*for*/
						rplus++;
					} /*if*/
				} /*for*/
				col += cplus;
				elt2 = elt2->next;
			} /*while*/

			row += rplus;
		} /*if*/

		elt1 = elt1->next;
	} /*while*/

	/* Calculate L-infinity maximum matrix of ic_matrix 
	---------------------------------------------------*/
	if (check_cond_num_E && num_eqns != 0) {
		norm1 = array_max_norm(ic_matrix,1,num_eqns,1,num_eqns);
	}

	/* Solve the system of linear eqns (courtesy of "Numerical Recipes")
	--------------------------------------------------------------------*/
#ifdef DEBUG 
	for (i=1; i<=num_eqns; i++){ 
		for (j=1; j<=num_eqns; j++) 
			fprintf(outMat,"%f\t",ic_matrix[j][i]); 
		fprintf(outMat,"\n"); 
	} 
	for (i=1; i<=num_eqns; i++) 
		fprintf(outVec,"%f\n",b_vector[i]); 
#endif 

        time_t start, time_lu, time_bs;
        time(&start);
	
        d_ludcmp(ic_matrix, num_eqns, pivot, &d); 
        time(&time_lu);
        time_t elapsed = difftime(time_lu, start);
        printf("\nLU in %3.3f s\n", elapsed);
	
        d_lubksb(ic_matrix, num_eqns, pivot, b_vector); 
        time(&time_bs);
        elapsed = difftime(time_bs, time_lu);
        printf("\nBacksolve in %3.3f s\n", elapsed);


#ifdef DEBUG 
	fprintf(outMat,"\n");fprintf(outMat,"\n"); 
	fprintf(outVec,"\n");fprintf(outVec,"\n"); 
	for (i=1; i<=num_eqns; i++){ 
		for (j=1; j<=num_eqns; j++) 
			fprintf(outMat,"%f\t",ic_matrix[j][i]); 
		fprintf(outMat,"\n"); 
	} 
	for (i=1; i<=num_eqns; i++) 
		fprintf(outVec,"%f\n",b_vector[i]); 
#endif 


	/* Compute the matrix condition number if requested. 
	----------------------------------------------------*/
	if (check_cond_num_E && num_eqns != 0) {

		/* Finish inverting ic_matrix
		-----------------------------*/
		for (i=1; i <= num_eqns; i++) {
			for (j=1; j <= num_eqns; j++)
				column[j] = 0.0;
			column[i] = 1.0;
			d_lubksb(ic_matrix, num_eqns, pivot, column);
			for (j=1; j <= num_eqns; j++)
				ic_matrix_inv[j][i] = column[j];
		}

		/* Calculate L-infinity maximum norm of the inverted ic_matrix
		--------------------------------------------------------------*/
		norm2 = array_max_norm(ic_matrix_inv,1,num_eqns,1,num_eqns);

		cond_num_E = norm1 * norm2; 
 
#ifdef DEBUG 
	fprintf(outMat,"\n\nCondition number");fprintf(outMat,"\n"); 
	fprintf(outVec,"\n\nCondition number");fprintf(outVec,"\n"); 
	for (i=1; i<=num_eqns; i++){ 
		for (j=1; j<=num_eqns; j++) 
			fprintf(outMat,"%f\t",ic_matrix_inv[j][i]); 
		fprintf(outMat,"\n"); 
	} 
	fclose(outVec); 
	fclose(outMat); 
#endif 

	}

	/* Free memory used by vectors/matrices, EXCEPT for b_vector
	------------------------------------------------------------*/
	free_dmatrix(ic_matrix,1,num_eqns,1,num_eqns);
	free_dmatrix(ic_matrix_inv,1,num_eqns,1,num_eqns);
	free_dvector(column,1,num_eqns);
	free_ivector(pivot,1,num_eqns);
}


/*************************** Function: read_line ***************************
* Uses getwords() to read a line from the input file. Adjusts numwords
* appropriately if the line contains a comment character.  Returns number
* of words on line.
*
* Out:	line	- the line read from the input file
* In:	word[i]	- the ith word on the line (null-terminated string)
***************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      read_line(char *line, char *word[])
#else
int read_line(line, word)
char	*line;
char	*word[];
#endif

{
	int		numwords;
	/*char	error_msg[MAX_ERROR_MSG];*/
	int		i, j;
	int		exit;

	/* get line, exit function on EOF
	---------------------------------*/
	if ((numwords = getwords(ifp_E,line,MAXLINE,word,MAXWORDS,
		CONTINUE_CHAR)) < 0) {
		switch (numwords) {
			case GW_EOF_ERR:
				p_error("Unexpected EOF in getwords()\n",NULL);
			case GW_MALLOC_ERR:
				p_error("Memory allocation error in getwords()\n",NULL);
			case GW_MAXWORDS_ERR:
				p_error("Too many words error in getwords()\n",NULL);
		}
	} /*if*/

	/* if line contains a comment, adjust numwords accordingly
	----------------------------------------------------------*/
	exit = FALSE;
	for (i=0; !exit && (i < numwords); i++) {

		/* Loop over characters in word.  If comment character is found,
		   replace with '\0' and exit loop.
		----------------------------------------------------------------*/
		j = 0;
		while (word[i][j] != '\0') {
			if (word[i][j] == COMMENT_CHAR) {
				word[i][j] = '\0';
				exit = TRUE;
				break;
			}
			j++;
		}
	}
	if (exit) 
		numwords = (j == 0) ? i-1 : i;

	return(numwords);
}


/************************* Function: get_double_var **************************
* Function for assigning values to double variables.  Returns TRUE if an
* assignment is made, or FALSE otherwise.
*
* In:		var_name	- name used for this variable in the input file
*			word		- array of words given in the input line
*			numwords	- number of words on the input line
* In/Out:	var			- variable to which a value is to be assigned
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      get_double_var(double *var, char *var_name, char *word[], int numwords)
#else
int get_double_var(var, var_name, word, numwords)
double	*var;
char	*var_name;
char	*word[];
int		numwords;
#endif

{
	if (!strcmp(word[CONST_NAME_POS],var_name)) {
		if (numwords == CONST_NUM_PARAMS) {
			*var = atof(word[CONST_VALUE_POS]);
		}
		return(TRUE);
	}
	return(FALSE);
}


/************************* Function: get_boolean_var *************************
* Function for assigning values to boolean variables.  Returns TRUE if an
* assignment is made, or FALSE otherwise.
*
* In:		var_name		- name used for this variable in the input file
*			true_string		- true value string for this var in input file
*			false_string	- false value string for this var in input file
*			word			- array of words given in the input line
*			numwords		- number of words on the input line
*			line			- line read from the input file
* In/Out:	var				- variable to which a value is to be assigned
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      get_boolean_var(int *var, char *var_name, char *true_string,char *false_string, char *word[], int	numwords, char *line)
#else
int get_boolean_var(var, var_name, true_string, false_string, word, numwords,
	line)
int		*var;
char	*var_name;
char	*true_string;
char	*false_string;
char	*word[];
int		numwords;
char	*line;
#endif

{
	char	error_msg[MAX_ERROR_MSG];

	if (!strcmp(word[CONST_NAME_POS],var_name)) {
		if (numwords == CONST_NUM_PARAMS) {
			if (!strcmp(word[CONST_VALUE_POS],true_string))
				*var = TRUE;
			else if (!strcmp(word[CONST_VALUE_POS],false_string))
				*var = FALSE;
			else {
				sprintf(error_msg, "%s requires a value of \"%s\" or \"%s\"",
					var_name, true_string, false_string);
				p_error(error_msg,line);
			}
		}
		return(TRUE);
	}
	return(FALSE);
}


/************************** Function: get_text_var ***************************
* Function for assigning values to test variables.  Returns TRUE if an
* assignment is made, or FALSE otherwise.
*
* In:		var_name	- name used for this variable in the input file
*			word		- array of words given in the input line
*			numwords	- number of words on the input line
*			line		- line read from the input file
* In/Out:	var			- variable to which a value is to be assigned
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int      get_text_var(char **var, char *var_name, char *word[], int numwords,char *line)
#else
int get_text_var(var, var_name, word, numwords, line)
char	**var;
char	*var_name;
char	*word[];
int		numwords;
char	*line;
#endif

{
	char	error_msg[MAX_ERROR_MSG];
	
	if (!strcmp(word[CONST_NAME_POS],var_name)) {
		if (numwords == CONST_NUM_PARAMS) {
			*var = (char *) malloc((size_t) strlen(word[CONST_VALUE_POS])+1);
			if (!*var) {
				sprintf(error_msg,"Cannot allocate memory for %s",var_name);
				p_error(error_msg,line);
			}
			strcpy(*var,word[CONST_VALUE_POS]);
		}
		return(TRUE);
	}
	return(FALSE);
}


/************************* Function: get_vert_info *************************
* Reads vertex info and sets up new member in linked list of vertices.
*
* In:		numwords		- number of words on input line
*			word			- array of words given in the input line
*			line			- the input line
* In/Out:	current_vert	- the current (last) vertex in the linked list
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     get_vert_info(vert_t **current_vert, int numwords, char *word[],char *line)
#else
void get_vert_info(current_vert, numwords, word, line)
vert_t	**current_vert;
int		numwords;
char	*word[];
char	*line;
#endif

{
	int		i;
	char	error_msg[MAX_ERROR_MSG];

	/* Check for proper number of parameters
	----------------------------------------*/
	if (numwords != V_NUM_PARAMS) {
		sprintf(error_msg,
			"Too %s parameters specified to define vertex",
			(numwords < V_NUM_PARAMS) ? "few" : "many");
		p_error(error_msg,line);
	}

	/* Allocate memory for vertex
	-----------------------------*/
	if (first_vert_E == NULL) {
		first_vert_E = (vert_t *) calloc((size_t) 1,sizeof(vert_t)); 
		if (!first_vert_E)
			p_error("Cannot allocate memory (calloc) for vertex",
				line);
		*current_vert = first_vert_E;
	} else {
		(*current_vert)->next = (vert_t *)
			calloc((size_t) 1,sizeof(vert_t));
		if (!(*current_vert)->next)
			p_error("Cannot allocate memory (calloc) for vertex",
				line);
		*current_vert = (*current_vert)->next;
	}

	/* Get vertex name
	------------------*/
	(*current_vert)->name = (char *) malloc((size_t) 
		strlen(word[V_NAME_POS])+1);
	if (!(*current_vert)->name)
		p_error("Cannot allocate memory for vertex name",
		 line);
	strcpy((*current_vert)->name,word[V_NAME_POS]);

	/* Get the coodinate system
	----------------------------*/
	if (((*current_vert)->csys = find_csys(word[V_CSYS_POS]))
		== NULL) {
		p_error("Undefined coordinate system",
			line);
	}

	/* Read the vertex position
	---------------------------*/
	for (i=0; i < 3; i++)
		(*current_vert)->x[i] = atof(word[V_X_POS+i]);

	/* Transform vertex position vector to global coords
	----------------------------------------------------*/
	transform_position_vector(INVERSE_TRANS,
		(*current_vert)->csys->origin,
		(*current_vert)->csys->local_rot,
		(*current_vert)->x);
}


/************************* Function: get_obj_info ***************************
* Reads object info and sets up new member in linked list of objects.
*
* In:		numwords		- number of words on input line
*			word			- array of words given in the input line
*			line			- the input line
* In/Out:	current_obj		- the current (last) object in the linked list
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     get_obj_info(obj_t **current_obj, int numwords, char *word[], char *line)
#else
void get_obj_info(current_obj, numwords, word, line)
obj_t	**current_obj;
int		numwords;
char	*word[];
char	*line;
#endif

{
	int		i;
	char	temp_char;

	/* Check for proper number of parameters
	----------------------------------------*/
	if (numwords != OBJ_MIN_NUM_PARAMS &&
		numwords != OBJ_MIN_NUM_PARAMS+2) {
		p_error(
			"Incorrect number of parameters specified to define object",
			line);
	}

	/* Allocate memory for object
	-----------------------------*/
	if (first_obj_E == NULL) {
		first_obj_E = (obj_t *) calloc((size_t) 1,sizeof(obj_t));
		if (!first_obj_E)
			p_error("Cannot allocate memory (calloc) for object",
				line);
		*current_obj = first_obj_E;
	} else {
		(*current_obj)->next = (obj_t *)
			calloc((size_t) 1,sizeof(obj_t));
		if (!(*current_obj)->next)
			p_error("Cannot allocate memory (calloc) for object",
				line);
		*current_obj = (*current_obj)->next;
	}

	/* Get object name
	------------------*/
	(*current_obj)->name = (char *) malloc((size_t) 
		strlen(word[OBJ_NAME_POS])+1);
	if (!(*current_obj)->name)
		p_error("Cannot allocate memory for object name",
		line);
	strcpy((*current_obj)->name,word[V_NAME_POS]);

	if (numwords > OBJ_MIN_NUM_PARAMS) {

		/* Get the print options
		------------------------*/
		i = 0;
		(*current_obj)->print[DISPL] = 
			(*current_obj)->print[STRESS] = FALSE;
		while ((temp_char = word[OBJ_PRINT_OPS_POS][i]) != '\0') {
			switch (temp_char) {
				case BVECTOR_CHAR:
					(*current_obj)->print[DISPL] = TRUE;
					break;
				case TRACTION_CHAR:
					(*current_obj)->print[STRESS] = TRUE;
					break;
				default:
					p_error("Invalid object print option",line);
			}
			i++;
		}

		/* Get the coordinate system for element positions
		--------------------------------------------------*/
		if (((*current_obj)->pos_csys = 
			find_csys(word[OBJ_POS_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",line);
		}

	}

	(*current_obj)->first_elt = NULL;
}
	

/************************* Function: get_elt_info ************************
* Reads element info and sets up new member in linked list of elements.
*
* In:		numwords		- number of words on input line
*			word			- array of words given in the input line
*			line			- the input line
* In/Out:	current_elt		- the current (last) element in the linked list
*			current_obj		- the object to which current_elt belongs
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     get_elt_info(elt_t **current_elt, obj_t *current_obj, int numwords,char *word[], char *line)
#else
void get_elt_info(current_elt, current_obj, numwords, word, line)
elt_t	**current_elt;
obj_t	*current_obj;
int		numwords;
char	*word[];
char	*line;
#endif

{
	int		num_vertices;
	int		num_params;
	int		i, j;
	char	error_msg[MAX_ERROR_MSG];
	vert_t	 *vert1;
	vert_t	 *vert2;
	double  dx[3];
	double	trend;
	double	plunge;
	double	rot2[3][3];
	double	rot1[3][3];
	double	trac_bc_adjust[3];
	int		seg;
	vert_t	*vert[3];
	double	vert_x[3];
	double		vector1[3];
	double		vector2[3];
	double	normal_vector[3];
	double		x[3][3];
	double		global_x[3][3];
	disloc_seg_t	*disloc_seg;
	static char	*elt_csys_name;


	/* Increment num_elts_E
	-----------------------*/
	num_elts_E++;

	/* Get the number of vertices
	-----------------------------*/
	num_vertices = atoi(word[E_NUM_VERT_POS]);
	if (num_vertices < 3) 
		p_error("Element must have at least three vertices",line);

	/* Check for proper number of parameters
	----------------------------------------*/
	num_params = E_MIN_NUM_PARAMS+(num_vertices-3);
	if (numwords != num_params) {
		sprintf(error_msg,
			"Too %s parameters specified to define %d-sided element",
			((numwords < num_params) ? "few" : "many"),num_vertices);
		p_error(error_msg,line);
	}

	/* Allocate memory for this element
	-----------------------------------*/
	if (first_elt_E == NULL) {
		first_elt_E = (elt_t *) malloc(sizeof(elt_t));
		if (!first_elt_E)
			p_error("Cannot allocate memory for element",line);
		*current_elt = first_elt_E;
		(*current_elt)->next = NULL;
	} else {
		(*current_elt)->next = (elt_t *) malloc(sizeof(elt_t));
		if (!(*current_elt)->next)
			p_error("Cannot allocate memory for element",line);
		*current_elt = (*current_elt)->next;
		(*current_elt)->next = NULL;
	}

	/* Set object pointers to first and last elements
	-------------------------------------------------*/
	if (first_obj_E == NULL)
		p_error("No objects defined. Element must be part of an object",line);
	if (current_obj->first_elt == NULL)
		current_obj->first_elt = *current_elt;
	current_obj->last_elt = *current_elt;

	/* Set the number of vertices
	-----------------------------*/
	(*current_elt)->num_vertices = num_vertices;

	/* Allocate memory for dislocation segment array
	------------------------------------------------*/
	(*current_elt)->disloc_seg = (disloc_seg_t *)
		calloc((size_t) num_vertices, sizeof(disloc_seg_t));
	if (!(*current_elt)->disloc_seg)
		p_error("Cannot allocate memory for dislocation segment array",
			line);

	/* Get dislocation segment vertices
	-----------------------------------*/
	for (i = 0; i < num_vertices; i++) {
		j = (i == (num_vertices-1)) ? 0 : i+1;
		if (i != 0) {
			(*current_elt)->disloc_seg[i].vert[0] =
			(*current_elt)->disloc_seg[i-1].vert[1];
		} else {
			if (((*current_elt)->disloc_seg[i].vert[0] =
				find_vert(word[E_VERTEX_POS+i])) == NULL) {
				p_error("Undefined vertex",line);
			}
		}
		if (((*current_elt)->disloc_seg[i].vert[1] =
			find_vert(word[E_VERTEX_POS+j])) == NULL) {
			p_error("Undefined vertex",line);
		}
	}


	/* Inititalize rotation matrices
	--------------------------------*/
	initialize_matrix(rot2,0.0);
	initialize_matrix(rot1,0.0);

	initialize_vector((*current_elt)->elt_csys.origin,0.0);

	/* Loop over the dislocation segments
	-------------------------------------*/
	disloc_seg = (*current_elt)->disloc_seg;
	for (seg = 0; seg < (*current_elt)->num_vertices; seg++) {

		/* Calculate this dislocation segment's first vertex's contribution
		   to the element center
		-------------------------------------------------------------------*/
		for (i=0; i < 3; i++) {
			(*current_elt)->elt_csys.origin[i] += 
				disloc_seg[seg].vert[0]->x[i] / (*current_elt)->num_vertices;
		}

		/* Determine the two vertices for this dislocation segment
		----------------------------------------------------------*/
		vert1 = disloc_seg[seg].vert[0];
		vert2 = disloc_seg[seg].vert[1];

		/* Calculate trend and plunge of this dislocation segment
		--------------------------------------------------------*/
		subtract_vectors(vert2->x,vert1->x,dx);
		trend = PI/2.0 - safe_atan2(dx[1],dx[0]);
		disloc_seg[seg].trend = trend;
		plunge = -safe_atan(dx[2],sqrt(dx[0]*dx[0] + dx[1]*dx[1]));
		disloc_seg[seg].plunge = plunge;

		/* Calculate the segment-local (Comninou & Dunders) to global
		   coordinates rotation matrix
		------------------------------------------------------------*/
		rot2[0][0] = 1.0;
		rot2[1][1] = rot2[2][2] = -1.0;
		rot1[0][0] = rot1[1][1] = sin(trend);
		rot1[1][0] = -(rot1[0][1] = cos(trend));
		rot1[2][2] = 1.0;
		matrix_mult(rot2,rot1,disloc_seg[seg].local_rot);

	} /*for*/

	/* Compute element's local coordinate system
	--------------------------------------------*/
	if (elt_csys_name == NULL) {
		elt_csys_name = (char *) malloc((size_t) strlen(ELT_CSYS_NAME)+1);
		if (!elt_csys_name)
			p_error("Cannot allocate memory for elt-local coord sys name",
				NULL);
		strcpy(elt_csys_name,ELT_CSYS_NAME);
	}
	(*current_elt)->elt_csys.name = elt_csys_name;

	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			global_x[i][j] = (i == j) ? 1.0 : 0.0;
		}
	}
	for (i=0; i < 3; i++) {
		vert[i] = disloc_seg[(i*(*current_elt)->num_vertices)/3].vert[0];
	}
	subtract_vectors(vert[1]->x,vert[0]->x,vector1);
	subtract_vectors(vert[2]->x,vert[0]->x,vector2);
	normalize_vector(vector1);
	normalize_vector(vector2);

	cross_product(vector1,vector2,x[2]);
	if (vector_magnitude(x[2]) < TINY_ANGLE)
		p_error(
			"Cannot calc element normal. Elt must have a very odd shape.",
			line);
	normalize_vector(x[2]);
	cross_product(global_x[2],x[2],x[1]);
	if (vector_magnitude(x[1]) < TINY_ANGLE)
		copy_vector(global_x[1],x[1]);
	normalize_vector(x[1]);
	cross_product(x[1],x[2],x[0]);
	normalize_vector(x[0]);

	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			(*current_elt)->elt_csys.local_rot[i][j] =
			dot_product(x[i],global_x[j]);
		}
	}

	/* Check that all vertices are co-planar
	----------------------------------------*/
	for (seg=0; seg < num_vertices; seg++) {
		copy_vector((*current_elt)->disloc_seg[seg].vert[0]->x,vert_x);
		transform_position_vector(FORWARD_TRANS,
			(*current_elt)->elt_csys.origin,
			(*current_elt)->elt_csys.local_rot,vert_x);
		if (fabs(vert_x[2]) >
			fabs(sqrt(vert_x[0]*vert_x[0]+vert_x[1]*vert_x[1])/COPLANAR_LIMIT))
			p_error("Vertices are not co-planar",line);
	}

	/* Get the bc coodinate system
	------------------------------*/
	if (!strcmp(word[E_BC_CSYS_POS],ELT_CSYS_NAME)) {
		(*current_elt)->bc_csys = &((*current_elt)->elt_csys);
	} else if (((*current_elt)->bc_csys =
		find_csys(word[E_BC_CSYS_POS])) == NULL) {
		p_error("Undefined coordinate system", line);
	}

	/* Read the boundary condition types
	------------------------------------*/
	for (i=0; i < 3; i++) {
		if (word[E_BC_TYPE_POS][i] == BVECTOR_CHAR)
			(*current_elt)->bc_type[i] = BVECTOR_BC;
		else if (word[E_BC_TYPE_POS][i] == TRACTION_CHAR)
			(*current_elt)->bc_type[i] = TRACTION_BC;
		else {
			p_error("Invalid boundary condition type",
				line);
		}
	}
	
	/* Calculate adjustment to traction BCs due to remote stresses
	--------------------------------------------------------------*/
	scalar_vector_mult(-1.0,x[2],normal_vector);
	cauchy(rem_stress_E,normal_vector,trac_bc_adjust);
	rotate_vector(FORWARD_ROT,(*current_elt)->bc_csys->local_rot,trac_bc_adjust);

	/* Read the boundary condition values
	-------------------------------------*/
	for (i=0; i < 3; i++) {
		(*current_elt)->bc[i] = atof(word[E_BC_POS+i]);

		/* Adjust traction BC components for remote stress
		--------------------------------------------------*/
		if ((*current_elt)->bc_type[i] == TRACTION_BC)
			(*current_elt)->bc[i] -= trac_bc_adjust[i];
	}

	/* Calculate the projection of each unit component of the element
	   burgers vector into a segment-local coordinates burgers vector
	   (for calculating influence coefficients)
	-----------------------------------------------------------------*/
	for (seg = 0; seg < (*current_elt)->num_vertices; seg++) {
		for (i=0; i < 3; i++) {
			initialize_vector(disloc_seg[seg].elt_b[i],0.0);
			disloc_seg[seg].elt_b[i][i] = -1.0;
			rotate_vector(INVERSE_ROT,
				(*current_elt)->bc_csys->local_rot,
				disloc_seg[seg].elt_b[i]);
			rotate_vector(FORWARD_ROT,disloc_seg[seg].local_rot,
				disloc_seg[seg].elt_b[i]);
		}
	}
}


/************************* Function: open_temp_files *************************
* Opens the temporary output files for displacement, strain, principal strain,
* stress, and principal stress data.
*
* In:	print	- array of flags giving the data to be printed (determines
*				  which temp files need to be opened)
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     open_temp_files(int print[])
#else
void open_temp_files(print)
int	print[];
#endif


{
	int		i;

	/* Open temporary files
	-----------------------*/
	for (i=0; i < NUM_PR_OPTS; i++) {
		tempfp_E[i] = NULL;
		if (print[i]) {
			if ((tempfp_E[i] = tmpfile()) == NULL) {
				p_error("Cannot open temporary file",NULL);
			} /*if*/
		} /*if*/
	} 

}


/************************ Function: close_temp_files ************************
* Closes the temporary output files for displacement, strain, principal strain,
* stress, and principal stress data.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     close_temp_files(void)
#else
void close_temp_files()
#endif


{
	int		i;

	/* Close temporary files
	-----------------------*/
	for (i=0; i < NUM_PR_OPTS; i++) {
		if (tempfp_E[i] != NULL) {
			if (fclose(tempfp_E[i]) == EOF) {
				p_error("Error closing temporary file",NULL);
			} /*if*/
			tempfp_E[i] = NULL;
		} /*if*/
	}

}


/********************** Function: copy_temp_files ***************************
* Copues the temporary output files for displacement, strain, principal strain,
* stress, and principal stress data to the main output file.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     copy_temp_files(void)
#else
void copy_temp_files()
#endif

{
		int		i;
		FILE	*tfp;
		int		c;

		/* Rewind temp files and copy to main output file
		-------------------------------------------------*/
		for (i=0; i < NUM_PR_OPTS; i++) {
			if ((tfp = tempfp_E[i]) != NULL) {
				rewind(tfp);
				while ((c = fgetc(tfp)) != EOF)
					fputc(c,ofp_E);
			}
		}
}


/************************** Function: print_elt_geometry *********************
* Loops over all objects, printing the geometry of each element (name &
* coordinates of each vertex) to the output file. 
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     print_elt_geometry(void)
#else
void print_elt_geometry()
#endif
{
	csys_t	*elt_geom_csys;
	obj_t	*current_obj;
	elt_t	*current_elt;
	int		elt_num;
	int		i;
	double	x[3];
	double	first_vert_x[3];

	/* Get element geometry coordinate system
	-----------------------------------------*/
	if (elt_geom_csys_name_E == NULL) {
		if ((elt_geom_csys = find_csys(GLOBAL_NAME)) == NULL) {
			p_error(
				"Cannot find default coordinate system for elt_geom_csys",
				NULL);
		}
	} else {
		if ((elt_geom_csys = find_csys(elt_geom_csys_name_E))
			== NULL) {
			p_error("Coord sys given for elt_geom_csys was never defined",
				NULL);
		}
	}

	/* Print titles
	--------------*/
	fprintf(ofp_E,"\n\n===================================================\n");
	fprintf(ofp_E,    "ELEMENT GEOMETRY (Organized by object)\n");
	fprintf(ofp_E,    "COORD SYS: %s\n",elt_geom_csys->name);
	fprintf(ofp_E,    "===================================================\n");

	/* Loop over elements
	--------------------*/
	current_obj = first_obj_E;
	current_elt = first_elt_E;
	while (current_elt != NULL) {

		fprintf(ofp_E,"\n");

		if (current_elt == current_obj->first_elt) {
			elt_num = 1;
			fprintf(ofp_E,"OBJECT: %s\n\n",current_obj->name);
			fprintf(ofp_E,ELT_GEOM_LABELS);
			fprintf(ofp_E,ELT_GEOM_UNDLNS);
		}

		for (i=0; i < current_elt->num_vertices; i++) {
			copy_vector(current_elt->disloc_seg[i].vert[0]->x,x);
			transform_position_vector(FORWARD_TRANS,
				elt_geom_csys->origin, elt_geom_csys->local_rot, x);
			if (i == 0) 
				copy_vector(x,first_vert_x);
			fprintf(ofp_E,ELT_GEOM_FMT,elt_num,
				current_elt->disloc_seg[i].vert[0]->name,x[0],x[1],x[2]);
		}
		fprintf(ofp_E,ELT_GEOM_FMT,elt_num,
			current_elt->disloc_seg[0].vert[0]->name,first_vert_x[0],
				first_vert_x[1],first_vert_x[2]);

		if (current_elt == current_obj->last_elt)
			current_obj = current_obj->next;

		elt_num++;
		current_elt = current_elt->next;
	}
}


/*********************** Function: array_max_norm ****************************
* Returns the maximum matrix norm (L-infinity norm) of the matrix a.
* Adapted from a similar function by Ken C. Cruikshank.
*
* In:	a			- the matrix for which the norm is to be found
*		start_row	- index of the matrix's start row
*		end_row		- index of the matrix's end row
*		start_col	- index of the matrix's start column
*		end_col		- index of the matrix's end column
******************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
double   array_max_norm(double **a, int start_row, int end_row, int start_col,int end_col)
#else
double array_max_norm(a, start_row, end_row, start_col, end_col)
double	**a;
int		start_row;
int		end_row;
int		start_col;
int		end_col;
#endif

{
	int		i,j;
	double	norm;
	double	row_sum;

	norm = 0.0;
	for (i = start_row; i <= end_row; i++) {
		row_sum = 0.0;
		for (j = start_col; j <= end_col; j++) {
			row_sum += fabs(a[i][j]);
		}
		norm = MAX(norm,row_sum);
	}
	return(norm);
}


/*********************** Function: get_program_args *************************
* Prompts the user for input and output file names.  Only used if the
* symbolic constant FPROMPT is defined at compile time.
*****************************************************************************/
#ifdef FPROMPT
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     get_program_args(void)
#else
void get_program_args()
#endif

{
	int		numwords;
	char	*word[MAXWORDS];
	char	line[MAXLINE];

	/* Get the input file name
	--------------------------*/
	printf("\n INPUT FILE: ");
	numwords = read_line(line,word);
	switch (numwords) {
		case 0:
			printf(  "             (Using default value)\n");
			break;
		case 1:
			if (strlen(word[0]) > MAXFILE-1) {
				p_error("File name too long",NULL);
			}
			strcpy(infile_E,word[0]);
			break;
		default:
			p_error("Invalid file name",NULL);
	}

	/* Get the output file name
	---------------------------*/
	printf("OUTPUT FILE: ");
	numwords = read_line(line,word);
	switch (numwords) {
		case 0:
			printf(  "             (Using default value)\n");
			break;
		case 1:
			if (strlen(word[0]) > MAXFILE-1) {
				p_error("File name too long",NULL);
			}
			strcpy(outfile_E,word[0]);
			break;
		default:
			p_error("Invalid file name",NULL);
	}
	printf("\n");
}
#endif
