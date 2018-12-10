#ifndef HEAD_PREPROC
#define HEAD_PREPROC

#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include "mkl.h"
//#include "cgnslib.h"

#define STLEN 128
#define _FLOAT double

#define _SQRT sqrt
#define _FABS fabs
#define _POW pow
#define max(a,b) ((a)>(b) ? (a):(b))
#define min(a,b) ((a)<(b) ? (a):(b))
#define sign(a)  ((a)>=(0) ? (1):(-1))
#define trig(a)  ((a)>(0.0) ? (1):(0))
#define fsign(a) ((a)>=(0.0) ? (1.0):(-1.0))
#define ftrig(a) ((a)>(0.0) ? (1.0):(0.0))

#define NO_CLRUN 999
#define MAX_NITC 10
#define MAX_MGLEVEL 15
#define MAX_IPSTEP 101
#define MAX_ICOMP 10

/*Boundary condition types*/
#define _NAB 60   /*Not a boundary - for chimera purposes*/
#define _CON 70
#define _CHI 75
#define _SUI 79
#define _FAR 80
#define _WAL 90
#define _SYM 100

#define PI   (4.*atan(1.))
#define INF  1.0E300 /*Highest number available with double precision*/

/*Flow type*/
#define DROPLETS           20
#define AIR                21
#define THERMO             22
#define GEO                23

/*Boundary update models*/
#define RESIDUALS          0
#define CORRECTIONS        1
#define VARIABLES          2
#define MULTIGRID          3
#define VISCOSITIES        4
#define SA_TURBULENCE      5
#define KW_TURBULENCE      6
#define GRT_TURBULENCE     7
#define COPYVAR            8
#define ROUGH_TURBULENCE   9
#define SCALAR             10
#define MATRIX             11
#define INFWING            12


/*Banana and bapples*/
#define CELL_CENTERED      30
#define FACE_CENTERED      31
#define VERTEX_CENTERED    32
#define _FC_DISSIP         33
#define _CC_CONV           34
#define SOURCE             35
#define MESH               36
#define DEBUG              37

/*CHIMERA DEFINITION*/
#define CCW                1
#define CW                 -1

#define OV_HOLECUT         -1
#define OV_NOOVERLAP       0
#define OV_COMPUTED        1
#define OV_INTERPOLATED    2
#define OV_BUFFER          3


/*CGNS DEFINITION*/
#define INIT_FILE          40
#define UPD_FILE           41

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

/* Dimension definition as used by CGNS. 2 for 2D, 3 for 3D */
#define CELL_DIM 2
#define PHYS_DIM 2

/*****************************/
/*   structure prototypes    */
/*****************************/
typedef struct S_bag S_bag;
typedef struct S_mesh S_mesh;
typedef struct S_block S_block;
typedef struct S_box S_box;
typedef struct S_set S_set;
typedef struct S_faces S_faces;
typedef struct S_BCface S_BCface;
typedef struct S_subBCface S_subBCface;
typedef struct S_WALLface S_WALLface;
typedef struct S_INTERNALmesh S_INTERNALmesh;



/*********************************************************/
/* Division of the structures in NSCODE is as following  */
/*                                                       */
/* bag->blocks                                           */
/*        ->mesh.                      (MESH branch)     */
/*        ->faces->BCface->subBCface.  (TOPOLOGY branch) */
/*                                                       */
/* This division allows for a clear separation between   */
/* the mesh info and the topology info                   */
/*********************************************************/

struct S_bag
{
   /*run environment*/
   int ncpus;                 /*maximum number of threads enabled by the code*/
   
   /*input*/
   char title[STLEN];            /*case title*/
   char ctrlfilename[STLEN];     /*input file*/
   char topofilename[STLEN];     /*topology file*/
   char meshfilename[STLEN];     /*mesh file*/
   
   char version[STLEN];          /*NSCODE Version number in long tag format #.#.#.#-branch-node_short*/
   char datetime[STLEN];         /*date and time information*/
   char airflowfilename[STLEN];  /*airflow file*/
   int  nmesh;                   /*number of mesh used by chimera method*/
   char **topolist;              /*list of topology files to read for chimera method*/
   char **meshlist;              /*list of mesh files to read for chimera method*/
   int icing;                    /*Icing use: 1 if for icing, else 0 */

   /* output */
   int outputformat;             /*short/long ouput format*/
   int write_output_every;       /*Number of iterations between each writing of output files. 0 means no writing until the end of the simulation. */
   //Variable do deactivate NSCODE output
   int output;
   int output_type;              /*Type of output: 0-tecplot ; 1-CGNS ; 2-Both*/
   int output_location;          /*Location of the flow solution: 0-Vertex-centered ; 1-Cell-centered ; 2 - debug (vertex-centered)*/
   int cgns_level;               /*Level on which CGNS file is written, typically 0, but may be higher if only a coarser grid is run*/
   char cgnsfilename[STLEN];     /*CGNS output file*/
   
   /*chimera*/
   int chimera;                 /*Chimera method flag*/
   _FLOAT crit_number_per_box;  /*controls the number of boxes of the main boxing function*/
   int crit_number_subbox;      /*controls the number of boxes of the recursive boxing function*/
   int max_box_level;           /*controls the maximum refinment levels that can be used for the box (-1 for infinite)*/
   int max_box_level_used;
   int nsets;
   int *chierarchy;           /*hierarchy of the meshes*/
   int chlayer_size;          /*size of layer between overlapping meshes*/
   S_set *firstset;
   int ch_criteria;           /*criteria used to determine the boundary between two chimera blocks and donor cell amongst all overlaping cells*/
   int chinterpfct;           /*int used to determine the function that will be used for the interpolations*/
   int trimethod;             /*int used to determine the function that will be used for the calculation of a point inside a triangle*/

   
   /*profilling*/
   int level_boxing_memory[100]; /* Number of integer indices stored per boxing level */
   int total_boxing_memory;      /* Total number of integer indices stored in boxing */
   _FLOAT initialisation_time;   /* time required before iterations start */
   _FLOAT chimera_time;          /* time added by chimera pre-processing */
   _FLOAT box_time;              /* time required by thew boxing function */
   _FLOAT SO_time;               /* time required by the search_overlap function */
   _FLOAT internal_time;         /* time required by the internal mesh function */
   _FLOAT hole_time;             /* time required by the hole cutting function */
   _FLOAT wall_time;             /* time required by the local wall distance function */
   _FLOAT init_time;             /* time required by the init interpolation function */
   _FLOAT buffer_time;           /* time required by the add buffer function */

   /*blocks*/
   int nblks,prev_blk_n;      /*number of blocks*/
   S_block *firstBlock;       /*pointer to first element of the block chain list*/

	/* Wall faces */
	int *nwalfs; 												  /* Table of number of wall faces for each elements*/	
   int *walf_type;                                   /* Boolean containing the type of WALLface (closed - 1 or open - 0) */
	int *iTE[MAX_MGLEVEL],*jTE[MAX_MGLEVEL];		     /* Table of indices (i,j) of trailing edge for each elements (and each level) */	
	int *W_nbpt[MAX_MGLEVEL];                         /* Total number of cells on an icomp */
	struct S_WALLface **firstWALLface[MAX_MGLEVEL];   /* Table of pointers to first element of each wall face chain list (for each level) */

   /* internal mesh*/	 
   struct S_INTERNALmesh *firstINTERNALmesh[MAX_IPSTEP][MAX_MGLEVEL];  /* internal mesh of the geometry, for hole cutting*/
 
   /* flow & geometry properties */
   _FLOAT alpha;                                   /* geometry angle */
   _FLOAT roinf,ro_water,altitude;                 /* density and altitude*/
   _FLOAT tinf,mu_inf,ppinf;                       /* temperature, viscosity and pressure */
   _FLOAT Cpa;                                     /* Air heat capacity */
   _FLOAT mach,reynolds,prandtl_l,prandtl_t;       /* flow properties */
   _FLOAT MVD,LWC;                                 /* droplets properties */    
   _FLOAT xref,yref,cmac;                          /* forces/moments reference points */
   _FLOAT xtrans[MAX_ICOMP],ytrans[MAX_ICOMP];     /* Airfoil translation from (0,0) for icing multi-timesteps */
   _FLOAT beta,tu;
   _FLOAT ks;                                      /* Smooth roughness height */
   _FLOAT dt;                                      /* Simulation associated time interval in unsteady or quasi-steady simulations */
   int max_icomp,min_icomp;                        /* Ending and starting indice of the icomp elements*/
   int ppalt;													/* pressure or altitude input: altitude -> 0, pressure ->1*/

   /* Multi-grid */
   int coarsest_mg_level;
   int finest_mg_level;
   int coarsening;                  /* determines which type of grid coarsening to use */

   /* Infinite wing property */
   _FLOAT sweep;                 /* Sweep angle considered in 2.5D process by adding geometric sweep and flow sideslip angles */ 
   _FLOAT geom_sweep;            /* Geometric sweep angle */
   _FLOAT sideslip;              /* Flow sideslip angle. Positive sideslip adds to sweep effect. */
   _FLOAT eta   ; // eta = 0 => apply root correction ; eta = 1 => no root correction

   /*alpha ramping*/
   int    Rampalpha;          /* boolean that activates the ramping function for alpha */
   int    Rampmach;           /* boolean that activates the ramping function for mach */
   int    Ramp_iter_ini;      /* iteration at which ramping starts */
   int    Ramp_iter_end;      /* iteration at which ramping ends */
   _FLOAT Ramp_alpha_ini;     /* initial angle */
   _FLOAT Ramp_mach_ini;      /* initial mach */
   _FLOAT alpha0;             /* final angle */
   _FLOAT mach0;              /* final mach */

   /* Constant cl run*/
   _FLOAT cstcl_kp,cstcl_ki,cstcl_kd;

   /* Transition point calculation */
   _FLOAT dkdn_threshold;

   /* constants */
   _FLOAT gamma,eps,grav,Rair;
      
   /*NLFD*/
   int ipstep,curstep,nlfd_nout;       /* number of NLFD resolutions,curent timestep of the nlfd resolution, number of timestep to print for the output */
   int nlfd_initperturbation;          /* 0- no initialisation perturbation    1- perturbated initialisation*/
   int nlfd_VTP, nlfd_VTP_ini;
   _FLOAT nlfd_VTP_gain;
   _FLOAT nlfd_period,nlfd_omega;      /* period length of the nlfd resolution and corresponding omega value */
   _FLOAT nlfd_PeriodGrad;

   _FLOAT HB_TSV_vis;
   int HB_TSV_Harm;

   _FLOAT Phase, Phase_o, Phase_oo, frequency_o, frequency_oo, frequency_ooo;
   /* BiCGstab */
   int ipar[128],job[7];
   _FLOAT dpar[128];

};


struct S_block
{
	int block_id;                       /* determines the id of the block */
	int set_id;                         /* Differentiates blocks of the same mesh to lighten chimera preprocessing */ 
	struct S_set *set;                  /* Pointer to the appropriate set of the block */
	S_mesh *mesh[MAX_MGLEVEL];          /* MESHES: table containing pointers to the meshes (organized in MG levels) */
	S_faces *faces[MAX_MGLEVEL];        /* TOPOLOGY: table containing pointers to the faces (organized in MG levels) */
	S_box   *boxing[MAX_IPSTEP*MAX_MGLEVEL];       /* BOXING : contains the boxing of the block to accelerate grid search used in chimera (organized in MG levels) */
   char zonename[STLEN];               /* solely used for cgns output */
	struct S_block *nextBlock;          /* blocks are organized in chain list */
};

struct S_box                   /* structure storing boxing informations */
{
   int block_id;                       /* determines on which block the BCface is */
   int level;                          /* determines the level of multi-grid */
   int box_level;                      /* refinment level of the boxing */
   int nx,ny;                          /* number of elements of the boxing */
   _FLOAT dx,dy,xmin,xmax,ymin,ymax;   /* geometric definition of the boxing */
   int *sizelist;                      /* contains the start point of each box in boxlist */
   int *boxlist;                       /* contains the indices list contained in each box */
   
   struct S_box **fine_box;            /* each box may be refined with another boxing */
};

struct S_set                   /* structure storing mesh sets information */
{
   int set_id;                         /* id of the set */
   int nblks;                          /* number of blocks in the set */
   int *list_blks;                     /* list of the blocks of the set */
   _FLOAT Walextr[MAX_IPSTEP*4];       /* x/y parameters of the bounding box of the set */
   _FLOAT xref,yref;
   _FLOAT thetaref;
   struct S_set *nextset;              /* pointer to the next set of the linked list */
};

   /*note that S_faces should be declared first to respect structure ordering*/
   /*but for compilation purposes, S_BCface must be declared first*/
struct S_BCface                /* structure storing all the boundary conditions of a block face */
{
   int block_id;                       /* determines on which block the BCface is */
   int face_id;                        /* determines on which face the BCface is */
   int level;                          /* determines the level of multi-grid */
   int BCtype;                         /* determines the type of boundary condition */
   int BCtype2;                        /* determines priority between two BCface of the same type */
   int nsubBCface;                     /* number of subBCface on current BCface */
   struct S_subBCface *firstsubBCface; /* subBCfaces are organized in chain list */
};

struct S_faces                 /*structure storing all the topology of a given block on a given level*/
{
   int block_id;                       /* determines on which block the face is */
   int level;                          /* determines the level of multi-grid */
   int BCorder[4];                     /* table containing the order of priority of the faces */
   S_BCface BCface[4];                 /* table containing the faces */
   int nholecut[MAX_IPSTEP];           /* number of holecut cells in chimera method*/
   int ninterp[MAX_IPSTEP];            /* number of interpolated cells in chimera method */
   int *list_interp[MAX_IPSTEP];       /* table containing interpolated cells in chimera method */
   int *cell_type[MAX_IPSTEP];         /* table containing interpolated cells type (halo or not) in chimera method */
};

struct S_subBCface             /* structure caracterizing a single boundary condition */
{
   int block_id;                       /* determines on which block the subBCface is */
   int face_id;                        /* determines on which face the subBCface is */
   int level;                          /* determines the level of multi-grid */
   int type;                           /* determines the type of boundary */
   int icomp;                          /* determines to which part of the geometry the subBCface is related */
   int iflip;                          /* for 3D considerations */
   int indexBC[4];                     /* starting and ending indices of the subBCface */
   int indexCBC[4];                    /* starting and ending indices of the connected subBCface */
   int cblock_id;                      /* connecting block id */
   int cface_id;                       /* connecting face id */
   int iis,iie,incii;                  /* initialized in initial thermo field */
   int nbpt;                           /* initialized in initial thermo field // himax+1 or hjmax+1 depending on the type of face*/
   char BCname[STLEN];                 /* solely used for cgns output */
   _FLOAT *Cf,*Cp;                     /* friction and pressure coefficients */
   _FLOAT *Cfx,*Cfy,*Cfz;              /* Skin friction line */
   _FLOAT *Cpx,*Cpy,*Cpz;              /* Pressure line */
   _FLOAT *y_plus,*k_plus;             /* KW model turbulence variables */
   _FLOAT *s,*ds;                      /* Curvilinear distance around airfoil */
   _FLOAT *curv;                       /* 2D surface curvature */
   _FLOAT *weight;                     /* patch-grids panel weight for aerodynamic forces calculation */

   /// ICING variables ///
   _FLOAT *bb,*m_imp;                  /* beta : impingement fraction on airfoil */
   _FLOAT *Ts,*Ts0,*Ts00;              /* Surface temperature of ice or water*/
   _FLOAT *m_in_w,*m_in_e;             /* Incoming water runback mass from north and south */
   _FLOAT *m_out,*m_out0;              /* Total outgoing water runback mass */
   _FLOAT *m_ice,*m_ice0;              /* Ice mass */
   _FLOAT *m_st,*m_st0;                /* Water mass staying in control volume between each timestep */
   _FLOAT *m_ev,*m_su;                 /* Evaporative or sublimative water mass */
   _FLOAT *Ri_mm,*Ri_ee;               /* Mass and energy equations residual */
   _FLOAT *fw,*fe;                     /* Stress coefficients */
   _FLOAT *f;                          /* Freezing fraction */
   _FLOAT *hc;                         /* Convective heat coefficient */
   _FLOAT *ks;                         /* Roughness height */
   _FLOAT *Ue;                         /* Airspeed at end of boundary */
   _FLOAT *Qf;                         /* Aerodynamic heating flux */
   _FLOAT *h00,*h0,*h;                 /* Rectangular ice height and iterative ice height */
   _FLOAT *bsx,*bsy;                   /* bissectrix components */
   _FLOAT *XS;                         /* new geometry points coordinates and derivatives */
   
   /// Extended Messinger ///
   _FLOAT *iglaze;                     /* 1 if glaze, 0 if rime, -1 if dry*/
   _FLOAT *hitot;                      /* totat height of ice on all layers */
   _FLOAT *tg,*hig;                    /* time and ice height at which glaze ice first occur */

   /// SHALLOW WATER VARIABLES ///
   _FLOAT *Cft;                        /* Tangential friction coefficient */
   _FLOAT *uf,*hf,*hf0,*hf00;          /* water film mean velocity (uf: averaged along the height) and height (hf,hf0,hf00) */
   _FLOAT *mu_w;                       /* local viscosity of water */
   _FLOAT *U,*F;                       /* U and F matrices in hyperbolic system of eqns */
   _FLOAT *iF,*iFp,*iFm;               /* intercells fluxes for the hyperbolic system of eqns*/
   _FLOAT *dU,*dU0;                    /* Variation of U (Delta U)*/
   _FLOAT *S;                          /* Source Terms*/
    
   /// MESH ///
   S_mesh *mesh;
   struct S_subBCface *nextsubBCface;  /* subBCfaces are organized in chain list */
};

struct S_WALLface                     /*structure storing the limit points of a wall subBCface and connecting wall subBCfaces*/
{	
   int iistart;                        /* Absolute indice at which the Wall section starts on the icomp frm the TE */
   int nbpt;                           /* Vertex number of points on a WALLface (himax+1 or hjmax+1) */
   int indexWALL[4];                   /* starting and ending indices of the WALLface */
   int wall_id;                        /* determines the id of the WALLface */
   struct S_subBCface *subBCface;      /* subBCfaces contain all information */
   struct S_WALLface *nextWALLface;    /* WALLfaces are organized in chain list in geometrical order around elements*/

};

struct S_INTERNALmesh
{
   int nbnodes;                              /* Number of nodes of the internal mesh */
   int nbcells;                              /* Number of cells of the internal mesh */
   int set_id;                               /* set identification of the internal mesh */
   _FLOAT *x, *y;                            /* x and y coordinates of the nodes */
   _FLOAT *area;                             /* Area of the cells of the internal mesh */
   int **cells;                              /* Triangulation of the cells */
   struct S_INTERNALmesh *nextINTERNALmesh;
};

struct S_mesh
{
   // indexing variables //

   int rimax,rjmax;     /*restricted domain cc imax and jmax*/
   int himax,hjmax;     /*halo cc imax and jmax*/
   int inci,incj;       /*address increments in i and j*/
   int nbpt;            /*number of point of the mesh (cv)*/
   int ncell;           /*number of cell of the mesh (cc)*/
   int imax,jmax;       /*these terms really are never used in the code... they are the same as rimax,rjmax*/

   // geometrical variables //

   _FLOAT ****Xnlfd,**x_o,**y_o,**x_oo,**y_oo,**x0,**y0,**dvi,**dvj,**dvi_o,**dvj_o;            /*mesh coordinates and the volume sweep by cell edges*/
   _FLOAT **area,**area_o,**area_oo,**area_s;             /*cell area*/
   _FLOAT ****SIIXnlfd;  /*face i and j projections, in x and y directions respectively*/

   // flow variables //

   _FLOAT ****Wnlfd,****Whatnlfd;//***nlfd_ro,***nlfd_uu,***nlfd_vv,***nlfd_pp;    /*primitive variables for the NLFD approach*/
   _FLOAT ****MUUnlfd;                                  /*mesh speed, used for ALE mesh movement*/
   _FLOAT ***dSAdT_nlfd, ***dSAdTjacob_nlfd;
   _FLOAT ****FFnlfd;                        /*primitive variables*/
   _FLOAT **rocv,**uucv,**vvcv,**ppcv;                    /*primitive variables cell-vertex*/
   _FLOAT ****W0nlfd;//***nlfd_ro0,***nlfd_ru0,***nlfd_rv0,***nlfd_re0;/*conservative variables rk(0) for the NLFD approach*/
   _FLOAT ****W1nlfd;/*conservative variables rk(k-1) for the NLFD approach*/
   _FLOAT ****Rnlfd;                /*inviscid Residual, artificial Residual, viscous residual*/
   _FLOAT ****tmp_Wnlfd;            /*artificial Residual*/
   _FLOAT **dw1c,**dw2c,**dw3c,**dw4c,**dw5c,**dw6c;                    /*dummy arrays for variables or corrections*/
   _FLOAT ****dWnlfd;                        /*LUSGS corrections*/
   _FLOAT ****specnlfd;      /*spectral radii cell centered (speci,specj) and face centered (speciFC, specjFC)*/
   _FLOAT ***dtnlfd;                                 /*time step for flow equations*/  
   _FLOAT ****W4MGnlfd;                /*store multigrid primitive variables*/

   // infinite wing variables
   _FLOAT **ww_rhs,**ww_imp;
   _FLOAT **wwcv;

    // Special variables for the modified tridiagonal solver used for infinit wing
   _FLOAT **di_old,**dj_old;

   // icing droplets variables
   _FLOAT **aa,**ud,**vd,**wd;                     /*primitive droplets variables*/  
   _FLOAT **aacv,**udcv,**vdcv,**wdcv;             /*primitive droplets variables cell-vertex*/
   _FLOAT **tmp_s0,**tmp_s1,**tmp_s2,**tmp_s3;     /*droplets' source terms*/
   _FLOAT **dwd1c,**dwd2c,**dwd3c,**dwd4c;         /*dummy arrays for droplets variables or corrections*/
   _FLOAT **Sspec1,**Sspec2,**Sspecm;              /*droplets' source jacobian spectral radius */
   _FLOAT **aa4mg,**ud4mg,**vd4mg,**wd4mg;         /*store multigrid primitive droplets variables*/
   _FLOAT *grad_phi;                               /*Stores gradient components for a variable */
   _FLOAT ****dAi,****dAj;                          /*Face Centered Convective flux Jacobian*/

   // geometrical variables //

   _FLOAT ****XCnlfd;                  /*mesh cell center coordinates*/
   _FLOAT *y_plus,*k_plus;             /*y^+ and k^+ (normalized wall distance and roughness height)*/
   _FLOAT ****SIInlfd;                 /*cell distances in i and j directions*/
   _FLOAT **l;                         /*wall distance*/
   _FLOAT ****Ji,****Jj;               /*Coordinate jacobian determinant*/
   _FLOAT **det_Ji,**det_Jj;           /*Coordinate jacobian determinant(face center)*/
   _FLOAT *dxy,*Adxy;                  /* Coordinates difference matrix for gradient reconstruction used by droplets implicit solver */
   _FLOAT **Walextr;                   /*Geometric extremums of the wall boundary for current mesh*/
      
   // flow variables //

   _FLOAT ****dUdXnlfd;                           /*derivatif at vertex*/

//   INFINIT WING
     
   _FLOAT ****tmp2_Wnlfd;        /*viscid residual*/
   _FLOAT ****MUnlfd;               /*laminar nad eddy viscosities*/
   _FLOAT ****VORTnlfd,****STRAINnlfd;          /* rotation rate and strain rate tensor*/
   _FLOAT ****dUdXCnlfd;        /*velocity gradients at cell centers*/
   _FLOAT ****SODnlfd;                                /*vorticity, strain rate and dilatation magnitudes*/
   _FLOAT ****SPECVnlfd;         /*viscous spectral radii cell centered*/
   _FLOAT **specviFC,**specvjFC;                  /*viscous spectral radii face centered*/
   _FLOAT **ro_o,**uu_o,**vv_o,**pp_o,**ww_o,**sa_o;            /*primitive variables at old steps*/
   _FLOAT **ro_oo,**uu_oo,**vv_oo,**pp_oo,**ww_oo,**sa_oo;        /*primitive variables at old steps*/
   _FLOAT **ro_s,**uu_s,**vv_s,**pp_s,**ww_s,**sa_s;            /*primitive variables source terms*/

   // turbulence variables (spalart-allmaras) //

   _FLOAT **sscv,***sa_imp;                 /*spalart-allmaras variables*/
   _FLOAT ***dSA;
   //_FLOAT **sa_o,**sa_oo;                         /*spalart-allmaras variables at old steps*/

   // turbulence variables (kw-sst) // 

   _FLOAT ***kw,***kw_rhs,***kw_imp,****kw_imp_coupled,**F1,***sigma_kw;         /*kw variables*/
   _FLOAT ***kw_o,***kw_oo;                       /*kw variables at old steps*/
   _FLOAT ***kw_grt_rhs,****kw_grt_imp;
   // transition variables (gamma Re-teta) // 

   _FLOAT ***grt,***grt_rhs,***grt_imp,****grt_imp_coupled;           /*gamma Re_teta transition variables*/
   _FLOAT ***grt_o,***grt_oo; 
   // roughess equation //

   _FLOAT **rough,**rough_rhs,**rough_imp;        /*roughness equation variables*/

   // output variables //

   _FLOAT ***test,***testcv,****nlfd_outW,****nlfd_outX;
   _FLOAT ****nlfd_outMu, ***nlfd_outSAustd, ****nlfd_outMvv;        /*test variables for cell center and vertex (output)*/
   _FLOAT ***nlfd_outcksi;

   // Jacobians //

   _FLOAT ****Ai,****Aj;                          /*Face Centered Convective flux Jacobian*/
   _FLOAT ****Aci,****Acj;                        /*Cell Centered Convective flux Jacobian*/
   _FLOAT ****Avi,****Avj;                        /*Viscous flux Jacobian*/
   _FLOAT ****Avi_num,****Avj_num;                /*Numerical flux Jacobian*/

/* Shahin's New Variables */

   // geometrical variables // 

   _FLOAT **sx2,**sy2;                            /*face i projections for i box*/ 
   _FLOAT **sx3,**sy3;                            /*face j projections for i box*/ 
   _FLOAT **sx4,**sy4;                            /*face i projections for j box*/ 
   _FLOAT **sx5,**sy5;                            /*face j projections for j box*/ 

   // flow variables //

   _FLOAT **mui, **muj;                           /*viscosity at center of vertex*/

/* Kazem's New Variables */
   
   /*Chimera variables*/
   
   int ***NLFDksi;                          /*used to determine overlap type visually*/
   _FLOAT ***NLFDcksi,**cksicv;            /*used to determine overlap type computationally*/
   _FLOAT ****XCHnlfd;
   _FLOAT **ll;
   int ***NLFDdonor_block;          /*block containing the donor cell for the interpolation*/
   int ****NLFDdonor_cells;        /*donor cells used to interpolate the value of the cell value*/
   _FLOAT ****NLFDdonor_weight;   /*ponderation of each donor cells for the interpolation of the cell value*/

   /* Implicit Euler scheme*/
   int *JAIA,*dist;
   _FLOAT *work,*XB;
   _FLOAT **val;

   /* plasma source terms */
   _FLOAT **source_x,**source_y;
   
   /* Low Mach Preconditioning */
   _FLOAT *beta2;                                  /* Low Mach Preconditioning various variables */
   _FLOAT  ****PRECOND;
   _FLOAT  ****PRECONDI;

   /* Harmonic Balance */
   _FLOAT **HB_D, **HB_TSV;

};



/*****************************/
/* declare the nsc structure */
/*****************************/
extern S_bag *bag;


/*****************************/
/*   declare the functions   */
/*****************************/
extern void	      add_buffer(int level);
extern void       adjust_BCindex(int block_id,int level,int mod);
extern _FLOAT	  *allocate_1d_array(int ntot,char *str);
extern _FLOAT	 **allocate_2d_arrays(int imax,int jmax,char *str);
extern int	    **allocate_2d_int(int imax,int jmax,char *str);
extern _FLOAT  ***allocate_3d_arrays(int imax,int jmax,int kmax, char *str);
extern _FLOAT  ***allocate_3d_arrays3(int imax,int jmax,int kmax, char *str);
extern int     ***allocate_3d_int(int imax,int jmax,int kmax, char *str);
extern int    ****allocate_4d_int(int imax,int jmax,int kmax,int lmax,char *str);
extern _FLOAT ****allocate_4d_arrays(int imax,int jmax,int kmax,int lmax,char *str);
extern _FLOAT ****allocate_4d_arrays2(int imax,int jmax,int kmax,int lmax,char *str);
extern _FLOAT ****allocate_4d_arrays3(int imax,int jmax,int kmax,int lmax,char *str);
extern void       allocate_wall_variables(int nbpt,S_subBCface *wksubBCface);
extern void       BC_CON_cv(S_subBCface *subBCface,_FLOAT *var,_FLOAT *cvar,int nvar, int nhalo);
extern void	      BC_SYM_mesh(S_subBCface *subBCface,_FLOAT *var,int nvar,int cv,int nhalo);
extern void	      BC_WAL_mesh(S_subBCface *subBCface,_FLOAT *var,int nvar,int nhalo);
extern int        smbicgstab(int *iter,_FLOAT *resid,int n,char *matdescra,_FLOAT *A,int lval,int ndiag,int *dist,char *predescra,_FLOAT *b,_FLOAT *x);
extern void	      boxing(int level);
extern void	      checkmg(int blk);
extern void	      checktopo(int level);
extern void	      CHI_donor(int level);
extern void	      chi_interp(void);
extern void	      chi_wall_extremums(void);
extern void	      chimera_output(void);
extern void	      chimera4mg(int level);
extern void	      create_topology(int init);
extern void	      curv_distance(int level);
extern void       curv_distance(int level);
extern void	      d2wall(int level, int block_id);
extern void	      d2wall_local(int level);
extern void	      deallocate_2d_arrays(_FLOAT **tmp);
extern void       deallocate_2d_int(int **tmp);
extern void	      deallocate_3d_arrays(_FLOAT ***tmp);
extern void	      deallocate_4d_arrays(_FLOAT ****tmp,int imax,int jmax,int kmax,int lmax,char *str);
extern void	      deallocate_boxing(int step,int level,int set_id);
extern void       deallocate_chi_preprocess();
extern void       eign(int n,_FLOAT **D,_FLOAT *Lambda);
extern void       fgemat2spdiamat(int nb,int ndiag,_FLOAT **A,_FLOAT **val,int *dist);
extern S_block    *find_block(int id);
extern S_WALLface *find_WALLface(int level,int icomp,int id);
extern void       gradient_geometry_matrix(int level);
extern void       hole_cutting(int level);
extern void	      init_interp(int level);
extern void	      init_set();
extern void	      initial_system(void);
extern void       internal_mesh(int level);
extern void	      interp_test(int level);
extern void	      intpbilin(_FLOAT *x,_FLOAT *y,int ia0,int inci,int incj,_FLOAT xin,_FLOAT yin,_FLOAT *COEFS,int blk,int donor_block,char *origin);
extern void	      intpbilin_old(_FLOAT *x,_FLOAT *y,int ia0,int inci,int incj,_FLOAT xin,_FLOAT yin,_FLOAT *COEFS,int blk,int donor_block,char *origin);
//extern void	    intpbilin_v3(_FLOAT *x,_FLOAT *y,int ia0,int inci,int incj,_FLOAT xin,_FLOAT yin,_FLOAT *COEFS,int blk,int donor_block,char *origin);
extern void	      intpinvd(_FLOAT *x,_FLOAT *y,int ia0,int inci,int incj,_FLOAT xin,_FLOAT yin,_FLOAT *COEFS, int blk,int donor_block,char *origin);
extern void	      intpquad(_FLOAT *x,_FLOAT *y,int ia0,int inci,int incj,_FLOAT xin,_FLOAT yin,_FLOAT *COEFS,int blk,int donor_block,char *origin);
extern void	      intptetra(_FLOAT *x,_FLOAT *y,int ia0,int inci,int incj,_FLOAT xin,_FLOAT yin,_FLOAT *COEFS,int blk,int donor_block,char *origin);
extern void       invert_mat(_FLOAT **a, int N, int tracking_num);
extern void       lubksb(_FLOAT **a, int n, int *indx, _FLOAT *b);
extern void       ludcmp(_FLOAT **a, int n, int *indx, _FLOAT *d,int tracking_num);
extern void       matrix_mult(_FLOAT **rtnmat,_FLOAT **A,_FLOAT **B,int n);
extern void	      mesh4cgrids(int step,int block_id);
extern void	      mesh4halos(int step,int level,int set);
extern void	      metric(int step,int level,int set);
extern void	      modBCindex4cc(int block_id);
extern void	      modBCindex4halos(int block_id);
extern void	      modBCindexEdge(S_BCface *faceA, int dirA, int minmaxA, S_BCface *faceB, int dirB, int minmaxB);
extern void	      move_geo(int step,int inset, _FLOAT indx, _FLOAT indy, _FLOAT indtheta);
extern void       multip_mat(int n,_FLOAT **L,_FLOAT **R);
extern void       multip_mat_vect(int n,_FLOAT **A,_FLOAT *B);
extern void	      new_block(int id);
extern S_box	  *new_box(int block_id,int level,int nx,int ny,_FLOAT xmin,_FLOAT xmax,_FLOAT ymin,_FLOAT ymax,int box_level);
extern S_faces   *new_faces(int block_id,int level);
extern void       new_INTERNALmesh(int icomp, int step, int level);
extern S_mesh	  *new_mesh(int imax,int jmax,int level);
extern void	      new_set(int id);
extern void	      new_subBCface(int block_id,int level,int face_id);
extern void	      new_WALLface(S_subBCface *wksubBCface);
extern void       output_cgns_chimera(void);
extern void       output_cgns_grid(int setid, int step);
extern void       output_cgns_structure(int iop);
extern void       output_cgns_topo(void);
extern void       output_internal_mesh(int level);
extern int        point_inside(int method, _FLOAT xp, _FLOAT yp, _FLOAT x0, _FLOAT y0, _FLOAT x1, _FLOAT y1, _FLOAT x2, _FLOAT y2);
extern void       printerror(char *str);
extern void       printwarning(char *str);
extern void       print_mesh(int level);
extern void	      print_set_info();
extern void	      print_topo(void);
extern void	      prolong_chi(int level);
extern void       ramp_alpha(int iop,int init);
extern void	      readctrl(void);
extern void	      readmultimesh(char meshfilename[STLEN]);
extern void       readplasma();
extern void	      readtopo(char topofilename[STLEN],int set);
extern void	      recursive_boxing(S_box *wkbox,int step,int box_id,_FLOAT xmin,_FLOAT xmax,_FLOAT ymin,_FLOAT ymax);
extern void       reset_4d_arrays(int imax,int jmax, int kmax, int lmax, _FLOAT ****Mat);
extern void	      search_overlap(int level);
extern void	      subBCface_reordering(int blk);
extern void	      topo4cgrids(int block_id);
extern void       tridiagonal(int il,int iu,_FLOAT *b,_FLOAT *d,_FLOAT *a,_FLOAT *c);
extern void       tridiagonal2(int jl, int ju, int il,int iu,_FLOAT **aa,_FLOAT **bb,_FLOAT **cc,_FLOAT **dd);
extern void       tridiagonal2_mod(int jl, int ju, int il,int iu,_FLOAT **aa,_FLOAT **bb,_FLOAT **cc,_FLOAT **dd,_FLOAT **dd_mod);
extern void       tridiagonal_block(int il,int iu,int k,_FLOAT ***a,_FLOAT ***b,_FLOAT ***c,_FLOAT **d);
extern void       updatemultimesh(char meshfilename[STLEN]);
extern void	      wall_sections(int level);
extern void       wall_patching(int level);


#endif

