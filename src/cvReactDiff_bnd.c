/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2010/12/01 22:51:32 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODE.

 * OLD
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).

 * NEW
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * mesh, with simple polynomial initial profiles.

 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVBAND band linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Gib: effectively the values for i=0, i=MX+1, j=0 and j=MY+1 are to set 0.
 * First step is to allow these boundary values to take another (constant) value.
 * Next step is to accommodate an irregular region.
 * This requires a bi-mapping:
 *   eqtn. number ke <=> (i,j)
 * and a way to detect a boundary site:
 *   bdry(i,j) = bnd_value or NOT_BDRY
 * where NOT_BDRY is a number flagging the location as interior, e.g. pow(2,20)
 * If there is only one boundary value (same value everywhere on the boundary)
 * then it would be simpler to make (i,j) map to an impossible eqtn. number,
 * e.g. (i,j) -> -1.
 * Note: the advective terms have been removed.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Header files with a description of contents used in cvbanx.c */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_band.h>        /* prototype for CVBand */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_band.h>  /* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS and EXP */

#define USE_JACK FALSE				/* this has virtually no effect on speed */

#define TWO  RCONST(2.0)
#define RTOL    RCONST(1.0e-5)    /* scalar relative tolerance (1.0e-5)*/
#define FLOOR   RCONST(1.0)     /* value of C1 or C2 at which tolerances (100.0)*/
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */

/* User-defined vector access macro IJth */

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. 
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the macro call vdata = NV_DATA_S(v),
   where v is an N_Vector. 
   The variables are ordered by the y index j, then by the x index i. */

#define IJth(vdata,i,j,data) (vdata[(j-1) + (i-1)*data->MY])
#define GMAP2(i,j,data) ((data)->gmap[(j) + (i)*((data)->MY+2)]);
#define GMAP3(i,j,k,data) ((data)->gmap[(k) + (j)*((data)->MZ+2) + (i)*((data)->MZ+2)*((data)->MY+2)]);

/* Type : UserData (contains grid constants) */
typedef struct {
	int NVARS;
	int ndim;
	int MX;
	int MY;
	int MZ;
	int NG;
	int NEQ;
	int lbw, ubw;
	realtype *vbnd;
	realtype xdcoef, ydcoef, zdcoef;	//, hacoef;		// NOTE: currently same diff. coeff. for each constituent!
	int *gmap;	/* kg = gmap[y + x*(MY+2)] x=0,..,MX+1; y=0,..,MY+1 */
	int *xmap;	/* x = xmap[k] */
	int *ymap;	/* y = ymap[k] */
	int *zmap;	/* z = zmap[k] */
	int *jack;
//	realtype *C;
//	realtype *dCdt;
//	realtype *dfdC;
} *UserData;

/* Private Helper Functions */

static int SetupUserData(UserData data);
static void TestMap(UserData data);
//static void SetIC(int ndim, int mx, int my, int mz, int NG, int nvars, 
//	int *xmap, int *ymap, int *zmap, realtype *v);
//static void	MakeMaps(int ndim, int mx, int my, int mz, int *NG, int **xmap, int **ymap, int **zmap);
static void PrintHeader(realtype t0, realtype reltol, realtype abstol, realtype umax,  UserData data);
static void PrintOutput(realtype t, realtype umax, long int nst);
static void PrintFinalStats(void *cvode_mem);
static void ShowSol(N_Vector u);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector u, N_Vector fu, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

extern void react(realtype *C, realtype *dCdt);
extern void reactjac(realtype C[], realtype dfdC[], int nvars); 

int n_f = 0;
N_Vector u_save;

//----------------------------------------------------------------------------
// The mappings from grid number kg to (x,y,z) are passed in xmap, ymap, zmap.
// The number of active grid cells is NG, and the number of constituents is nvars.
// These map values and the grid dimensions mx, my, mz must be such that:
//		1 <= x <= mx
//		1 <= y <= my
//		1 <= z <= mz
// The boundary values of the constituents (for grid locations outside the
// active region) are passed in vbnd[].
// The starting values of the variables are in v[], which holds the solution
// when the execution is completed.
//----------------------------------------------------------------------------
int solve(int ndim, int mx, int my, int mz, int NG, int nvars, realtype *vbnd, realtype *cdiff, 
	int *xmap, int *ymap, int *zmap, realtype *v, realtype t0, realtype t1, realtype dtout, int nout)
{
	int k, ic, ierr;
	realtype reltol, abstol, t, tout, umax;
	N_Vector u;
	realtype *udata;
	UserData data;
	void *cvode_mem;
	int iout, flag;
	long int nst;
	clock_t begin, end;
//	double time_spent;

	u = NULL;
	data = NULL;
	cvode_mem = NULL;
	reltol = RTOL;		/* Set the tolerances */
	abstol = ATOL;

	data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
	if(check_flag((void *)data, "malloc", 2)) return(1);
	data->NVARS = nvars;
	data->vbnd = (realtype *)malloc(data->NVARS*sizeof(realtype));
//	data->C = (realtype *)malloc(data->NVARS*sizeof(realtype));
//	data->dCdt = (realtype *)malloc(data->NVARS*sizeof(realtype));
//	data->dfdC = (realtype *)malloc(data->NVARS*sizeof(realtype));
	for (ic=0; ic<data->NVARS; ic++) {
		data->vbnd[ic] = vbnd[ic];
	}
	data->ndim = ndim;
	data->MX = mx;
	data->MY = my;
	data->MZ = mz;
	data->xmap = xmap;
	data->ymap = ymap;
	data->zmap = zmap;
	data->NG = NG;
	data->NEQ = data->NG*data->NVARS;
	data->xdcoef = cdiff[0];	//CDIFF/(dx*dx);
	data->ydcoef = cdiff[1];	//CDIFF/(dy*dy);
	data->zdcoef = cdiff[2];	//CDIFF/(dz*dz);
	ierr = SetupUserData(data);
	if (ierr) {
		printf("Error: data setp failed\n");
		return(1);
	}
	printf("NEQ: %d\n",data->NEQ);

	TestMap(data);

	/* Create a serial vector */
	u = N_VNew_Serial(data->NEQ);  /* Allocate u vector */
	if(check_flag((void*)u, "N_VNew_Serial: u", 0)) return(1);

	udata = NV_DATA_S(u);
	for (k=0; k<data->NEQ; k++)
		udata[k] = v[k];
//	SetIC(u, data);  /* Initialize u vector */

	printf("ShowSol\n");
	ShowSol(u);

	printf("CVodeCreate\n");
	/* Call CVodeCreate to create the solver memory and specify the 
	* Backward Differentiation Formula and the use of a Newton iteration */
//	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
	if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	printf("CVodeInit\n");
	/* Call CVodeInit to initialize the integrator memory and specify the
	* user's right hand side function in u'=f(t,u), the inital time T0, and
	* the initial dependent variable vector u. */
	flag = CVodeInit(cvode_mem, f, t0, u);
	if(check_flag(&flag, "CVodeInit", 1)) return(1);

	printf("CVodeSStolerances\n");
	/* Call CVodeSStolerances to specify the scalar relative tolerance
	* and scalar absolute tolerance */
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

	printf("CVodeSetUserData\n");
	/* Set the pointer to user-defined data */
	flag = CVodeSetUserData(cvode_mem, data);
	if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

	printf("CVBand\n");
	/* Call CVBand to specify the CVBAND band linear solver */
	flag = CVBand(cvode_mem, data->NEQ, data->MY, data->MY);
	if(check_flag(&flag, "CVBand", 1)) return(1);

	printf(" CVDlsSetBandJacFn\n");
	/* Set the user-supplied Jacobian routine Jac */
	flag = CVDlsSetBandJacFn(cvode_mem, Jac);
	if(check_flag(&flag, "CVDlsSetBandJacFn", 1)) return(1);

	/* In loop over output points: call CVode, print results, test for errors */

	printf("Start solving\n");
	umax = N_VMaxNorm(u);
	PrintHeader(t0, reltol, abstol, umax, data);
	begin = clock();
	for (iout=1, tout=t1; iout <= nout; iout++, tout += dtout) {
		flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
		if(check_flag(&flag, "CVode", 1)) break;
		umax = N_VMaxNorm(u);
		flag = CVodeGetNumSteps(cvode_mem, &nst);
		check_flag(&flag, "CVodeGetNumSteps", 1);
		PrintOutput(t, umax, nst);
	}

//    end = clock();
//	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
//	printf("Solve time: %6.1f\n",time_spent);

	ShowSol(u);

	PrintFinalStats(cvode_mem);  /* Print some final statistics   */

	//  if (u_save == NULL) {
	//	u_save = N_VNew_Serial(data->NEQ);  /* Allocate u vector */
	//	if(check_flag((void*)u_save, "N_VNew_Serial: u_save", 0)) return(1);
	//  }
	//  printf("u_save: %p\n",u_save);
	//  Save(u, data);

	for (k=0; k<data->NEQ; k++)
		v[k] = udata[k];

	N_VDestroy_Serial(u);   /* Free the u vector */
	CVodeFree(&cvode_mem);  /* Free the integrator memory */
	/* Free the user data */
//	free(data->xmap);
//	free(data->ymap);
//	if (data->zmap != NULL)
//		free(data->zmap);
	data->xmap = NULL;
	data->ymap = NULL;
	data->zmap = NULL;
	free(data->vbnd);
//	free(data->C);
//	free(data->dCdt);
//	free(data->dfdC);
	free(data->gmap);
	free(data);
	return(0);
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
static void ShowSol(N_Vector u)
{
	int ke;
	realtype *udata;
	udata = NV_DATA_S(u);
	for (ke=0; ke<10; ke++)
		printf("%7.4f ",udata[ke]);
	printf("\n");
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

//----------------------------------------------------------------------------
// f routine. Compute du/dt = f(t,u). 
//----------------------------------------------------------------------------
static int f(realtype t, N_Vector u,N_Vector udot, void *user_data)
{
  realtype uij, udn, uup, ult, urt, ufd, ubk;
  realtype xdc, ydc, zdc, xdiff, ydiff, zdiff;	//, horac, hadv;
  realtype *C, *dCdt;
  realtype *udata, *dudata;
  realtype dumax;
  int NVARS;
  int i, j, k;
  int kg, ic, ke, kup, kdn, klt, krt, kbk, kfd, kjack;
  int  kgup, kgdn, kglt, kgrt, kgbk, kgfd;
  UserData data;

  n_f++;
  dumax = 0;
  //printf("f\n");
  udata = NV_DATA_S(u);
  dudata = NV_DATA_S(udot);

  /* Extract needed constants from data */

  data = (UserData) user_data;
  NVARS = data->NVARS;
  C = (realtype *)malloc(NVARS*sizeof(realtype));
  dCdt = (realtype *)malloc(NVARS*sizeof(realtype));
  xdc = data->xdcoef;
  ydc = data->ydcoef;
  zdc = data->zdcoef;
//  horac = data->hacoef;
 
  /* Loop over all eqtns */
  if (data->ndim == 2) {
	  for (kg=0; kg<data->NG; kg++) {
		  i = data->xmap[kg];
		  j = data->ymap[kg];
		  k = GMAP2(i,j,data);
//		  printf("i,j,k: %d %d %d\n",i,j,k);
		  kgdn = GMAP2(i,j-1,data);
		  kgup = GMAP2(i,j+1,data);
		  kglt = GMAP2(i-1,j,data);
		  kgrt = GMAP2(i+1,j,data);
//		  printf("%d %d %d %d   %d %d %d %d\n",kg,i,j,k,kgdn,kgup,kglt,kgrt);
		  for (ic=0; ic<NVARS; ic++) {
			  ke = kg*NVARS + ic;
			  C[ic] = udata[ke];
		  }
		  react(C,dCdt);
		  // Set diffusion and advection terms and load into udot
		  for (ic=0; ic<NVARS; ic++) {
			  ke = kg*NVARS + ic;
			  uij = udata[ke];
			  if (kgdn < 0)
				  udn = data->vbnd[ic];
			  else {
				  kdn = kgdn*NVARS + ic;
				  udn = udata[kdn];
			  }
			  if (kgup < 0)
				  uup = data->vbnd[ic];
			  else {
				  kup = kgup*NVARS + ic;
				  uup = udata[kup];
			  }
			  if (kglt < 0)
				  ult = data->vbnd[ic];
			  else {
				  klt = kglt*NVARS + ic;
				  ult = udata[klt];
			  }
			  if (kgrt < 0)
				  urt = data->vbnd[ic];
			  else {
				  krt = kgrt*NVARS + ic;
				  urt = udata[krt];
			  }
			  xdiff = xdc*(ult - TWO*uij + urt);
			  ydiff = ydc*(uup - TWO*uij + udn);
//			  hadv = horac*(urt - ult);
			  dudata[ke] = xdiff + ydiff + dCdt[ic];	// + hadv
		  }
//		  printf("did dudata loop: %d\n",ke);
	  }
  } else if (data->ndim == 3) {
	  for (kg=0; kg<data->NG; kg++) {
		  if (USE_JACK) {
			  kjack = 6*kg;
			  kgdn = data->jack[kjack++];
			  kgup = data->jack[kjack++];
			  kglt = data->jack[kjack++];
			  kgrt = data->jack[kjack++];
			  kgbk = data->jack[kjack++];
			  kgfd = data->jack[kjack++];
		  } else {
 			  i = data->xmap[kg];
			  j = data->ymap[kg];
			  k = data->zmap[kg];
			  kgdn = GMAP3(i,j-1,k,data);
			  kgup = GMAP3(i,j+1,k,data);
			  kglt = GMAP3(i-1,j,k,data);
			  kgrt = GMAP3(i+1,j,k,data);
			  kgbk = GMAP3(i,j,k-1,data);
			  kgfd = GMAP3(i,j,k+1,data);
		  }
		  for (ic=0; ic<NVARS; ic++) {
			  ke = kg*NVARS + ic;
			  C[ic] = udata[ke];
		  }
		  react(C,dCdt);
		  
		  // Set diffusion and advection terms and load into udot
		  for (ic=0; ic<NVARS; ic++) {
			  ke = kg*NVARS + ic;
			  uij = udata[ke];
			  if (kgdn < 0)
				  udn = data->vbnd[ic];
			  else {
				  kdn = kgdn*NVARS + ic;
				  udn = udata[kdn];
			  }
			  if (kgup < 0)
				  uup = data->vbnd[ic];
			  else {
				  kup = kgup*NVARS + ic;
				  uup = udata[kup];
			  }
			  if (kglt < 0)
				  ult = data->vbnd[ic];
			  else {
				  klt = kglt*NVARS + ic;
				  ult = udata[klt];
			  }
			  if (kgrt < 0)
				  urt = data->vbnd[ic];
			  else {
				  krt = kgrt*NVARS + ic;
				  urt = udata[krt];
			  }
			  if (kgbk < 0)
				  ubk = data->vbnd[ic];
			  else {
				  kbk = kgbk*NVARS + ic;
				  ubk = udata[kbk];
			  }
			  if (kgfd < 0)
				  ufd = data->vbnd[ic];
			  else {
				  kfd = kgfd*NVARS + ic;
				  ufd = udata[kfd];
			  }

			  xdiff = xdc*(ult - TWO*uij + urt);
			  ydiff = ydc*(udn - TWO*uij + uup);
			  zdiff = zdc*(ubk - TWO*uij + ufd);
//			  hadv = horac*(urt - ult);
			  dudata[ke] = xdiff + ydiff + zdiff + dCdt[ic];	// + hadv;
			  dumax = MAX(dumax,fabs(dudata[ke]));
		  }
	  }
  }
  free(C);
  free(dCdt);
	if (n_f%100 == 0) {
		kg = data->NG/2;
		ke = NVARS*kg;
		printf("At %d %d %d  dudata: %12.4e  dumax: %12.4e\n",data->xmap[kg],data->ymap[kg],data->zmap[kg],dudata[ke],dumax);
	}
//  printf("did f\n");
  return(0);
}


//----------------------------------------------------------------------------
// Jacobian routine. Compute J(t,u).
// The basic treatment of the diffusion Laplacian and the advective term is unchanged.
// Inter-constituent reactions add contributions relating the ic terms for a given kg.
//----------------------------------------------------------------------------
static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector u, N_Vector fu, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	int NVARS;
	int i, j, k, kg, ke, ic, kje, jc;
	int kup, kdn, klt, krt, kbk, kfd;
	int  kgup, kgdn, kglt, kgrt, kgbk, kgfd, kjack;
	realtype *kthCol, xdc, ydc, zdc;	//, horac;
	realtype *udata;
	realtype *C, *dfdC;
	realtype dfdC_elem;
	UserData data;
  
//	printf("Jac\n");
  /*
    The components of f = udot that depend on u(i,j) are:
	For 2D:
    f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
      df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
      df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
      df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
      df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
      df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
	Note that the 
	For 3D:
	f(i,j,k), f(i-1,j,k), f(i+1,j,k), f(i,j-1,k), f(i,j+1,k), f(i,j,k-1), f(i,j,k+1)
	with
      df(i,j,k)/du(i,j,k) = -2 (1/dx^2 + 1/dy^2 + 1/dz^2)
      df(i-1,j,k)/du(i,j,k) = 1/dx^2 + .25/dx  (if i > 1)
      df(i+1,j,k)/du(i,j,k) = 1/dx^2 - .25/dx  (if i < MX)
      df(i,j-1,k)/du(i,j,k) = 1/dy^2           (if j > 1)
      df(i,j+1,k)/du(i,j,k) = 1/dy^2           (if j < MY)
      df(i,j,k-1)/du(i,j,k) = 1/dz^2           (if k > 1)
      df(i,j,k+1)/du(i,j,k) = 1/dz^2           (if k < MZ)
  */

//	for (i=0; i<NVARS; i++) {
//		dfdC[i] = (realtype *)malloc(NVARS*sizeof(realtype));
//	}
	udata = NV_DATA_S(u);

	data = (UserData) user_data;
	NVARS = data->NVARS;
	C = (realtype *)malloc(NVARS*sizeof(realtype));
	dfdC = (realtype *)malloc(NVARS*NVARS*sizeof(realtype));
	xdc = data->xdcoef;
	ydc = data->ydcoef;
	zdc = data->zdcoef;
//	horac = data->hacoef;
  
	if (data->ndim == 2) {
		for (kg=0; kg<data->NG; kg++) {
			i = data->xmap[kg];
			j = data->ymap[kg];

			kgdn = GMAP2(i,j-1,data);
			kgup = GMAP2(i,j+1,data);
			kglt = GMAP2(i-1,j,data);
			kgrt = GMAP2(i+1,j,data);

			for (ic=0; ic<NVARS; ic++) {
				ke = kg*NVARS + ic;
				C[ic] = udata[ke];
			}
			reactjac(C,dfdC,NVARS);	// dfdC[k'][k] = df[k']/dC[k]

			for (ic=0; ic<NVARS; ic++) {
				ke = kg*NVARS + ic;
				// set the kth column of J
				kthCol = BAND_COL(J,ke);
//				BAND_COL_ELEM(kthCol,ke,ke) = -TWO*(ydc+xdc);
				for (jc=0; jc<NVARS; jc++) {
					kje = kg*NVARS + jc;
					dfdC_elem = dfdC[jc + ic*NVARS];
					if (jc == ic)
						BAND_COL_ELEM(kthCol,ke,ke) = -TWO*(ydc+xdc) + dfdC_elem;	//dfdC[ic][ic];
					else
						BAND_COL_ELEM(kthCol,kje,ke) = dfdC_elem;	//dfdC[ic][jc];
				}
				if (kgdn >= 0) {
					kdn = kgdn*NVARS + ic;
					BAND_COL_ELEM(kthCol,kdn,ke)  = ydc;
				}
				if (kgup >= 0) {
					kup = kgup*NVARS + ic;
					BAND_COL_ELEM(kthCol,kup,ke)  = ydc;
				}
				if (kglt >= 0) {
					klt = kglt*NVARS + ic;
					BAND_COL_ELEM(kthCol,klt,ke) = xdc;	// + horac;
				}
				if (kgrt >= 0) {
					krt = kgrt*NVARS + ic;
					BAND_COL_ELEM(kthCol,krt,ke) = xdc;	// - horac;
				}
			}
		}
	} else if (data->ndim == 3) {
	  
		for (kg=0; kg<data->NG; kg++) {
			if (USE_JACK) {
				kjack = 6*kg;
				kgdn = data->jack[kjack++];
				kgup = data->jack[kjack++];
				kglt = data->jack[kjack++];
				kgrt = data->jack[kjack++];
				kgbk = data->jack[kjack++];
				kgfd = data->jack[kjack++];
			} else {
				i = data->xmap[kg];
				j = data->ymap[kg];
				k = data->zmap[kg];
				kgdn = GMAP3(i,j-1,k,data);
				kgup = GMAP3(i,j+1,k,data);
				kglt = GMAP3(i-1,j,k,data);
				kgrt = GMAP3(i+1,j,k,data);
				kgbk = GMAP3(i,j,k-1,data);
				kgfd = GMAP3(i,j,k+1,data);
			}

			for (ic=0; ic<NVARS; ic++) {
				ke = kg*NVARS + ic;
				C[ic] = udata[ke];
			}
			reactjac(C,dfdC,NVARS);	// dfdC[k'][k] = df[k']/dC[k]

			for (ic=0; ic<NVARS; ic++) {
				ke = kg*NVARS + ic;

			  // set the kth column of J
				kthCol = BAND_COL(J,ke);
//				BAND_COL_ELEM(kthCol,ke,ke) = -TWO*(zdc+ydc+xdc);
				for (jc=0; jc<NVARS; jc++) {
					kje = kg*NVARS + jc;
					dfdC_elem = dfdC[jc + ic*NVARS];
					if (jc == ic)
						BAND_COL_ELEM(kthCol,ke,ke) = -TWO*(zdc+ydc+xdc) + dfdC_elem;	//dfdC[ic][ic];
					else
						BAND_COL_ELEM(kthCol,kje,ke) = dfdC_elem;	//dfdC[ic][jc];
				}
				if (kgdn >= 0) {
					kdn = kgdn*NVARS + ic;
					BAND_COL_ELEM(kthCol,kdn,ke)  = ydc;
				}
				if (kgup >= 0) {
					kup = kgup*NVARS + ic;
					BAND_COL_ELEM(kthCol,kup,ke)  = ydc;
				}
				if (kglt >= 0) {
					klt = kglt*NVARS + ic;
					BAND_COL_ELEM(kthCol,klt,ke) = xdc;  //+ horac;
				}
				if (kgrt >= 0) {
					krt = kgrt*NVARS + ic;
					BAND_COL_ELEM(kthCol,krt,ke) = xdc;  //- horac;
				}
				if (kgbk >= 0) {
					kbk = kgbk*NVARS + ic;
					BAND_COL_ELEM(kthCol,kbk,ke) = zdc;
				}
				if (kgfd >= 0) {
					kfd = kgfd*NVARS + ic;
					BAND_COL_ELEM(kthCol,kfd,ke) = zdc;
				}
			}
		} 
	}
	free(C);
	free(dfdC);
//	printf("did Jac\n");
	return(0);
}

 //--------------------------------------------------------------
 // Private helper functions
 //--------------------------------------------------------------



//------------------------------------------------------------------------------
// Now there are NVARS variables at each (x,y,z), therefore the 
// constituent index is also needed.
// Grid index kg is found from:
//	n2 = y + x*(data->MY+2)
//	n3 = z + y*(data->MZ+2) + x*(data->MZ+2)*(data->MY+2)
// Variable index for constituent ic, where ic=0,..,NVARS-1
// is ke = kg*NVARS + ic
//------------------------------------------------------------------------------
static int SetupUserData(UserData data)
{
	int x, y, z, i, j, k, kg, n;
	int kup, kdn, klt, krt, kbk, kfd, kjack;
	int lbw, ubw;

	printf("SetupUserData\n");

//	data->hacoef = CADV/(TWO*dx);

// Generate gmap
	if (data->ndim == 2) {
		n = (data->MX+2)*(data->MY+2);
	} else if (data->ndim == 3) {
		n = (data->MX+2)*(data->MY+2)*(data->MY+2);
	}
	data->gmap = (int *)malloc(n*sizeof(int));
	for (k=0; k<n; k++)
		data->gmap[k] = -1;

	lbw = 0;
	ubw = 0;
	if (data->ndim == 2) {
//		printf("gmap: NG: %d MX: %d MY: %d\n",data->NG,data->MX,data->MY);
		for (kg=0; kg<data->NG; kg++) {
			x = data->xmap[kg];
			y = data->ymap[kg];
			k = y + x*(data->MY+2);
			data->gmap[k] = kg;
//			printf("kg, x, y, k: %d %d %d %d\n",kg,x,y,k);
		} 

		for (kg=0; kg<data->NG; kg++) {
			i = data->xmap[kg];
			j = data->ymap[kg];
			kdn = GMAP2(i,j-1,data);
			if (kdn > 0) {
				ubw = MAX(kg-kdn, ubw);
				lbw = MAX(kdn-kg, lbw);
//				printf("kdn: %d %d %d   %d  %d\n",kg,i,j,kdn,lbw);
			}
			kup = GMAP2(i,j+1,data);
			if (kup > 0) {
				ubw = MAX(kg-kup, ubw);
				lbw = MAX(kup-kg, lbw);
//				printf("kup: %d %d %d   %d  %d\n",kg,i,j,kup,lbw);
			}
			klt = GMAP2(i-1,j,data);
			if (klt > 0) {
				ubw = MAX(kg-klt, ubw);
				lbw = MAX(klt-kg, lbw);
//				printf("klt: %d %d %d   %d  %d\n",kg,i,j,klt,lbw);
			}
			krt = GMAP2(i+1,j,data);
			if (krt > 0) {
				ubw = MAX(kg-krt, ubw);
				lbw = MAX(krt-kg, lbw);
//				printf("krt: %d %d %d   %d  %d\n",kg,i,j,krt,lbw);
			}
		}
		printf("lbw, ubw: %d %d\n",lbw,ubw);
	} else if (data->ndim == 3) {
		for (kg=0; kg<data->NG; kg++) {
			x = data->xmap[kg];
			y = data->ymap[kg];
			z = data->zmap[kg];
			k = z + y*(data->MZ+2) + x*(data->MY+2)*(data->MZ+2);
			data->gmap[k] = kg;
		} 		
		if (USE_JACK) {
			data->jack = (int *)malloc(6*data->NG*sizeof(int));
		}
		kjack = 0;
		for (kg=0; kg<data->NG; kg++) {
			i = data->xmap[kg];
			j = data->ymap[kg];
			k = data->zmap[kg];
			kdn = GMAP3(i,j-1,k,data);
			if (kdn > 0) {
				ubw = MAX(kg-kdn, ubw);
				lbw = MAX(kdn-kg, lbw);
			}
			kup = GMAP3(i,j+1,k,data);
			if (kup > 0) {
				ubw = MAX(kg-kup, ubw);
				lbw = MAX(kup-kg, lbw);
			}
			klt = GMAP3(i-1,j,k,data);
			if (klt > 0) {
				ubw = MAX(kg-klt, ubw);
				lbw = MAX(klt-kg, lbw);
			}
			krt = GMAP3(i+1,j,k,data);
			if (krt > 0) {
				ubw = MAX(kg-krt, ubw);
				lbw = MAX(krt-kg, lbw);
			}
			kbk = GMAP3(i,j,k-1,data);
			if (kbk > 0) {
				ubw = MAX(kg-kbk, ubw);
				lbw = MAX(kbk-kg, lbw);
			}
			kfd = GMAP3(i,j,k+1,data);
			if (kfd > 0) {
				ubw = MAX(kg-kfd, ubw);
				lbw = MAX(kfd-kg, lbw);
			}
			if (USE_JACK) {
				data->jack[kjack++] = kdn;
				data->jack[kjack++] = kup;
				data->jack[kjack++] = klt;
				data->jack[kjack++] = krt;
				data->jack[kjack++] = kbk;
				data->jack[kjack++] = kfd;
			}
		}
		printf("NG: %d  lbw, ubw: %d %d\n",data->NG,lbw,ubw);	
	}
	data->lbw = lbw;
	data->ubw = ubw;
	return(0);
}

//----------------------------------------------------------------------------
// To test generation of gmap from xmap, ymap, zmap.
// This is successful, showing that gmap[] can be constructed from xmap, ymap, zmap.
// This means that we need pass only xmap, ymap and zmap for an aribitrary geometry.
//----------------------------------------------------------------------------
static void TestMap( UserData data)
{
	int *gmap_check;
	int kg, x, y, z, k;

	if (data->ndim == 2) {
		gmap_check = (int *)malloc((data->MX+2)*(data->MY+2)*sizeof(int));
		for (k=0; k<(data->MX+2)*(data->MY+2); k++)
			gmap_check[k] = -1;
		for (kg=0; kg<data->NG; kg++) {
			x = data->xmap[kg];
			y = data->ymap[kg];
			k = y + x*(data->MY+2);
			gmap_check[k] = kg;
		} 
		for (k=0; k<(data->MX+2)*(data->MY+2); k++) {
			if (gmap_check[k] != data->gmap[k]) {
				printf("gmap check fails: %d %d %d\n",k,gmap_check[k],data->gmap[k]);
				exit(1);
			}
		}
	} else if (data->ndim == 3) {
		gmap_check = (int *)malloc((data->MX+2)*(data->MY+2)*(data->MZ+2)*sizeof(int));
		for (k=0; k<(data->MX+2)*(data->MY+2)*(data->MZ+2); k++)
			gmap_check[k] = -1;
		for (kg=0; kg<data->NG; kg++) {
			x = data->xmap[kg];
			y = data->ymap[kg];
			z = data->zmap[kg];
			k = z + y*(data->MZ+2) + x*(data->MY+2)*(data->MZ+2);
			gmap_check[k] = kg;
		} 
		for (k=0; k<(data->MX+2)*(data->MY+2)*(data->MZ+2); k++) {
			if (gmap_check[k] != data->gmap[k]) {
				printf("gmap check fails: %d %d %d\n",k,gmap_check[k],data->gmap[k]);
				exit(1);
			}
		}
	}
	printf("gmap check OK\n");
}



/* Print first lines of output (problem description) */

static void PrintHeader(realtype t0, realtype reltol, realtype abstol, realtype umax,  UserData data)
{
  printf("\n2-D Advection-Diffusion Equation\n");
  printf("Mesh dimensions = %d X %d\n", data->MX, data->MY);
  printf("Total system size = %d\n", data->NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters: reltol = %Lg   abstol = %Lg\n\n",
         reltol, abstol);
  printf("At t = %Lg      max.norm(u) =%14.6Le \n", T0, umax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters: reltol = %lg   abstol = %lg\n\n",
         reltol, abstol);
  printf("At t = %lg      max.norm(u) =%14.6le \n", t0, umax);
#else
  printf("Tolerance parameters: reltol = %g   abstol = %g\n\n", reltol, abstol);
  printf("At t = %g      max.norm(u) =%14.6e \n", T0, umax);
#endif

  return;
}

/* Print current value */

static void PrintOutput(realtype t, realtype umax, long int nst)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %4.2Lf   max.norm(u) =%14.6Le   nst = %4ld\n", t, umax, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %4.2f   max.norm(u) =%14.6le   nst = %4ld\n", t, umax, nst);
#else
  printf("At t = %4.2f   max.norm(u) =%14.6e   nst = %4ld\n", t, umax, nst);
#endif

  return;
}

/* Get and print some final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  int flag;
  long int nst, nfe, nsetups, netf, nni, ncfn, nje, nfeLS;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 nni, ncfn, netf);

  return;
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
