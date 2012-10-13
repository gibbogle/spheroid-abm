// Main program to test solver
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Problem Constants */

// double precision case
typedef double realtype;

#define TRUE 1
#define FALSE 0

#define SPHERE TRUE
#define RADIUS 15
#define K01	  0.001
//#define CADV  0.0			 /* coeff of advection term */

#define XMIN         0           /* grid boundaries in x  */
#define XMAX         20.0
#define YMIN         30.0        /* grid boundaries in y  */
#define YMAX         50.0
#define ZMIN		 0
#define ZMAX		 20.0

static void SetIC(int ndim, int mx, int my, int mz, int NG, int nvars, int *xmap, int *ymap, int *zmap, realtype *v);
static void	MakeMaps(int ndim, int mx, int my, int mz, int *NG, int **xmap, int **ymap, int **zmap);
void React(realtype *C, realtype *dCdt);
void ReactJac(realtype C[], realtype dfdC[], int nvars);

extern int Solve(int ndim, int mx, int my, int mz, int NG, int nvars, realtype *vbnd, realtype *cdiff, 
	int *xmap, int *ymap, int *zmap, realtype *v, realtype t0, realtype t1, realtype dtout, int nout);

//------------------------------------------------------------------------------
// This function is to generate test geometries.
// In actual use, NG, xmap, ymap, zmap will be passed from the invoking program.
// ndim, MX, MY, MZ will also be supplied.
// Note that all "inside" gridcells must have 1 <= x <= MX etc.
//------------------------------------------------------------------------------
static void	MakeMaps(int ndim, int MX, int MY, int MZ, int *NG, int **xmap, int **ymap, int **zmap)
{
	int x, y, z, kg, nin, n2, n3;
	realtype x0, y0, z0, r2;
	int *inside, in;

	if (ndim == 2) {
		if (SPHERE) {
			x0 = (MX + 1)/2.0;
			y0 = (MY + 1)/2.0;
			inside = (int *)malloc((MX+2)*(MY+2)*sizeof(int));
			nin = 0;
			for (x=0; x<MX+2; x++) {
				for (y=0; y<MY+2; y++) {
					r2 = (x - x0)*(x-x0) + (y-y0)*(y-y0);
					if (r2 < RADIUS*RADIUS) {
						inside[y + x*(MY+2)] = TRUE;
						nin++;
					} else {
						inside[y + x*(MY+2)] = FALSE;
					}
				}
			}
		} else {
			nin = MX*MY;
		}
		*xmap = (int *)malloc(nin*sizeof(int));
		*ymap = (int *)malloc(nin*sizeof(int));
		*zmap = NULL;
		*NG = nin;
//		data->gmap = (int *)malloc((data->MX+2)*(data->MY+2)*sizeof(int));
		kg = -1;
		for (x=0; x<MX+2; x++) {
			for (y=0; y<MY+2; y++) {
				n2 = y + x*(MY+2);
				in = FALSE;
				if (SPHERE) {
					if (inside[n2]) in = TRUE;
				} else {
					if (x >= 1 && x <= MX && y >= 1 && y <= MY) in = TRUE;
				}
				if (in) {
					kg++;
					(*xmap)[kg] = x;
					(*ymap)[kg] = y;
				}
			}
		}
		// kg is the grid index in the list.
		if (SPHERE) {
			free(inside);
		}
	} else if (ndim == 3) {
		if (SPHERE) {
			x0 = (MX + 1)/2.0;
			y0 = (MY + 1)/2.0;
			z0 = (MZ + 1)/2.0;
			inside = (int *)malloc((MX+2)*(MY+2)*(MZ+2)*sizeof(int));
			nin = 0;
			for (x=0; x<MX+2; x++) {
				for (y=0; y<MY+2; y++) {
					for (z=0; z<MZ+2; z++) {
						r2 = (x - x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0);
						if (r2 < RADIUS*RADIUS) {
							inside[z + y*(MZ+2) + x*(MZ+2)*(MY+2)] = TRUE;
							nin++;
						} else {
							inside[z + y*(MZ+2) + x*(MZ+2)*(MY+2)] = FALSE;
						}
					}
				}
			}
		} else {
			nin = MX*MY*MZ;
		}
		*xmap = (int *)malloc(nin*sizeof(int));
		*ymap = (int *)malloc(nin*sizeof(int));
		*zmap = (int *)malloc(nin*sizeof(int));
		*NG = nin;
//		data->gmap = (int *)malloc((data->MX+2)*(data->MY+2)*(data->MZ+2)*sizeof(int));
		kg = -1;
		for (x=0; x<MX+2; x++) {
			for (y=0; y<MY+2; y++) {
				for (z=0; z<MZ+2; z++) {
					n3 = z + y*(MZ+2) + x*(MZ+2)*(MY+2);
					in = FALSE;
					if (SPHERE) {
						if (inside[n3]) in = TRUE;
					} else {
						if (x >= 1 && x <= MX 
							&& y >= 1 && y <= MY
							&& z >= 1 && z <= MZ) in = TRUE;
					}
					if (in) {
						kg++;
						(*xmap)[kg] = x;
						(*ymap)[kg] = y;
						(*zmap)[kg] = z;
					}
				}
			}
		}
		if (SPHERE) {
			free(inside);
		}
	}
}

//----------------------------------------------------------------------------
// Set test-case initial conditions in v
//----------------------------------------------------------------------------
static void SetIC(int ndim, int mx, int my, int mz, int NG, int nvars, 
	int *xmap, int *ymap, int *zmap, realtype *v)
{
	int kg, ke;
	realtype x, y, z, dx, dy, dz;
	realtype xmid, ymid, zmid, ax, ay, az;

	dx = (XMAX-XMIN)/mx;  
	dy = (YMAX-YMIN)/my;
	xmid = 0.5*(XMIN + XMAX);
	ymid = 0.5*(YMIN + YMAX);
	if (ndim == 2) {
		for (kg=0; kg<NG; kg++) {
			x = XMIN + dx*(xmap[kg] - 0.5);
			y = YMIN + dy*(ymap[kg] - 0.5);
			ke = kg*nvars;
			if (x < xmid) {
				ax = (x-XMIN)/(xmid-XMIN);
			} else {
				ax = (XMAX-x)/(XMAX-xmid);
			}
			if (y < ymid) {
				ay = (y-YMIN)/(ymid-YMIN);
			} else {
				ay = (YMAX-y)/(YMAX-ymid);
			}
			v[ke  ] = 	ax*ay;
			v[ke+1] = 10*ax*ay;
		}
	} else if (ndim == 3) {	
		dz = (ZMAX-ZMIN)/mz;
		xmid = 0.5*(XMIN + XMAX);
		ymid = 0.5*(YMIN + YMAX);
		zmid = 0.5*(ZMIN + ZMAX);
		for (kg=0; kg<NG; kg++) {
			x = XMIN + dx*(xmap[kg] - 0.5);
			y = YMIN + dy*(ymap[kg] - 0.5);
			z = ZMIN + dz*(zmap[kg] - 0.5);
			ke = kg*nvars;
			if (x < xmid) {
				ax = (x-XMIN)/(xmid-XMIN);
			} else {
				ax = (XMAX-x)/(XMAX-xmid);
			}
			if (y < ymid) {
				ay = (y-YMIN)/(ymid-YMIN);
			} else {
				ay = (YMAX-y)/(YMAX-ymid);
			}
			if (z < zmid) {
				az = (z-ZMIN)/(zmid-ZMIN);
			} else {
				az = (ZMAX-z)/(ZMAX-zmid);
			}
			v[ke  ] = 	20*ax*ay*az;
			v[ke+1] = 10*ax*ay*az;
		}	
	}
}

//------------------------------------------------------------------------------
// Rates of change of consituents are determined from constituent concentrations
// (later we may wish to add sources and sinks)
//------------------------------------------------------------------------------
void React(realtype *C, realtype *dCdt)
{
	dCdt[0] = K01*C[0]*C[1];
	dCdt[1] = -K01*C[0]*C[1];
}

//------------------------------------------------------------------------------
// dfdC[k'][k] = df[k']/dC[k]
// df(jc)/dC(ic) = dfdC[k], index k = jc + ic*nvars
//------------------------------------------------------------------------------
void ReactJac(realtype C[], realtype dfdC[], int nvars)
{
	int ic, jc, k;

	//dfdC[0][0] = K01*C[1];
	//dfdC[0][1] = K01*C[0];
	//dfdC[1][0] = -K01*C[1];
	//dfdC[1][1] = -K01*C[0];
	jc = 0; ic = 0; k = jc + ic*nvars;
	dfdC[k] = K01*C[1];
	jc = 0; ic = 1; k = jc + ic*nvars;
	dfdC[k] = K01*C[0];
	jc = 1; ic = 0; k = jc + ic*nvars;
	dfdC[k] = -K01*C[1];
	jc = 1; ic = 1; k = jc + ic*nvars;
	dfdC[k] = -K01*C[0];
}


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main(void)
{
	clock_t begin, end;
	double time_spent;
	int ndim, mx, my, mz, NG, nvars, nout;
	int *xmap, *ymap, *zmap;
	realtype *v;
	realtype dx, cdiff[3], vbnd[2], t0, t1, dtout;
	realtype CDIFF = 1.0e-1;		 /* coeff of diffusion */

  	begin = clock();

	ndim = 3;
	mx = 30;
	my = 30;
	mz = 30;
	nvars = 2;
	vbnd[0] = 0;
	vbnd[1] = 0;
	dx = (XMAX-XMIN)/mx;  
	cdiff[0] = CDIFF/(dx*dx);
	cdiff[1] = CDIFF/(dx*dx);
	cdiff[2] = CDIFF/(dx*dx);

	t0 = 0;
	t1 = 10;
	dtout = t1;
	nout = 20;
	MakeMaps(ndim,mx,my,mz,&NG,&xmap,&ymap,&zmap);

	v = (realtype *)malloc(NG*nvars*sizeof(realtype));
	SetIC(ndim,mx,my,mz,NG,nvars,xmap,ymap,zmap,v);
	Solve(ndim,mx,my,mz,NG,nvars,vbnd,cdiff,xmap,ymap,zmap,v,t0,t1,dtout,nout);

    end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Solve time: %6.1f\n",time_spent);
	return(0);
}
