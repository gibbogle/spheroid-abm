#ifndef GLOBAL_H
#define GLOBAL_H

#include <QtGui>

// Note that in the Fortran DLL the chemokine constituent numbering starts at 1
#define CFSE 0
#define OXYGEN 1
#define GLUCOSE 2
#define TRACER 3
#define DRUG_A 4
#define DRUG_A_METAB_1 5
#define DRUG_A_METAB_2 6
#define DRUG_B 7
#define DRUG_B_METAB_1 8
#define DRUG_B_METAB_2 9
#define GROWTH_RATE 10      // we pretend that this is a concentration
#define CELL_VOLUME 11      // we pretend that this is a concentration

#define MAX_CELLS 200000
#define N_CELLINFO 7
//#define N_FACS_VARS 3

namespace Global
{
    extern int data1;
    extern int data2;

    extern int MAX_CHEMO;
    extern int NX, NY, NZ;
    extern double DELTA_T;
    extern double dfraction;
    extern int nt_vtk;
    extern int istep;
    extern bool leftb;

    extern int nvars_used;
    extern int GUI_to_DLL_index[32];
    extern int DLL_to_GUI_index[32];
    extern QString var_string[32];

    extern double *FACS_data;
    extern int nFACS_cells;
    extern int nFACS_dim;

    extern double *histo_data;
    extern int nhisto_boxes;
    extern int nhisto_dim;
    extern double histo_vmax[3*32];
    extern int histo_celltype;

    extern int summaryData[100];
    extern int i_hypoxia_cutoff;
    extern int i_growth_cutoff;

    extern double concData[4000];
    extern int conc_nc;
    extern double conc_dx;

    extern double volProb[100];
    extern int vol_nv;
    extern double vol_v0;
    extern double vol_dv;
    extern double oxyProb[100];
    extern int oxy_nv;
    extern double oxy_dv;

    extern int cell_list[N_CELLINFO*MAX_CELLS];
    extern int ncell_list;

    extern bool showingVTK;
    extern bool recordingVTK;
    extern bool showingFACS;
    extern bool recordingFACS;

} // namespace Global

#endif // GLOBAL_H
