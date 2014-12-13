#include "global.h"

namespace Global
{
    int MAX_CHEMO;
    int NX, NY, NZ;
    double DELTA_T;
    double dfraction;
    int nt_vtk;
    int istep;
    bool leftb;

    int nvars_used;
    int GUI_to_DLL_index[32];
    int DLL_to_GUI_index[32];
    QString var_string[32];

    double *FACS_data=NULL;
    int nFACS_cells=0;
    int nFACS_dim=0;

    double *histo_data=NULL;
    int nhisto_boxes=20;
    int nhisto_dim=0;
    double histo_vmax[3*32];
    int histo_celltype=0;

    int summaryData[100];
    int i_hypoxia_cutoff;
    int i_growth_cutoff;

    double concData[4000];
    int conc_nc;
    double conc_dx;

    double volProb[100];
    int vol_nv;
    double vol_v0;
    double vol_dv;
    double oxyProb[100];
    int oxy_nv;
    double oxy_dv;

    int cell_list[N_CELLINFO*MAX_CELLS];
    int ncell_list;

//    double *profile_x[20];
//    double *profile_y[20];
//    int profile_n[20];

    bool showingVTK;
    bool recordingVTK;
    bool showingFACS;
    bool recordingFACS;


} // namespace Global
