#ifndef TRANSFER_H
#define TRANSFER_H

#include <QMutex>

extern int showingVTK;

#define MAX_BC 100000
#define MAX_DC 500
#define MAX_BOND 1
#define NINFO 5

extern int VTKbuffer[100];
extern int BC_list[NINFO*MAX_BC];
extern int nBC_list;
extern int DC_list[NINFO*MAX_DC];
extern int nDC_list;
extern int bond_list[2*MAX_BOND];
extern int nbond_list;
extern QMutex mutex1, mutex2;

extern int summaryData[100];
extern int NX, NY, NZ;
extern int nt_vtk;
extern int istep;
extern bool leftb;

#endif // TRANSFER_H
