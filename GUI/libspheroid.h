#ifndef LIBSPHEROID_H
#define LIBSPHEROID_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void get_dll_build_version(char *, int *);
void execute(int *,char *, int *,char *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *, int *, int *, int *, double *, int *, int *, bool *, double *, double *);
void get_scene(int *, int *);
void get_summary(int *, int *, int *);
void get_concdata(int *, int *, double *, double *);
//void get_distdata(int *, DIST_SET *);
void get_volprob(int *, double *, double *, double*);
void get_oxyprob(int *, double *, double *, double *);
void get_nfacs(int *);
void get_facs(double *);
void get_histo(int, double *, double *, double *, double *, double *, double *);

void get_constituents(int *, int *, int *, char *, int *);

//void get_profile_oxygen(double *, double *, int *);

void get_string(char **);

//
//
#ifdef __cplusplus
}
#endif

#endif // LIBSPHEROID2_H
