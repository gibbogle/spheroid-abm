#ifndef LIBBCELL32_H
#define LIBBCELL32_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void execute(int *,char *, int *,char *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *,int *,int *, int *, double *, int *, bool *, double *);
void get_scene(int *, int *);
void get_summary(int *, int *);
void get_concdata(int *, double *, double *);
void get_volprob(int *, double *, double *, double*);
void get_oxyprob(int *, double *, double*);
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBBCELL32_H
