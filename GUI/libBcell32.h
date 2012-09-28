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
void get_dimensions(int *,int *,int *);
void get_scene(int *, int *, int *, int *, int *, int *);
void get_summary(int *);
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBBCELL32_H
