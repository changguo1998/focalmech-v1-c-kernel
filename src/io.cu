#include <stdio.h>
#include <stdlib.h>

#include "io.h"

void read_input(const char*    fname,
                int*           n_obs_traces,
                int*           n_mt,
                int*           n_srcloc,
                Trace**        traces,
                MomentTensor** mts,
                GFtrace**      gf_database) {
    int   i, ngfs;
    FILE* fp = NULL;

    fp = fopen(fname, "rb");

    fread(n_obs_traces, sizeof(int), 1, fp);
    fread(n_mt, sizeof(int), 1, fp);
    fread(n_srcloc, sizeof(int), 1, fp);
    ngfs = (*n_obs_traces) * (*n_srcloc);

    printf("%d trace, %d mt, %d srcloc\n", *n_obs_traces, *n_mt, *n_srcloc);

    *traces = (Trace*)malloc(*n_obs_traces * sizeof(Trace));
    for(i = 0; i < *n_obs_traces; i++) {
        // printf("read trace [%d]\n", i);
        read_trace(&(*traces)[i], fp);
    }

    *mts = (MomentTensor*)malloc((*n_mt) * sizeof(MomentTensor));
    for(i = 0; i < *n_mt; i++) {
        // printf("read mt %d\n", i);
        read_mt(&(*mts)[i], fp);
    }

    *gf_database = (GFtrace*)malloc(ngfs * sizeof(GFtrace));
    for(i = 0; i < ngfs; i++) {
        // printf("read srcloc %d\n", i);
        read_gf_trace(&(*gf_database)[i], fp);
    }
    fclose(fp);
}

void write_result(const char*  fname,
                  const int    n_obs_traces,
                  const int    n_mt,
                  const int    n_srcloc,
                  const float* misfit_pl2,
                  const int*   misfit_pshift,
                  const float* misfit_sl2,
                  const int*   misfit_sshift,
                  const float* misfit_polarity,
                  const float* misfit_psr) {
    int   misfit_size;
    FILE* fp = NULL;
    fp = fopen(fname, "wb");
    misfit_size = n_obs_traces * n_mt * n_srcloc;
    fwrite(&n_obs_traces, sizeof(int), 1, fp);
    fwrite(&n_mt, sizeof(int), 1, fp);
    fwrite(&n_srcloc, sizeof(int), 1, fp);
    fwrite(misfit_pl2, sizeof(float), misfit_size, fp);
    fwrite(misfit_pshift, sizeof(int), misfit_size, fp);
    fwrite(misfit_sl2, sizeof(float), misfit_size, fp);
    fwrite(misfit_sshift, sizeof(int), misfit_size, fp);
    fwrite(misfit_polarity, sizeof(float), misfit_size, fp);
    fwrite(misfit_psr, sizeof(float), misfit_size, fp);
    fclose(fp);
}
