#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "types.h"
#include "io.h"

void read_input(char*          fname,
                int64_t*       n_obs_traces,
                int64_t*       n_mt,
                int64_t*       n_srcloc,
                Trace**        traces,
                MomentTensor** mts,
                GFtrace**      gf_database) {
    int64_t i, ngfs;
    FILE*   fp = NULL;

    fp = fopen(fname, "rb");

    fread(n_obs_traces, sizeof(int64_t), 1, fp);
    fread(n_mt, sizeof(int64_t), 1, fp);
    fread(n_srcloc, sizeof(int64_t), 1, fp);
    ngfs = (*n_obs_traces) * (*n_srcloc);

    printf("%lld trace, %lld mt, %lld srcloc\n", *n_obs_traces, *n_mt, *n_srcloc);

    *traces = (Trace*)malloc(*n_obs_traces * sizeof(Trace));
    for(i = 0; i < *n_obs_traces; i++) {
        // printf("read trace [%lld]\n", i);
        read_trace(&(*traces)[i], fp);
    }

    *mts = (MomentTensor*)malloc((*n_mt) * sizeof(MomentTensor));
    for(i = 0; i < *n_mt; i++) {
        // printf("read mt %lld\n", i);
        read_mt(&(*mts)[i], fp);
    }

    *gf_database = (GFtrace*)malloc(ngfs * sizeof(GFtrace));
    for(i = 0; i < ngfs; i++) {
        // printf("read srcloc %lld\n", i);
        read_gf_trace(&(*gf_database)[i], fp);
    }
    fclose(fp);
}

void write_result(char*          fname,
                  const int64_t  n_obs_traces,
                  const int64_t  n_mt,
                  const int64_t  n_srcloc,
                  const float*   misfit_pl2,
                  const int64_t* misfit_pshift,
                  const float*   misfit_sl2,
                  const int64_t* misfit_sshift,
                  const float*   misfit_polarity,
                  const float*   misfit_psr) {
    int64_t misfit_size;
    FILE*   fp = NULL;
    fp = fopen(fname, "wb");
    misfit_size = n_obs_traces * n_mt * n_srcloc;
    fwrite(&n_obs_traces, sizeof(int64_t), 1, fp);
    fwrite(&n_mt, sizeof(int64_t), 1, fp);
    fwrite(&n_srcloc, sizeof(int64_t), 1, fp);
    fwrite(misfit_pl2, sizeof(float), misfit_size, fp);
    fwrite(misfit_pshift, sizeof(int64_t), misfit_size, fp);
    fwrite(misfit_sl2, sizeof(float), misfit_size, fp);
    fwrite(misfit_sshift, sizeof(int64_t), misfit_size, fp);
    fwrite(misfit_polarity, sizeof(float), misfit_size, fp);
    fwrite(misfit_psr, sizeof(float), misfit_size, fp);
    fclose(fp);
}
