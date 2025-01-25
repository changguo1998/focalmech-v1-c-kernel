#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "fm.h"

// #define DEBUG 1

#define NUM_THREADS 14

#define allocmisfit(name, type) type *name = NULL; \
    name = (type *)malloc(n_result * sizeof(type)); \
    memset(name, 0, n_result * sizeof(type))

int main() {
    int       n_obs_traces = 0, n_mt = 0, n_srcloc = 0;
    Trace*        traces = NULL;
    MomentTensor* mts = NULL;
    // green fun: n_srcloc x n_trace
    GFtrace* gf_database = NULL;

    omp_set_num_threads(NUM_THREADS);

    printf("Load data\n");
    read_input("input_omp.bin", &n_obs_traces, &n_mt, &n_srcloc, &traces, &mts, &gf_database);

    printf("Allocate memory\n");
    int n_result = n_obs_traces * n_mt * n_srcloc;
    printf("n_result = %d\n", n_result);

    // misfit: n_trace x n_mt x n_srcloc

    allocmisfit(misfit_pl2, float);
    allocmisfit(misfit_pshift, int);
    allocmisfit(misfit_sl2, float);
    allocmisfit(misfit_sshift, int);
    allocmisfit(misfit_polarity, float);
    allocmisfit(misfit_psr, float);

    int imt, itr, isrc, imisfit;
    printf("Calculate\n");

#pragma omp parallel for default(none) \
    private(imt, itr, isrc, imisfit) \
    shared(n_obs_traces, n_mt, n_srcloc, traces, mts, gf_database, misfit_pl2, \
    misfit_pshift, misfit_sl2, misfit_sshift, misfit_polarity, misfit_psr) \
    collapse(3)
    for(itr = 0; itr < n_obs_traces; itr++)
        for(isrc = 0; isrc < n_srcloc; isrc++)
            for(imt = 0; imt < n_mt; imt++) {
                imisfit = itr + n_obs_traces * (imt + n_mt * isrc);
#ifdef DEBUG
                printf("trace: %d, src: %d, mt: %d, gf: %d, result: %d\n", itr, isrc, imt, igf, imisfit);
#endif
                kernel(imisfit,
                    n_obs_traces, traces,
                    n_mt, mts,
                    n_srcloc, gf_database,
                    misfit_pl2, misfit_pshift,
                    misfit_sl2, misfit_sshift,
                    misfit_polarity, misfit_psr);
            }


    printf("Output\n");
    write_result("output_omp.bin", n_obs_traces, n_mt, n_srcloc,
        misfit_pl2, misfit_pshift,
        misfit_sl2, misfit_sshift,
        misfit_polarity,
        misfit_psr);

    printf("Post process\n");
    for(itr = 0; itr < n_obs_traces; itr++)
        free_trace(&(traces[itr]));
    free(traces);
    free(mts);
    for(int i = 0; i <  n_obs_traces*n_srcloc; i++)
        free_gf_trace(&(gf_database[i]));
    free(gf_database);
    free(misfit_pl2);
    free(misfit_pshift);
    free(misfit_sl2);
    free(misfit_sshift);
    free(misfit_polarity);
    free(misfit_psr);

    return 0;
}

#undef DEBUG
