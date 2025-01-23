#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fm.h"

// #define DEBUG 1

#define allocmisfit(name, type) type *name = NULL; \
    name = (type *)malloc(n_result * sizeof(type)); \
    memset(name, 0, n_result * sizeof(type))

int main() {
    int64_t       n_obs_traces = 0, n_mt = 0, n_srcloc = 0;
    Trace*        traces = NULL;
    MomentTensor* mts = NULL;
    // green fun: n_srcloc x n_trace
    GFtrace* gf_database = NULL;

    printf("Load data\n");
    read_input("input.bin", &n_obs_traces, &n_mt, &n_srcloc, &traces, &mts, &gf_database);

    printf("Allocate memory\n");
    int64_t n_result = n_obs_traces * n_mt * n_srcloc;
    printf("n_result = %lld\n", n_result);

    // misfit: n_trace x n_mt x n_srcloc

    allocmisfit(misfit_pl2, float);
    allocmisfit(misfit_pshift, int64_t);
    allocmisfit(misfit_sl2, float);
    allocmisfit(misfit_sshift, int64_t);
    allocmisfit(misfit_polarity, float);
    allocmisfit(misfit_psr, float);

    int64_t imt, itr, isrc, igf, imisfit;
    printf("Calculate\n");
    for(itr = 0; itr < n_obs_traces; itr++)
        for(isrc = 0; isrc < n_srcloc; isrc++)
            for(imt = 0; imt < n_mt; imt++) {
                igf = isrc + n_srcloc * itr;
                imisfit = itr + n_obs_traces * (imt + n_mt * isrc);
#ifdef DEBUG
                    printf("trace: %lld, src: %lld, mt: %lld, gf: %lld, result: %lld\n", itr, isrc, imt, igf, imisfit);
#endif
                kernel(traces[itr], mts[imt], gf_database[igf],
                    &(misfit_pl2[imisfit]),
                    &(misfit_pshift[imisfit]),
                    &(misfit_sl2[imisfit]),
                    &(misfit_sshift[imisfit]),
                    &(misfit_polarity[imisfit]),
                    &(misfit_psr[imisfit]));
            }

    printf("Output\n");
    write_result("output.bin", n_obs_traces, n_mt, n_srcloc,
        misfit_pl2, misfit_pshift,
        misfit_sl2, misfit_sshift,
        misfit_polarity,
        misfit_psr);

    printf("Post process\n");
    for(itr = 0; itr < n_obs_traces; itr++)
        free_trace(&(traces[itr]));
    free(traces);
    free(mts);
    for(igf = 0; igf < n_srcloc; igf++)
        free_gf_trace(&(gf_database[igf]));
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
