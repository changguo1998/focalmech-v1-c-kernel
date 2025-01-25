#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cuda.h>

#include "fm.h"

// #define DEBUG 1

#define allocmisfit(name, type)                                               \
  type *name = NULL, *name##_gpu = NULL;                                      \
  name = (type *)malloc (n_result * sizeof (type));                           \
  memset(name, 0, n_result * sizeof(type));                                   \
  cudaMalloc((type**)&(name##_gpu), n_result * sizeof (type));                \
  cudaMemcpy(name##_gpu, name, n_result * sizeof(type), cudaMemcpyHostToDevice)

#define copybackmisfit(name, type) cudaMemcpy(name, name##_gpu, n_result * sizeof(type), cudaMemcpyDeviceToHost)
#define freemisfit(name) free(name); cudaFree(name##_gpu)

int main() {

    printf("Begin\n");
    int iDeviceCount = 0;
    cudaError_t error = cudaGetDeviceCount(&iDeviceCount);

    if (error != cudaSuccess || iDeviceCount == 0){
        printf("No GPU found\n");
        exit(-1);
    }

    int n_obs_traces = 0, n_mt = 0, n_srcloc = 0;
    // int imt, itr, isrc, imisfit, igf;
    int itr, igf;
    Trace* traces = NULL, *traces_gpu = NULL;
    MomentTensor* mts = NULL, *mts_gpu = NULL;
    // green fun: n_srcloc x n_trace
    GFtrace* gf_database = NULL, *gf_gpu = NULL;

    printf("Load data\n");
    read_input("input_cuda.bin", &n_obs_traces, &n_mt, &n_srcloc, &traces, &mts,
        &gf_database);
    printf("Copy to GPU\n");
    // printf("    traces\n");
    // printf("        alloc\n");
    cudaMalloc((Trace**)&(traces_gpu), n_obs_traces * sizeof(Trace));
    // printf("        allocted at %p\n", traces_gpu);
    for(itr = 0; itr < n_obs_traces; itr++) copy_trace_to_gpu(&(traces[itr]), &(traces_gpu[itr]));
    // printf("    mts\n");
    cudaMalloc((MomentTensor**)&(mts_gpu), n_mt * sizeof(MomentTensor));
    cudaMemcpy(mts_gpu, mts, n_mt*sizeof(MomentTensor), cudaMemcpyHostToDevice);
    // printf("    gfs\n");
    cudaMalloc((GFtrace**)&(gf_gpu), (n_obs_traces*n_srcloc)*sizeof(GFtrace));
    for(igf = 0; igf < n_obs_traces*n_srcloc; igf ++) copy_gftrace_to_gpu(&(gf_database[igf]), &(gf_gpu[igf]));

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

    printf("Calculate\n");
//    for(itr = 0; itr < n_obs_traces; itr++)
//        for(isrc = 0; isrc < n_srcloc; isrc++)
//            for(imt = 0; imt < n_mt; imt++) {
//                imisfit = itr + n_obs_traces * (imt + n_mt * isrc);
//                kernel(imisfit,
//                    n_obs_traces, traces,
//                    n_mt, mts,
//                    n_srcloc, gf_database,
//                    misfit_pl2, misfit_pshift,
//                    misfit_sl2, misfit_sshift,
//                    misfit_polarity, misfit_psr);
//            }
    dim3 block(32);
    dim3 grid(n_result/32+1);
    kernel<<<grid,block>>>(n_obs_traces, traces_gpu, n_mt, mts_gpu, n_srcloc, gf_gpu,
                    misfit_pl2_gpu, misfit_pshift_gpu,
                    misfit_sl2_gpu, misfit_sshift_gpu,
                    misfit_polarity_gpu, misfit_psr_gpu);
    cudaDeviceSynchronize();
    // copy misfit to memory
    copybackmisfit(misfit_pl2, float);
    copybackmisfit(misfit_pshift, int);
    copybackmisfit(misfit_sl2, float);
    copybackmisfit(misfit_sshift, int);
    copybackmisfit(misfit_polarity, float);
    copybackmisfit(misfit_psr, float);
    // cudaMemcpy(misfit_pl2_gpu, misfit_pl2, n_result*sizeof(float), cudaMemcpyDeviceToHost);
    // cudaMemcpy(misfit_pshift_gpu, misfit_pshift, n_result*sizeof(int), cudaMemcpyDeviceToHost);
    // cudaMemcpy(misfit_sl2_gpu, misfit_sl2, n_result*sizeof(float), cudaMemcpyDeviceToHost);
    // cudaMemcpy(misfit_sshift_gpu, misfit_sshift, n_result*sizeof(int), cudaMemcpyDeviceToHost);
    // cudaMemcpy(misfit_polarity_gpu, misfit_polarity, n_result*sizeof(float), cudaMemcpyDeviceToHost);
    // cudaMemcpy(misfit_psr_gpu, misfit_psr, n_result*sizeof(float), cudaMemcpyDeviceToHost);

    printf("Output\n");
    write_result("output_cuda.bin", n_obs_traces, n_mt, n_srcloc,
        misfit_pl2, misfit_pshift,
        misfit_sl2, misfit_sshift,
        misfit_polarity, misfit_psr);

    printf("Post process\n");
    // traces
    for(itr = 0; itr < n_obs_traces; itr++){
        free_trace(&(traces[itr]));
        free_trace_on_gpu(&(traces_gpu[itr]));
    }
    free(traces);
    cudaFree(traces_gpu);
    // mt
    free(mts);
    cudaFree(mts_gpu);
    // gf
    for(int igf = 0; igf < n_obs_traces * n_srcloc; igf++) {
        free_gf_trace(&(gf_database[igf]));
        free_gftrace_on_gpu(&(gf_gpu[igf]));
    }
    free(gf_database);
    cudaFree(gf_gpu);
    //
    freemisfit(misfit_pl2);
    freemisfit(misfit_pshift);
    freemisfit(misfit_sl2);
    freemisfit(misfit_sshift);
    freemisfit(misfit_polarity);
    freemisfit(misfit_psr);

    cudaDeviceReset();

    return 0;
}

#undef DEBUG
