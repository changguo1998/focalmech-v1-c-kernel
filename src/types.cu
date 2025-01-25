#include <stdio.h>
#include <stdlib.h>
#include "types.h"

// #define DEBUG 1
// #define GPU

#define readint(var) fread(&(var), sizeof(int), 1, fp)
#define writeint(var) fwrite(&(var), sizeof(int), 1, fp)
#define readfloat(var) fread(&(var), sizeof(float), 1, fp)
#define writefloat(var) fwrite(&(var), sizeof(float), 1, fp)

void printfloatvec(int n, float* v) {
    int i;
    if(n == 0) return;
    if(n < 7) {
        for(i = 0; i < n - 1; i++) printf("%g, ", v[i]);
        printf("%g", v[n - 1]);
        return;
    }
    printf("[");
    for(i = 0; i < 3; i++) printf("%g, ", v[i]);
    printf("...");
    for(i = n - 3; i < n; i++) printf(", %g", v[i]);
    printf("]");
}

void read_trace(Trace* tr, FILE* fp) {
    readint(tr->n_samples);
    readint(tr->ref_i);
    readint(tr->pwin_i);
    readint(tr->pwin_l);
    readint(tr->swin_i);
    readint(tr->swin_l);
#ifdef DEBUG
    printf("(read_trace)tr: %p, n: %d, ref: %d, pwin_i: %d, pwin_l: %d, swin_i: %d, swin_l: %d, %p, ",
        tr, tr->n_samples, tr->ref_i, tr->pwin_i, tr->pwin_l, tr->swin_i, tr->swin_l, tr->obs);
#endif
    tr->obs = (float*)malloc(tr->n_samples * sizeof(float));
    fread(tr->obs, sizeof(float), tr->n_samples, fp);
#ifdef DEBUG
    printf("%p\n", tr->obs);
    printfloatvec(tr->n_samples, tr->obs);
    printf("\n");
#endif
}

void write_trace(const Trace* tr, FILE* fp) {
    writeint(tr->n_samples);
    writeint(tr->ref_i);
    writeint(tr->pwin_i);
    writeint(tr->pwin_l);
    writeint(tr->swin_i);
    writeint(tr->swin_l);
    fwrite(tr->obs, sizeof(float), tr->n_samples, fp);
}

void free_trace(Trace* tr) {
    free(tr->obs);
    tr->obs = NULL;
}

void copy_trace_to_gpu(Trace *tr_mem, Trace *tr_gpu){
    // printf("(copy_trace_to_gpu) \n");
    // printf("(copy_trace_to_gpu) memcpy Trace struct\n");
    cudaMemcpy(tr_gpu, tr_mem, sizeof(Trace), cudaMemcpyHostToDevice);
    size_t nb = tr_mem->n_samples * sizeof(float);
    float *p;
    // printf("(copy_trace_to_gpu) alloc obs array on GPU\n");
    cudaMalloc((float**)&p, nb);
    // printf("(copy_trace_to_gpu) memcpy obs array\n");
    cudaMemcpy(p, tr_mem->obs, nb, cudaMemcpyHostToDevice);
    cudaMemcpy(&(tr_gpu->obs), &p, sizeof(float*), cudaMemcpyHostToDevice);
}

void free_trace_on_gpu(Trace* tr_gpu){
    cudaFree(tr_gpu->obs);
}

void read_mt(MomentTensor* mt, FILE* fp) {
    readfloat(mt->m11);
    readfloat(mt->m22);
    readfloat(mt->m33);
    readfloat(mt->m12);
    readfloat(mt->m13);
    readfloat(mt->m23);
#ifdef DEBUG
    printf("%g %g %g %g %g %g\n", mt->m11, mt->m22, mt->m33, mt->m12, mt->m13, mt->m23);
#endif
}

void write_mt(const MomentTensor* mt, FILE* fp) {
    writefloat(mt->m11);
    writefloat(mt->m22);
    writefloat(mt->m33);
    writefloat(mt->m12);
    writefloat(mt->m13);
    writefloat(mt->m23);
}

void read_gf_trace(GFtrace* gf, FILE* fp) {
    readint(gf->n_samples);
    readint(gf->p_i);
    readint(gf->s_i);
#ifdef DEBUG
    printf("(read_gf_trace) npts: %d, pi: %d, si: %d\n", gf->n_samples, gf->p_i, gf->s_i);
#endif

    gf->g11 = (float*)malloc(gf->n_samples * sizeof(float));
    fread(gf->g11, sizeof(float), gf->n_samples, fp);
#ifdef DEBUG
    printf("(read_gf_trace) g11: ");
    printfloatvec(gf->n_samples, gf->g11);
    printf("\n");
#endif

    gf->g22 = (float*)malloc(gf->n_samples * sizeof(float));
    fread(gf->g22, sizeof(float), gf->n_samples, fp);
#ifdef DEBUG
    printf("(read_gf_trace) g22: ");
    printfloatvec(gf->n_samples, gf->g22);
    printf("\n");
#endif

    gf->g33 = (float*)malloc(gf->n_samples * sizeof(float));
    fread(gf->g33, sizeof(float), gf->n_samples, fp);
#ifdef DEBUG
    printf("(read_gf_trace) g33: ");
    printfloatvec(gf->n_samples, gf->g33);
    printf("\n");
#endif

    gf->g12 = (float*)malloc(gf->n_samples * sizeof(float));
    fread(gf->g12, sizeof(float), gf->n_samples, fp);
#ifdef DEBUG
    printf("(read_gf_trace) g12: ");
    printfloatvec(gf->n_samples, gf->g12);
    printf("\n");
#endif

    gf->g13 = (float*)malloc(gf->n_samples * sizeof(float));
    fread(gf->g13, sizeof(float), gf->n_samples, fp);
#ifdef DEBUG
    printf("(read_gf_trace) g13: ");
    printfloatvec(gf->n_samples, gf->g13);
    printf("\n");
#endif

    gf->g23 = (float*)malloc(gf->n_samples * sizeof(float));
    fread(gf->g23, sizeof(float), gf->n_samples, fp);
#ifdef DEBUG
    printf("(read_gf_trace) g23: ");
    printfloatvec(gf->n_samples, gf->g23);
    printf("\n");
#endif
}

void write_gf_trace(const GFtrace* gf, FILE* fp) {
    writeint(gf->n_samples);
    writeint(gf->p_i);
    writeint(gf->s_i);
    fwrite(gf->g11, sizeof(float), gf->n_samples, fp);
    fwrite(gf->g22, sizeof(float), gf->n_samples, fp);
    fwrite(gf->g33, sizeof(float), gf->n_samples, fp);
    fwrite(gf->g12, sizeof(float), gf->n_samples, fp);
    fwrite(gf->g13, sizeof(float), gf->n_samples, fp);
    fwrite(gf->g23, sizeof(float), gf->n_samples, fp);
}

void free_gf_trace(GFtrace* gf) {
    free(gf->g11);
    gf->g11 = NULL;
    free(gf->g22);
    gf->g22 = NULL;
    free(gf->g33);
    gf->g33 = NULL;
    free(gf->g12);
    gf->g12 = NULL;
    free(gf->g13);
    gf->g13 = NULL;
    free(gf->g23);
    gf->g23 = NULL;
}

#define copygfcomponent(c)\
cudaMalloc((float**)&p, nb);\
cudaMemcpy(p, gf_mem->g##c, nb, cudaMemcpyHostToDevice);\
cudaMemcpy(&(gf_gpu->g##c), &p, sizeof(float*), cudaMemcpyHostToDevice)

void copy_gftrace_to_gpu(GFtrace *gf_mem, GFtrace *gf_gpu){
    cudaMemcpy(gf_gpu, gf_mem, sizeof(GFtrace), cudaMemcpyHostToDevice);
    size_t nb = gf_mem->n_samples * sizeof(float);
    float *p;
    copygfcomponent(11);
    copygfcomponent(22);
    copygfcomponent(33);
    copygfcomponent(12);
    copygfcomponent(13);
    copygfcomponent(23);
}

void free_gftrace_on_gpu(GFtrace *gf_gpu){
    cudaFree(gf_gpu->g11);
    cudaFree(gf_gpu->g22);
    cudaFree(gf_gpu->g33);
    cudaFree(gf_gpu->g12);
    cudaFree(gf_gpu->g13);
    cudaFree(gf_gpu->g23);
}

#undef DEBUG
