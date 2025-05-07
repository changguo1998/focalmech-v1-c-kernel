#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>

#define LW 256

typedef float Float;

Float rrr(double m){
	return (double)rand() / (double)RAND_MAX * m;
}

Float gauss(double t, double mu, double sigma){
	return exp(-(t-mu)*(t-mu)/(sigma*sigma));
}

void calnm(Float* n, Float* nm, int i0, int l){
    // *n = nm[i0] * nm[i0];
	*n=0.0;
	int i;
	for(i=0; i<l; i++) *n += nm[i + i0]*nm[i + i0];
}

__device__
void calnm_gpu(Float* n, Float* nm, int i0, int l){
    // *n = nm[i0] * nm[i0];
	*n=0.0;
	int i;
	for(i=0; i<l; i++) *n += nm[i + i0]*nm[i + i0];
}

void cal_normalized_l2(Float* nm,
                       int    l,
                       int    i1,
                       Float* w1,
                       int    i2,
                       Float* g11,
                       Float* g22,
                       Float* g33,
                       Float* g12,
                       Float* g13,
                       Float* g23,
                       Float  m11,
                       Float  m22,
                       Float  m33,
                       Float  m12,
                       Float  m13,
                       Float  m23) {
    Float sum = 0, s1 = 0, s2 = 0, t = 0, g = 0;
    int   i;
    for(i = 0; i < l; i++) {
        s1 += w1[i1 + i] * w1[i1 + i];
        g = g11[i2 + i] * m11 +
            g22[i2 + i] * m22 +
            g33[i2 + i] * m33 +
            g12[i2 + i] * m12 +
            g13[i2 + i] * m13 +
            g23[i2 + i] * m23;
        s2 += g * g;
    }
    s1 = sqrtf(s1);
    s2 = sqrtf(s2);
    t = s1 * s2;
    for(i = 0; i < l; i++) {
        g = g11[i2 + i] * m11 +
            g22[i2 + i] * m22 +
            g33[i2 + i] * m33 +
            g12[i2 + i] * m12 +
            g13[i2 + i] * m13 +
            g23[i2 + i] * m23;
        sum += w1[i1 + i] * g;
    }
    *nm = sum / t;
}

__device__
void cal_normalized_l2_gpu(Float* nm,
                       int    l,
                       int    i1,
                       Float* w1,
                       int    i2,
                       Float* g11,
                       Float* g22,
                       Float* g33,
                       Float* g12,
                       Float* g13,
                       Float* g23,
                       Float  m11,
                       Float  m22,
                       Float  m33,
                       Float  m12,
                       Float  m13,
                       Float  m23) {
    Float sum = 0, s1 = 0, s2 = 0, t = 0, g = 0;
    int   i;
    for(i = 0; i < l; i++) {
        s1 += w1[i1 + i] * w1[i1 + i];
        g = g11[i2 + i] * m11 +
            g22[i2 + i] * m22 +
            g33[i2 + i] * m33 +
            g12[i2 + i] * m12 +
            g13[i2 + i] * m13 +
            g23[i2 + i] * m23;
        s2 += g * g;
    }
    s1 = sqrtf(s1);
    s2 = sqrtf(s2);
    t = s1 * s2;
    for(i = 0; i < l; i++) {
        g = g11[i2 + i] * m11 +
            g22[i2 + i] * m22 +
            g33[i2 + i] * m33 +
            g12[i2 + i] * m12 +
            g13[i2 + i] * m13 +
            g23[i2 + i] * m23;
        sum += w1[i1 + i] * g;
    }
    *nm = sum / t;
}

__global__
void knl(Float *wr, Float* syn, Float* res){
	int imisfit, bid, tid;
        bid = blockIdx.x;
        tid = threadIdx.x;
        imisfit = tid + bid * blockDim.x;

	if(imisfit >= (LW*2-1)) return;

	// cal_normalized_l2_gpu(&(res[imisfit]), LW,
	// 						imisfit, wr,
	// 						0, syn, syn, syn, syn, syn, syn,
	// 						   1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    calnm_gpu(&res[imisfit], wr, imisfit, LW);
}

int main(){
	int i, nbuf=LW*3;
	Float *ref, *ref_gpu, *syn, *syn_gpu, *buf, *buf_cpu, *buf_gpu, t;

	printf("allocate\n");
	ref = (Float *)malloc(nbuf*sizeof(Float));
    cudaMalloc((Float**)&ref_gpu, nbuf*sizeof(Float));

	syn = (Float *)malloc(LW*sizeof(Float));
    cudaMalloc((Float**)&syn_gpu, LW*sizeof(Float));

    buf = (Float *)malloc((LW*2-1)*sizeof(Float));
    buf_cpu = (Float *)malloc((LW*2-1)*sizeof(Float));
    cudaMalloc((Float**)&buf_gpu, (LW*2-1)*sizeof(Float));

	printf("initalize\n");
    for(i=0; i<nbuf; i++) ref[i] = rrr(0.01);

	for(i=0; i<LW; i++){
		t = gauss(i*(double)0.01, 1.28, 0.2);
        ref[i+LW] += t;
        syn[i] = t;
    }

	for(i=0; i<(LW*2-1); i++){
		buf[i] = 0.0;
        buf_cpu[i] = 0.0;
	}

	cudaMemcpy(ref_gpu, ref, nbuf*sizeof(Float), cudaMemcpyHostToDevice);
    cudaMemcpy(syn_gpu, syn, LW*sizeof(Float), cudaMemcpyHostToDevice);
	cudaMemcpy(buf_gpu, buf, (LW*2-1)*sizeof(Float), cudaMemcpyHostToDevice);

    // cpu
	printf("cpu\n");
	for(i=0; i<(LW*2-1); i++){
    	// cal_normalized_l2(&(buf_cpu[i]), LW,
		// 					i, ref,
		// 					0, syn, syn, syn, syn, syn, syn,
		// 					   1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        calnm(&buf_cpu[i], ref, i, LW);
	}

	// gpu
	printf("gpu\n");
	dim3 block(16);
    dim3 grid((LW*2-1)/16+1);
	knl<<<grid, block>>>(ref_gpu, syn_gpu, buf_gpu);
    cudaDeviceSynchronize();
	cudaMemcpy(buf, buf_gpu, (LW*2-1)*sizeof(Float), cudaMemcpyDeviceToHost);

	// compare
    FILE *fid=NULL;
	fid = fopen("test_nl2.txt", "w");
	printf("output\n");
	for(i=0; i<(LW*2-1); i++){
		t = buf_cpu[i] - buf[i];
		fprintf(fid, "%.16f - %.16f = %.16f\n", buf_cpu[i], buf[i], t);
	}
	fclose(fid);

	fid = fopen("test_nl2_raw.txt", "w");
	for(i=0; i<nbuf; i++) fprintf(fid, "%.16f\n", ref[i]);
	fprintf(fid, "==\n");
	for(i=0; i<LW; i++) fprintf(fid, "%.16f\n", syn[i]);
    fclose(fid);

	// free memory
	free(ref);
	free(syn);
    free(buf_cpu);
    free(buf);
    cudaFree(ref_gpu);
    cudaFree(syn_gpu);
    cudaFree(buf_gpu);
	return 0;
}
