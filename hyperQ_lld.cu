#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#define N 2

__global__ void Kernel1()
{
    if(threadIdx.x == 0 && blockIdx.x == 0)
        printf("Hello from 1\n");
}

__global__ void Kernel2()
{
    if(threadIdx.x == 0 && blockIdx.x == 0)
        printf("Hello from 2\n");
}

int main (int argc, char **argv) {
    cudaStream_t streams[N];
    for(int i = 0; i < N; i++)
        cudaStreamCreate(&streams[i]);
    printf("Start\n");
    Kernel1<<<128, 256, 0, streams[0]>>>();
    Kernel2<<<128, 256, 0, streams[1]>>>();
    printf("Finish\n");
	return 0;
}
