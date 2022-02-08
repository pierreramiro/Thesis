// ING618 Algoritmos Paralelos - Prefix Sum
// Kernel code from Mark Harris - Parallel Prefix Sum (Scan) with CUDA
// Optimal in W, limited to sizes up to 2048 elements.

// Scan03: Adding proper offset on shared memory to avoid bank conflicts.

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h> // for rand();
#include <time.h> // to use clock() functions

#define ngpu 1024
#define threadsPerBlock 512 // threads per block
#define iter 10000
#define iterCPU 100000

#define NUM_BANKS 32
#define LOG_NUM_BANKS 5
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_NUM_BANKS)

void iScan( float* output, float* input, int length);
void eScan( float* output, float* input, int length);

__global__ void eScanGPU(float *g_odata, float *g_idata, int n)
{
	__shared__ float temp[2*ngpu];// allocated on invocation
	int thid = threadIdx.x;
	int offset = 1;
	
	int ai = thid;
	int bi = thid + (n/2);
	int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
	int bankOffsetB = CONFLICT_FREE_OFFSET(bi);
	temp[ai + bankOffsetA] = g_idata[ai];
	temp[bi + bankOffsetB] = g_idata[bi];

//	temp[2*thid] = g_idata[2*thid]; // load input into shared memory
//	temp[2*thid+1] = g_idata[2*thid+1];
	
	for (int d = n>>1; d > 0; d >>= 1) // build sum in place up the tree
	{
		__syncthreads();
		if (thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;
			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);
//			int ai = offset*(2*thid+1)-1;
//			int bi = offset*(2*thid+2)-1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
	}
	if (thid==0) { temp[n - 1 + CONFLICT_FREE_OFFSET(n - 1)] = 0; }
	//	if (thid == 0) { temp[n - 1] = 0; } // clear the last element
	for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
		if (thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;
			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);
//			int ai = offset*(2*thid+1)-1;
//			int bi = offset*(2*thid+2)-1;

			float t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}
	__syncthreads();
	g_odata[ai] = temp[ai + bankOffsetA];
	g_odata[bi] = temp[bi + bankOffsetB];
//	g_odata[2*thid] = temp[2*thid]; // write results to device memory
//	g_odata[2*thid+1] = temp[2*thid+1];
}

int main()
{

	float *in,*outgpu,*outcpu;
	cudaError_t cudaerr;
	int z;

	in  = (float *)malloc(ngpu*sizeof(float)); // input data
	outgpu = (float *)malloc(ngpu*sizeof(float)); // output data
	outcpu = (float *)malloc(ngpu*sizeof(float)); // output data


	// Fill data
	for(z=0;z<ngpu;z++)
	{
		in[z]=rand()%8; // Numbers between 0 and 7
	}
	in[0]=0;

	printf("\n %f %f %f %f %f\n",in[0],in[1],in[2],in[3],in[4]);
	// Setup timing using cudaEvent
	cudaEvent_t start, stop;
	float gpu_time;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// GPU using naive kernel (based on bs)
	printf("\nGPU eScan\n");

	int numBlocks;
    float *d_src = NULL;
    float *d_dst = NULL;

	cudaMalloc((void **)(&d_src), sizeof(float) * ngpu); // Input data
	// Move padded input image from Host to Device
	cudaerr = cudaMemcpy(d_src,in,sizeof(float)*ngpu,cudaMemcpyHostToDevice);
	if (cudaerr!=0)	printf("ERROR copying in data to d_src (Host to Dev). CudaMalloc value=%i\n\r", cudaerr);
    cudaMalloc((void **)(&d_dst), sizeof(float) * ngpu); // Output data

	cudaFuncSetCacheConfig(eScanGPU, cudaFuncCachePreferL1);

	cudaEventRecord(start);

	// Launch kernels
	numBlocks = 1; // WARNING: This must be an integer!!!, if not, add more code

	//Warmup
	eScanGPU <<<numBlocks,threadsPerBlock>>>(d_dst,d_src,ngpu);

	for(z=0;z<iter;z++)
	{
		eScanGPU <<<numBlocks,threadsPerBlock>>>(d_dst,d_src,ngpu);
	}

	cudaEventRecord( stop );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &gpu_time, start, stop );
	printf("GPU Time:  %fms\n\r",gpu_time/iter);
	printf("eScanGPU Speed: %f MegaOps/s\n",ngpu/(gpu_time/iter)/1000);
	cudaerr = cudaMemcpy(outgpu,d_dst,sizeof(float)*ngpu,cudaMemcpyDeviceToHost);
	if (cudaerr!=0)	printf("ERROR copying d_dst to outgpu (Dev to Host). CudaMalloc value=%i\n\r", cudaerr);

	clock_t startCPU;
	clock_t finishCPU;

	printf("\nCPU using eScan:\n");
	startCPU = clock();
	for(z=0;z<iterCPU;z++)
	{
		eScan(outcpu,in,ngpu);
	}
	finishCPU = clock();
	printf("CPU serial: %fms\n", (double)(finishCPU - startCPU)/1000/iterCPU);/// CLK_TCK);
	printf("eScanCPU Speed: %f MegaOps/s\n",ngpu/((double)(finishCPU - startCPU))/1000*iterCPU);
	// verify gpu vs cpu results
	for (z=0;z<ngpu;z++)
	{
		if (outgpu[z] != outcpu[z])
		{
			//error += abs(filteredImage[z] - filteredImageSerial[z]);
			printf("ERROR between CPU and GPU Scan on index: %i\n",z);
			printf("CPU: %f %f,%f %f\n", outcpu[z],outcpu[z+1],outcpu[z+2],outcpu[z+3]);
			printf("GPU: %f %f,%f %f\n", outgpu[z],outgpu[z+1],outgpu[z+2],outgpu[z+3]);
		
		}
	}


	printf("\n All DONE, press any key to end");


    cudaFree(d_src);
    cudaFree(d_dst);
	free(in);
	free(outcpu);
	free(outgpu);

    return 0;
}

void iScan( float* output, float* input, int length)
{

	output[0] = input[0]; 
	for(int z = 1; z < length; ++z)
	{
		output[z] = input[z] + output[z-1];
	}
}

void eScan( float* output, float* input, int length)
{

	output[0] = 0;
	for(int z = 1; z < length; ++z)
	{
		output[z] = input[z-1] + output[z-1];
	}
}
