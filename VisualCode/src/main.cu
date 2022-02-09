#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "ouster_reconstruction.h"
//CUDA libraries
#define threadsPerBlock 16  //
#define numBlocks (1024/threadsPerBlock) //(n_AZBLK/1024)
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#ifndef __CUDACC__  
	#define __CUDACC__
	#include <device_functions.h>
#endif

/**
 * \brief Devices functions. 
 */
__device__ void mult_matrix_dev(double* A,unsigned int m,unsigned int n,double*B,unsigned int l,double*C) {
    for (unsigned int i = 0; i < l;i++) {
        for (unsigned int j = 0; j < m; j++) {
            C[j*l+i]=0;
            for (unsigned int k = 0; k < n; k++) {
                C[j * l +i] += A[j*n+k]*B[k*l + i];
            }
        }
    }
}
__device__ void rot_x_axis_dev(double* XYZ_points,double angle ){
    double temp[3], rot_matrix[9] = { 1,0,0,0,cos(angle),-sin(angle),0,sin(angle),cos(angle) };
    temp[0]=XYZ_points[0];
    temp[1]=XYZ_points[1];
    temp[2]=XYZ_points[2];
    mult_matrix_dev(rot_matrix, 3, 3, temp, 1, XYZ_points);
}

__device__ double eq_line_dev(double m,double x,double xb,double yb) {
    double y= m*(x-xb)+yb;
    return y;
}
/**
 * \brief GenerateSphereCUDA. genera esfera 
 */
__global__ void GenerateSphereCUDA(double* Point_Cloud){
    double beam_altitude_angles[n_beams]= {15.379*D180_MPI,13.236*D180_MPI,11.128*D180_MPI,9.03*D180_MPI,6.941*D180_MPI,4.878*D180_MPI,2.788*D180_MPI,0.705*D180_MPI,-1.454*D180_MPI,-3.448*D180_MPI,-5.518*D180_MPI,-7.601*D180_MPI,-9.697*D180_MPI,-11.789*D180_MPI,-13.914*D180_MPI,-16.062*D180_MPI};
    double beam_azimuth_angles[n_beams] = { -1.24*D180_MPI, -1.2145*D180_MPI, -1.1889*D180_MPI, -1.1634*D180_MPI, -1.1379*D180_MPI, -1.1123*D180_MPI, -1.0868*D180_MPI, -1.0613*D180_MPI, -1.0357*D180_MPI, -1.0102*D180_MPI, -0.98467*D180_MPI, -0.95913*D180_MPI, -0.9336*D180_MPI, -0.90807*D180_MPI, -0.88253*D180_MPI, -0.857*D180_MPI };
    //Obtenemos el ID del thread
    int thid = threadIdx.x + blockIdx.x * blockDim.x;//thid value from 0 to 1023
	//Generamos el primer azimuth 
    if (thid<n_beams){
        //ubicamos el punto del azimut referencial en el plano XZ
        Point_Cloud[3 * thid + 0] = Radius_sphere* cos(beam_altitude_angles[thid] - MPI_2); //x
        Point_Cloud[3 * thid + 1] = 0;                                         //y
        Point_Cloud[3 * thid + 2] = Radius_sphere * sin(beam_altitude_angles[thid] - MPI_2);//z
        //Realizamos la rotacion del punto con respecto al eje x debido al desfase
        rot_x_axis_dev(&Point_Cloud[3*thid],beam_azimuth_angles[thid]);
        /*Creamos los azimuts que inician en cada sector*/
        /*mirror points from quarter Donut*/
        Point_Cloud[(thid + n_AZBLK / 4 * n_beams)*3+0] = Point_Cloud[3 * thid + 0];
        Point_Cloud[(thid+ n_AZBLK / 4 * n_beams)*3+1] = Point_Cloud[3 * thid + 2];
        Point_Cloud[(thid + n_AZBLK / 4 * n_beams)*3+2] = -Point_Cloud[3 * thid + 1];
        /*mirror points from midle Donut*/
        Point_Cloud[(thid + n_AZBLK / 2 * n_beams)*3+0] = Point_Cloud[3 * thid + 0];
        Point_Cloud[(thid + n_AZBLK / 2 * n_beams)*3+1] = -Point_Cloud[3 * thid + 1];
        Point_Cloud[(thid + n_AZBLK / 2 * n_beams)*3+2] = -Point_Cloud[3 * thid + 2];
        /*mirror points from 3 quater Donut*/
        Point_Cloud[(thid + n_AZBLK * 3 / 4 * n_beams)*3 +0] = Point_Cloud[3 * thid + 0];
        Point_Cloud[(thid + n_AZBLK * 3 / 4 * n_beams)*3 +1] = -Point_Cloud[3 * thid + 2];
        Point_Cloud[(thid + n_AZBLK * 3 / 4 * n_beams)*3 +2] = Point_Cloud[3 * thid + 1];
    }
    __syncthreads();
    /*Definimos la matrix de rotación para los azimuth*/
    /*Procedemos a realizar el barrido de cada sector para obtener la donnut referencial*/
    double temp[3];
    if(thid>0){
        double rot_matrix[9] = { 1,0,0,0,cos(angle_between_azimuths*(double)thid),-sin(angle_between_azimuths*(double)thid),0,sin(angle_between_azimuths*(double)thid),cos(angle_between_azimuths*(double)thid) };
        double XYZ[3];
        for (int j = 0; j < n_beams;j++) {
            //Obtain point from referencial azimuth
            XYZ[0] = Point_Cloud[j*3 + 0];
            XYZ[1] = Point_Cloud[j*3 + 1];
            XYZ[2] = Point_Cloud[j*3 + 2];
            //rotate that point
            mult_matrix_dev (rot_matrix,3,3,XYZ,1,temp);
            //Set the new azimuth
            Point_Cloud[(thid * n_beams + j)*3 + 0] = temp[0];
            Point_Cloud[(thid * n_beams + j)*3 + 1] = temp[1];
            Point_Cloud[(thid * n_beams + j)*3 + 2] = temp[2];
        }
    }
    __syncthreads();
    //Procedemos a rotar la Donut generada
    for (unsigned int i = 1; i < n_donuts; i++) {
        double rot_motor_matrix[9]={ cos(rot_angle*(double)i),-sin(rot_angle*(double)i),0,sin(rot_angle*(double)i),cos(rot_angle*(double)i) ,0,0,0,1 };
        //multiplicamos a todos los n_point_perdonut de la Donut referencial con la matriz de rotación respectiva
        for(int j=0;j<n_beams;j++){
            mult_matrix_dev(rot_motor_matrix, 3, 3, &Point_Cloud[(thid*n_beams+j) * 3], 1, temp);
            Point_Cloud[(i * n_points_perDonut + thid*n_beams+j) * 3 + 0] = temp[0];
            Point_Cloud[(i * n_points_perDonut + thid*n_beams+j) * 3 + 1] = temp[1];
            Point_Cloud[(i * n_points_perDonut + thid*n_beams+j) * 3 + 2] = temp[2];        
        }
    }
    __syncthreads();
}
/**
 * \brief SupressOverlapCUDA. OR 
 */
__global__ void SupressOverlapCUDA(double* Point_Cloud){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;
    //Declare temporary variables
    double x,y,y_temp,x1,y1,m;//xn,yn
    //Set vertical limits
    double X_L1=Point_Cloud[0],X_Ln=Point_Cloud[(n_beams-1)*3+0];
    //Array which contains lineś parameters
    double L[(n_donuts-2)*5];
    //Declare parameters from Donut_2 to Donut_n-1. Ya que para la supresión solo
    //se necesita los parametros de la Donut anterior y no necesitamos el de la
    //última Donut
    for (int i = 1; i < n_donuts-1; i++)
    {
        //Hallamos dos puntos de la recta proyectada en el plano XY
        x1=Point_Cloud[i*n_points_perDonut*3+0];
        y1=Point_Cloud[i*n_points_perDonut*3+1];
        x=Point_Cloud[i*n_points_perDonut*3+n_beams*3+0];
        y=Point_Cloud[i*n_points_perDonut*3+n_beams*3+1];
        //Calculamos la pendiente
        m=(y-y1)/(x-x1);
        //Guardamos los valores de un punto de la recta del beam_0
        L[(i-1)*5+0]=x1;
        L[(i-1)*5+1]=y1;
        //Guardamos los valores de un punto de la recta del beam_n
        L[(i-1)*5+2]=Point_Cloud[i*n_points_perDonut*3+n_beams*3-3];//xn
        L[(i-1)*5+3]=Point_Cloud[i*n_points_perDonut*3+n_beams*3-3+1];//yn
        //Guardamos el valor de la pendiente hallada
        L[(i-1)*5+4]=m;
    }
    //Supress redundant data only for Donut 2
    int i=1;
    for (unsigned int j = 0; j < n_beams; j++){
        //Calculamos si la coordenada del punto x
        x=Point_Cloud[(i*n_points_perDonut+thid*n_beams+j)*3+0];
        //Analizamos si se encuentra en la zona de supresión
        if (X_Ln<=x&&x<=X_L1){
            Point_Cloud[(i*n_points_perDonut+thid*n_beams+j)*3]=0;
            Point_Cloud[(i*n_points_perDonut+thid*n_beams+j)*3+1]=0;
            Point_Cloud[(i*n_points_perDonut+thid*n_beams+j)*3+2]=0;
            //points_deleted=points_deleted+1;
        }
    }
    //Supress redundant for the rest of the Donuts
    //Creamos variable booleana para saber la zona del plano en donde se encuentra el punto
    bool left_side;
    unsigned int offset;
    for (unsigned int i = 2; i < n_donuts; i++){
        for (unsigned int j = 0; j < n_beams; j++){
            //Hallamos las coordenadas del punto a analizar
            offset=(i*n_points_perDonut+thid*n_beams+j)*3;
            x=Point_Cloud[offset];
            y=Point_Cloud[offset+1];
            //Evaluamos si se encuentra en la zona referencial
            if (X_Ln<=x){
                if(x<=X_L1){
                    Point_Cloud[offset]=0;
                    Point_Cloud[offset+1]=0;
                    Point_Cloud[offset+2]=0;
                    continue;
                }else{
                    //Se encuentra del lado derecho
                    left_side=false;
                }
            }else{
                //Se encuentra del lazo izquierdo
                left_side=true;
            }
            //Calculamos el valor de y_temp el cual limitará la zona
            y_temp=eq_line_dev(L[(i-2)*5+4],x,L[(i-2)*5+left_side*2],L[(i-2)*5+left_side*2+1]);
            //Le colocamos un signo negativo, o no, para poder realizar un único condicional para ambos casos
            y_temp=y_temp*(1.0-2*left_side);
            y=y*(1.0-2*left_side);
            //Evaluamos la condición de supresión
            if (y>=y_temp){
                //Eliminamos los puntos
                Point_Cloud[offset]=0;
                Point_Cloud[offset+1]=0;
                Point_Cloud[offset+2]=0;
                //points_deleted=points_deleted+1;
            }
        }   
    }
    __syncthreads();
}
/**
 * \brief OneDonutFillCUDA.  
 */
__global__ void ODF_part1(double* Point_Cloud,unsigned int* T,unsigned int* T_temp,unsigned int* count_array){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int v0,v1,v2;
    //Definimos los vértices
    //Realizamos la malla triangular para la Donut referencial
    for (unsigned int j = 0; j < n_beams-1; j++)
    {   
        v0=thid*n_beams+j;
        v2=v0+1;
        v1=(v0+n_beams+1)&mask;
        T[thid*(n_beams-1)*3*2+j*6]=v0;
        T[thid*(n_beams-1)*3*2+j*6+1]=v1;
        T[thid*(n_beams-1)*3*2+j*6+2]=v2;
        v2=v1;
        v1=v2-1;
        T[thid*(n_beams-1)*3*2+j*6+3]=v0;
        T[thid*(n_beams-1)*3*2+j*6+4]=v1;
        T[thid*(n_beams-1)*3*2+j*6+5]=v2;
    }
    //__syncthreads(); Not necessary
    //En base a la malla referencial hallamos las demás superficies
    double xp,yp,zp;
    unsigned int count=0,offset,temp_vex,n_triangles_perThreadandDonut=2*(n_beams-1);
    for (unsigned int i = 1; i < n_donuts; i++){
        //Analizamos cada vertice del tríangulo
        for (unsigned int j = 0; j < n_triangles_perThreadandDonut; j++){
             //Analizamos el punto del vertice v0
            offset=(thid*n_triangles_perThreadandDonut+j)*3;
            temp_vex=(T[offset]+i*n_points_perDonut);
            xp=Point_Cloud[temp_vex*3+0];
            yp=Point_Cloud[temp_vex*3+1];
            zp=Point_Cloud[temp_vex*3+2];
            if ((xp!=0)||(yp!=0)||(zp!=0)){
                //analizamos el punto del vertice v1
                temp_vex=(T[offset+1]+i*n_points_perDonut);
                xp=Point_Cloud[temp_vex*3+0];
                yp=Point_Cloud[temp_vex*3+1];
                zp=Point_Cloud[temp_vex*3+2];
                if ((xp!=0)||(yp!=0)||(zp!=0)){
                    //analizamos el punto del vertice v2
                    temp_vex=(T[offset+2]+i*n_points_perDonut);
                    xp=Point_Cloud[temp_vex*3+0];
                    yp=Point_Cloud[temp_vex*3+1];
                    zp=Point_Cloud[temp_vex*3+2];
                    if ((xp!=0)||(yp!=0)||(zp!=0)){
                        //Si todo lo anterior se cumple, guardamos el triángulo
                        T_temp[(n_donuts-1)*n_triangles_perThreadandDonut*3*thid+count*3+2]=temp_vex;
                        T_temp[(n_donuts-1)*n_triangles_perThreadandDonut*3*thid+count*3+1]=(T[offset+1]+i*n_points_perDonut);
                        T_temp[(n_donuts-1)*n_triangles_perThreadandDonut*3*thid+count*3]=(T[offset]+i*n_points_perDonut);
                        count++;
                    }
                }
            }
        }
    }
    count_array[thid]=count;
    __syncthreads();
}
#define ngpu 1024
#define NUM_BANKS 32
#define LOG_NUM_BANKS 5
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_NUM_BANKS)
__global__ void eScanGPU(unsigned int *g_odata, unsigned int *g_idata)
{
	__shared__ unsigned int temp[2*ngpu];// allocated on invocation
	int thid = threadIdx.x;
	int offset = 1,n=ngpu;
	
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

			unsigned int t = temp[ai];
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

__global__ void ODF_part2(unsigned int* T,unsigned int* T_temp,unsigned int* count_array,unsigned int* index_offset){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;   
    unsigned int count,offset,n_triangles_perDonut=n_AZBLK*2*(n_beams-1),n_triangles_perThreadandDonut=2*(n_beams-1);
    count=count_array[thid];
    //Ha este punto, cada hilo contiene una cierta cantidad de n triangulos que han de ser colocadas en el array original
    if(thid==0)
    offset=0;
    else
    offset=index_offset[thid]*3;
    //copy triangles 
    for(int z=0;z<count*3;z++){
        T[n_triangles_perDonut*3+offset+z]=T_temp[(n_donuts-1)*n_triangles_perThreadandDonut*3*thid+z];
    }
    __syncthreads();
}

/**
 * \brief Generate surface.  
 */
void Generate_surfaceGPU(double* Point_Cloud,unsigned int* T,unsigned int *pointer_n_triangles){
    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      CUDA        ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    cudaError_t cudaerr;
    double* Sphere_Cloud_dev;
    cudaMalloc((void**)(&Sphere_Cloud_dev), sizeof(double) * n_total_points * 3);
    //1st step. GenerateSphere
    GenerateSphereCUDA<<<4, 256 >>> (Sphere_Cloud_dev);
    //2nd step. Overlap removing
    SupressOverlapCUDA<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    double* Sphere_Cloud;
    Sphere_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    cudaerr=cudaMemcpy(Sphere_Cloud, Sphere_Cloud_dev, sizeof(double) *n_total_points * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to SphereCloud. CudaMalloc value=%i\n\r",cudaerr);
    unsigned int *OneDonutMesh_dev,*T_temp_dev,*count_array_dev,n_triangles_perDonut=n_AZBLK*2*(n_beams-1);
    cudaMalloc((void**)(&OneDonutMesh_dev), sizeof(unsigned int) * OneDonutFill_triangles * 3);
    cudaMalloc((void**)(&T_temp_dev), sizeof(unsigned int) * n_triangles_perDonut* (n_donuts-1) * 3);
    cudaMalloc((void**)(&count_array_dev), sizeof(unsigned int) * 1024);
    //3rd step. First part of the ODF
    ODF_part1<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev,OneDonutMesh_dev,T_temp_dev,count_array_dev);
    unsigned int *index_offset_array_dev;
    cudaMalloc((void**)(&index_offset_array_dev), sizeof(unsigned int) * 1024);
    //4th step. eScan GPU
    eScanGPU<<<1, 512 >>> (index_offset_array_dev,count_array_dev);//try with CPU. Try with masking. Chapter: toolkit pin-memory. Optimized method. Chapter: asynch transfer
    //5th step. last part ODF
    ODF_part2<<<numBlocks, threadsPerBlock >>> (OneDonutMesh_dev,T_temp_dev,count_array_dev,index_offset_array_dev);
    unsigned int *T_temp;
    T_temp=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    cudaerr=cudaMemcpy(T_temp,OneDonutMesh_dev, sizeof(unsigned int) *OneDonutFill_triangles * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to T_temp. CudaMalloc value=%i\n\r",cudaerr);
	/////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      CPU        ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //6th step. last Fill
    TwoandTri_Donut_Fill(Sphere_Cloud,&T_temp[OneDonutFill_triangles*3],&T_temp[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T_temp[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);
    cudaDeviceSynchronize();
    cudaFree(Sphere_Cloud_dev);
    cudaFree(OneDonutMesh_dev);
    cudaFree(T_temp_dev);
    cudaFree(count_array_dev);
    cudaFree(index_offset_array_dev);
   
    unsigned int temp_vex,n_triangles=0;
    double xp,yp,zp;
    for (unsigned int i = 0; i < n_total_triangles; i++){
        //Analizamos el punto del vertice v0
        temp_vex=T_temp[i*3];
        xp=Point_Cloud[temp_vex*3+0];
        yp=Point_Cloud[temp_vex*3+1];
        zp=Point_Cloud[temp_vex*3+2];
        if ((xp!=0)&&(yp!=0)&&(zp!=0)){
            //analizamos el punto del vertice v1
            temp_vex=T_temp[i*3+1];
            xp=Point_Cloud[temp_vex*3+0];
            yp=Point_Cloud[temp_vex*3+1];
            zp=Point_Cloud[temp_vex*3+2];
            if ((xp!=0)&&(yp!=0)&&(zp!=0)){
                //analizamos el punto del vertice v2
                temp_vex=T_temp[i*3+2];
                xp=Point_Cloud[temp_vex*3+0];
                yp=Point_Cloud[temp_vex*3+1];
                zp=Point_Cloud[temp_vex*3+2];
                if ((xp!=0)&&(yp!=0)&&(zp!=0)){
                    //Si todo lo anterior se cumple, guardamos el triángulo
                    T[n_triangles*3+2]=temp_vex;
                    T[n_triangles*3+1]=T_temp[i*3+1];
                    T[n_triangles*3]=T_temp[i*3];
                    n_triangles++;
                }
            }
        }
    }
    pointer_n_triangles[0]=n_triangles; 
    free(Sphere_Cloud);
    free(T_temp);
}

/******************************************************************/
/*************************       MAIN     *************************/
/******************************************************************/
#define TESTING 0
int main()
{
#if TESTING == 0
    #define iter 500.0
    /*Allocate memory*/
    double* Point_Cloud;
    Point_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    //Leemos del csv los datos reales
    FILE* archivo;
    archivo = fopen("../files/MinaData.csv", "r");
    char buffer[200];
    char* token;
    //Saltamos la primera línea
    fgets(buffer,sizeof(buffer),archivo);
    for (unsigned int i = 0; i < n_total_points; i++){
        fgets(buffer,sizeof(buffer),archivo);
        token = strtok(buffer,",");
        Point_Cloud[i*3+0]=atof(token);
        token = strtok(NULL,",");
        Point_Cloud[i*3+1]=atof(token);
        token = strtok(NULL,",\n");
        Point_Cloud[i*3+2]=atof(token);
    }    
    fclose(archivo);    
    
    /*****************************************************/
    /********************       CPU     ******************/
    /*****************************************************/
    unsigned int *T_cpu,n_triangles_real_data_cpu;
    T_cpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    clock_t startCPU;
	clock_t finishCPU;
    printf ("/********************solo CPU*********************/:\n");
	startCPU = clock();
    for (int i=0;i<iter;i++){
        Generate_surface(Point_Cloud,T_cpu,&n_triangles_real_data_cpu);
    }
	finishCPU = clock();
	printf("numero de triangulos: %d\n",n_triangles_real_data_cpu);
    printf("CPU: %fms\n", ((double)(finishCPU - startCPU))*1000/(double)CLOCKS_PER_SEC/iter);
	
    /*****************************************************/
    /********************       GPU     ******************/
    /*****************************************************/
    unsigned int *T_gpu,n_triangles_real_data_gpu;
    T_gpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    clock_t startGPU;
	clock_t finishGPU;
    printf ("/********************CPU y GPU*********************/:\n");
	startGPU = clock();
    for (int i=0;i<iter;i++){
        Generate_surfaceGPU(Point_Cloud,T_gpu,&n_triangles_real_data_gpu);
    }
	finishGPU = clock();
	printf("numero de triangulos: %d\n",n_triangles_real_data_gpu);
    printf("GPU: %fms\n", ((double)(finishGPU - startGPU))*1000/(double)CLOCKS_PER_SEC/iter);
	//creamos archivo para ver results
    archivo = fopen("../files/CUDAMesh.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < n_triangles_real_data_gpu; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_gpu[i*3+0], T_gpu[i * 3 + 1], T_gpu[i * 3 + 2]);
    }
    fclose(archivo);
    //------------------------------------------
	//----------Generate the DXF file-----------
	//------------------------------------------
	//Open the DXF file
    archivo = fopen("../files/CUDAMinaSurface.dxf", "w");
    //assert(archivo);
    //header
    fprintf(archivo, "0\nSECTION\n2\nENTITIES\n0\n");
	float x0,y0,z0,x1,y1,z1,x2,y2,z2;
	for (int i = 0; i < n_triangles_real_data_gpu; i++)
    {
        // get the coordinates of each point from the triangle
        x0 = Point_Cloud[T_gpu[i* 3+0]*3 + 0];
        y0 = Point_Cloud[T_gpu[i* 3+0]*3 + 1];
        z0 = Point_Cloud[T_gpu[i* 3+0]*3 + 2];
        
        x1 = Point_Cloud[T_gpu[i* 3+1]*3 + 0];
        y1 = Point_Cloud[T_gpu[i* 3+1]*3 + 1];
        z1 = Point_Cloud[T_gpu[i* 3+1]*3 + 2];
        
        x2 = Point_Cloud[T_gpu[i* 3+2]*3 + 0];
        y2 = Point_Cloud[T_gpu[i* 3+2]*3 + 1];
        z2 = Point_Cloud[T_gpu[i* 3+2]*3 + 2];
        //create new 3DFACE element
        fprintf(archivo, "3DFACE\n8\n1\n");
        fprintf(archivo, " 62\n %d\n", 142);//corresponding color of the autocad pallete
        fprintf(archivo, "10\n %.4f\n 20\n %.4f\n 30\n %.4f\n", x0, y0, z0);
        fprintf(archivo, "11\n %.4f\n 21\n %.4f\n 31\n %.4f\n", x1, y1, z1);
        fprintf(archivo, "12\n %.4f\n 22\n %.4f\n 32\n %.4f\n", x2, y2, z2);
        fprintf(archivo, "13\n %.4f\n 23\n %.4f\n 33\n %.4f\n", x2, y2, z2);
        fprintf(archivo, "0\n");
    }
    fprintf(archivo, "ENDSEC\n 0\nEOF\n");
    fclose(archivo);

    return 0;
#else
    #define iter 100.0
    /*************************************************************************************/
    /***********************************    CPU     **************************************/
    /*************************************************************************************/
    double* Point_Cloud;
    Point_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    unsigned int *T;
    T=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));

    clock_t start,stop;
    double timeCPU=0;
    printf ("/********************solo CPU*********************/:\n");
    start=clock();
    for (int z=0;z<iter;z++)
    Generate_sphere(Point_Cloud);
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("GS time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    
    start=clock();
    for (int z=0;z<iter;z++)
    Supress_redundant_data(Point_Cloud);
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("OR time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    
    start=clock();
    for (int z=0;z<iter;z++)
    One_Donut_Fill(Point_Cloud,T);
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("ODF time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    
    start=clock();
    for (int z=0;z<iter;z++)
    TwoandTri_Donut_Fill(Point_Cloud,&T[OneDonutFill_triangles*3],&T[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);    
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("LastFill time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    printf("total time: %fms\n",timeCPU);
        
    free(Point_Cloud);
    free(T);
    /*************************************************************************************/
    /***********************************    GPU     **************************************/
    /*************************************************************************************/    
    double* Point_Cloud_dev,*Point_Cloud_gpu;
    cudaMalloc((void**)(&Point_Cloud_dev), sizeof(double) * n_total_points * 3);
    Point_Cloud_gpu = (double*)malloc(n_total_points * 3 *sizeof(double));

    cudaError_t cudaerr;
    cudaEvent_t start_gpu, stop_gpu;
    float timeGPU,totalGPU=0;
    cudaEventCreate(&start_gpu);
    cudaEventCreate(&stop_gpu);

    printf ("/*********************CPU y GPU********************/\n");
    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++)
    GenerateSphereCUDA<<<4, 256 >>> (Point_Cloud_dev);
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    printf("GS time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;
    cudaDeviceSynchronize();
    
    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++)
    SupressOverlapCUDA<<<numBlocks, threadsPerBlock  >>> (Point_Cloud_dev);
    cudaerr=cudaMemcpy(Point_Cloud_gpu, Point_Cloud_dev, sizeof(double) *n_total_points * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to Point_Cloud_gpu. CudaMalloc value=%i\n\r",cudaerr);
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    printf("OR time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;

    unsigned int *T_gpu,*OneMesh_dev,*OneMesh_temp_dev,n_triangles_perDonut=n_AZBLK*2*(n_beams-1);
    unsigned int *count_array_dev,*index_offset_array_dev;
    
    cudaMalloc((void**)(&OneMesh_dev), sizeof(unsigned int) * n_total_triangles * 3);
    cudaMalloc((void**)(&OneMesh_temp_dev), sizeof(unsigned int) * n_triangles_perDonut* (n_donuts-1) * 3);
    cudaMalloc((void**)(&count_array_dev), sizeof(unsigned int) * 1024);
    cudaMalloc((void**)(&index_offset_array_dev), sizeof(unsigned int) * 1024);
    T_gpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));

    //3rd step. First part of the ODF
    cudaDeviceSynchronize();
    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++){
        ODF_part1<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev,OneMesh_dev,OneMesh_temp_dev,count_array_dev);
        cudaDeviceSynchronize();
        //4th step. eScan GPU
        eScanGPU<<<1, 1024/2 >>> (index_offset_array_dev,count_array_dev);
        cudaDeviceSynchronize();
        //5th step. last part ODF
        ODF_part2<<<numBlocks, threadsPerBlock >>> (OneMesh_dev,OneMesh_temp_dev,count_array_dev,index_offset_array_dev);
        cudaDeviceSynchronize();
        cudaerr=cudaMemcpy(T_gpu,OneMesh_dev, sizeof(unsigned int) *OneDonutFill_triangles * 3, cudaMemcpyDeviceToHost);
    
    }
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    if (cudaerr != 0)	printf("ERROR copying to T_gpu. CudaMalloc value=%i\n\r",cudaerr);
    printf("ODF time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;
  
    start=clock();
    for (int z=0;z<iter;z++)
    TwoandTri_Donut_Fill(Point_Cloud_gpu,&T_gpu[OneDonutFill_triangles*3],&T_gpu[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T_gpu[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);    
    stop=clock();
    totalGPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("LastFill time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    
    printf("total time: %fms\n",totalGPU);

    
    FILE* archivo;
    archivo = fopen("../testfiles/CUDASphere_cloud.csv", "w+");
    fprintf(archivo, "X, Y, Z\n");
    for (unsigned int i=0; i < n_total_points; i++) {
        fprintf(archivo,"%.4f, %.4f, %.4f\n", Point_Cloud_gpu[i*3+0], Point_Cloud_gpu[i * 3 + 1], Point_Cloud_gpu[i * 3 + 2]);
    }
    fclose(archivo);
    archivo = fopen("../testfiles/CUDAOneMesh.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < n_total_triangles; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_gpu[i*3+0], T_gpu[i * 3 + 1], T_gpu[i * 3 + 2]);
    }
    fclose(archivo);
    printf("fin\n");
    return;
  



    //Finalmente, liberamos el resto de memoria
    cudaFree(Point_Cloud_dev);
    cudaFree(OneMesh_dev);
    cudaFree(OneMesh_temp_dev);
    cudaFree(count_array_dev);
    cudaFree(index_offset_array_dev);
    free(Point_Cloud_gpu);
    free(T_gpu); 
    return 0;
#endif
}
