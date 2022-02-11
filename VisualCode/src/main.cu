#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "ouster_reconstruction.h"
//CUDA libraries
#define threadsPerBlock 8  //
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

__global__ void cudaGenerateAZBLK(double* Point_Cloud){
    double beam_altitude_angles[n_beams]= {15.379*D180_MPI,13.236*D180_MPI,11.128*D180_MPI,9.03*D180_MPI,6.941*D180_MPI,4.878*D180_MPI,2.788*D180_MPI,0.705*D180_MPI,-1.454*D180_MPI,-3.448*D180_MPI,-5.518*D180_MPI,-7.601*D180_MPI,-9.697*D180_MPI,-11.789*D180_MPI,-13.914*D180_MPI,-16.062*D180_MPI};
    double beam_azimuth_angles[n_beams] = { -1.24*D180_MPI, -1.2145*D180_MPI, -1.1889*D180_MPI, -1.1634*D180_MPI, -1.1379*D180_MPI, -1.1123*D180_MPI, -1.0868*D180_MPI, -1.0613*D180_MPI, -1.0357*D180_MPI, -1.0102*D180_MPI, -0.98467*D180_MPI, -0.95913*D180_MPI, -0.9336*D180_MPI, -0.90807*D180_MPI, -0.88253*D180_MPI, -0.857*D180_MPI };
    int thid = threadIdx.x + blockIdx.x * blockDim.x;//thid value from 0 to 1023
	//Generamos el primer azimuth 
    if (thid<n_beams*3){
        int index_nBeams=thid/3;
        int index_XYZ=thid-index_nBeams*3;
        double XYZ[3];
        //ubicamos el punto del azimut referencial en el plano XZ
        XYZ[0] = Radius_sphere* cos(beam_altitude_angles[index_nBeams] - MPI_2); //x
        XYZ[1] = 0;                                                              //y
        XYZ[2] = Radius_sphere * sin(beam_altitude_angles[index_nBeams] - MPI_2);//z
        //Realizamos la rotacion del punto con respecto al eje x debido al desfase
        rot_x_axis_dev(XYZ,beam_azimuth_angles[index_nBeams]);
        Point_Cloud[3 * index_nBeams + 0] = XYZ[0];
        Point_Cloud[3 * index_nBeams + 1] = XYZ[1];
        Point_Cloud[3 * index_nBeams + 2] = XYZ[2];
        /*Creamos los azimuts que inician en cada sector*/
        if (index_XYZ==0){
            Point_Cloud[(index_nBeams + n_AZBLK / 4 * n_beams)*3+0] = XYZ[0];
            Point_Cloud[(index_nBeams + n_AZBLK / 2 * n_beams)*3+0] = XYZ[0];
            Point_Cloud[(index_nBeams + n_AZBLK * 3 / 4 * n_beams)*3 +0] = XYZ[0];
        }else if (index_XYZ==1){
            Point_Cloud[(index_nBeams+ n_AZBLK / 4 * n_beams)*3+1] = XYZ[2];
            Point_Cloud[(index_nBeams + n_AZBLK / 2 * n_beams)*3+1] = -XYZ[1];
            Point_Cloud[(index_nBeams + n_AZBLK * 3 / 4 * n_beams)*3 +1] = -XYZ[2];
        }else{
            Point_Cloud[(index_nBeams + n_AZBLK / 4 * n_beams)*3+2] = -XYZ[1];
            Point_Cloud[(index_nBeams + n_AZBLK / 2 * n_beams)*3+2] = -XYZ[2];
            Point_Cloud[(index_nBeams + n_AZBLK * 3 / 4 * n_beams)*3 +2] = XYZ[1];
        }
    }
}

__global__ void cudaGenerateDonut(double* Point_Cloud){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;//thid value from 0 to 1023
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
}

__global__ void cudaGenerateSphere(double* Point_Cloud){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;//thid value from 0 to 1023
    double temp[3];
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
__global__ void cudaODF_part1(double* Point_Cloud,unsigned int* T){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;
    if(thid<n_AZBLK*(n_beams-1)){
        int index_AZBLK=thid/(double)(n_beams-1);
        int index_nBeams=(thid-index_AZBLK*(n_beams-1));
        //T[thid*3]=thid;
        //T[thid*3+1]=index_AZBLK;
        //T[thid*3+2]=index_nBeams;
        
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6]=thid;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+1]=index_AZBLK;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+2]=index_nBeams;
        
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+3]=thid;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+4]=index_AZBLK;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+5]=index_nBeams;

    
    }else{
    //if(thid<n_AZBLK*(n_beams-1)){
        int index_AZBLK=thid/(double)(n_beams-1);
        int index_nBeams=(thid-index_AZBLK*(n_beams-1));
        unsigned int v0,v1,v2;
        //Definimos los vértices
        //Realizamos la malla triangular para la Donut referencial
        v0=index_AZBLK*n_beams+index_nBeams;
        v2=v0+1;
        v1=(v0+n_beams+1)&mask;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6]=v0;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+1]=v1;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+2]=v2;
        v2=v1;
        v1=v2-1;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+3]=v0;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+4]=v1;
        T[index_AZBLK*(n_beams-1)*3*2+index_nBeams*6+5]=v2;
    }
}

__global__ void cudaODF_part2(double* Point_Cloud,unsigned int* T,unsigned int* flag_array){
    int thid = threadIdx.x + blockIdx.x * blockDim.x;
    if (thid<n_triangles_perDonut*(n_beams-1)){ //thid from  0 to 184320
        int index_Donut=thid/(double)n_triangles_perDonut+1;
        int index_triPerDonut=thid-(index_Donut-1)*n_triangles_perDonut;
        double xp,yp,zp;
        unsigned int flag=0;
        unsigned int temp_vex=(T[index_triPerDonut*3]+index_Donut*n_points_perDonut);
        xp=Point_Cloud[temp_vex*3+0];
        yp=Point_Cloud[temp_vex*3+1];
        zp=Point_Cloud[temp_vex*3+2];
            
        if ((xp!=0)||(yp!=0)||(zp!=0)){
            //analizamos el punto del vertice v1
            temp_vex=(T[index_triPerDonut*3+1]+index_Donut*n_points_perDonut);
            xp=Point_Cloud[temp_vex*3+0];
            yp=Point_Cloud[temp_vex*3+1];
            zp=Point_Cloud[temp_vex*3+2];
            if ((xp!=0)||(yp!=0)||(zp!=0)){
                //analizamos el punto del vertice v2
                temp_vex=(T[index_triPerDonut*3+2]+index_Donut*n_points_perDonut);
                xp=Point_Cloud[temp_vex*3+0];
                yp=Point_Cloud[temp_vex*3+1];
                zp=Point_Cloud[temp_vex*3+2];
                if ((xp!=0)||(yp!=0)||(zp!=0)){
                    //Si todo lo anterior se cumple, guardamos el triángulo
                    flag=1;
                }
            }
        }
        flag_array[thid-n_triangles_perDonut]=flag;
    }  
}
void cpuODF_part3 (double* Point_Cloud,unsigned int *T,unsigned int *flag_array){
    //analizamos los flags
    int count=0;
    for (unsigned int z = 0; z < n_triangles_perDonut*(n_donuts-1); z++)
    {
        if (flag_array[z]==1){
            count++;
        }
    }
    printf ("%d\n",count);
    

}

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
    unsigned int count,offset,n_triangles_perThreadandDonut=2*(n_beams-1);
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
void Generate_surfaceGPU(double* Point_Cloud,unsigned int* T,double* Sphere_Cloud,unsigned int* T_Sphere,unsigned int *pointer_n_triangles,
                        double* Sphere_Cloud_dev,unsigned int *OneDonutMesh_dev,unsigned int *T_temp_dev,unsigned int *count_array_dev,unsigned int *index_offset_array_dev){
    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      CUDA        ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    cudaError_t cudaerr;
    //1st step. GenerateSphere
    cudaGenerateAZBLK<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    cudaGenerateDonut<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    cudaGenerateSphere<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    //2nd step. Overlap removing
    SupressOverlapCUDA<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    cudaerr=cudaMemcpy(Sphere_Cloud, Sphere_Cloud_dev, sizeof(double) *n_total_points * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to SphereCloud. CudaMalloc value=%i\n\r",cudaerr);
    //6th step. last Fill
    TwoandTri_Donut_Fill(Sphere_Cloud,&T_Sphere[OneDonutFill_triangles*3],&T_Sphere[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T_Sphere[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);
    //3rd step. First part of the ODF
    ODF_part1<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev,OneDonutMesh_dev,T_temp_dev,count_array_dev);
    //4th step. eScan GPU
    eScanGPU<<<1, 512 >>> (index_offset_array_dev,count_array_dev);//try with CPU. Try with masking. Chapter: toolkit pin-memory. Optimized method. Chapter: asynch transfer
    /*
    unsigned int count_array[1024],index_offset_array[1024];
    cudaerr=cudaMemcpy(count_array, count_array_dev, sizeof(unsigned int) *1024, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to count_array. CudaMalloc value=%i\n\r",cudaerr);
    index_offset_array[0]=0;
    for (int z=1;z<1024;z++){
        index_offset_array[z]=index_offset_array[z-1]+count_array[z-1];
    }
    cudaerr=cudaMemcpy(index_offset_array_dev, index_offset_array, sizeof(unsigned int) *1024, cudaMemcpyHostToDevice);
    if (cudaerr != 0)	printf("ERROR copying to index_offset_array_dev. CudaMalloc value=%i\n\r",cudaerr);
    */
     
    
    //5th step. last part ODF
    ODF_part2<<<numBlocks, threadsPerBlock >>> (OneDonutMesh_dev,T_temp_dev,count_array_dev,index_offset_array_dev);
    cudaerr=cudaMemcpy(T_Sphere,OneDonutMesh_dev, sizeof(unsigned int) *OneDonutFill_triangles * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to T_Sphere. CudaMalloc value=%i\n\r",cudaerr);
	/////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      CPU        ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    cudaDeviceSynchronize();
    unsigned int temp_vex,n_triangles=0;
    double xp,yp,zp;
    for (unsigned int i = 0; i < n_total_triangles; i++){
        //Analizamos el punto del vertice v0
        temp_vex=T_Sphere[i*3];
        xp=Point_Cloud[temp_vex*3+0];
        yp=Point_Cloud[temp_vex*3+1];
        zp=Point_Cloud[temp_vex*3+2];
        if ((xp!=0)&&(yp!=0)&&(zp!=0)){
            //analizamos el punto del vertice v1
            temp_vex=T_Sphere[i*3+1];
            xp=Point_Cloud[temp_vex*3+0];
            yp=Point_Cloud[temp_vex*3+1];
            zp=Point_Cloud[temp_vex*3+2];
            if ((xp!=0)&&(yp!=0)&&(zp!=0)){
                //analizamos el punto del vertice v2
                temp_vex=T_Sphere[i*3+2];
                xp=Point_Cloud[temp_vex*3+0];
                yp=Point_Cloud[temp_vex*3+1];
                zp=Point_Cloud[temp_vex*3+2];
                if ((xp!=0)&&(yp!=0)&&(zp!=0)){
                    //Si todo lo anterior se cumple, guardamos el triángulo
                    T[n_triangles*3+2]=temp_vex;
                    T[n_triangles*3+1]=T_Sphere[i*3+1];
                    T[n_triangles*3]=T_Sphere[i*3];
                    n_triangles++;
                }
            }
        }
    }
    pointer_n_triangles[0]=n_triangles; 
}


/******************************************************************/
/*************************       MAIN     *************************/
/******************************************************************/
#define TESTING 1
bool malloc_already=false;
int main()
{
#if TESTING == 0
    #define iter 1000.0
    /*Allocate memory*/
    double* Point_Cloud,*Sphere_Cloud_cpu;
    Point_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    Sphere_Cloud_cpu = (double*)malloc(n_total_points * 3 *sizeof(double));
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
    unsigned int *T_cpu,*T_Sphere_cpu,n_triangles_real_data_cpu;
    T_cpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    T_Sphere_cpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    clock_t startCPU;
	clock_t finishCPU;
    printf ("/********************solo CPU*********************/:\n");
	startCPU = clock();
    for (int i=0;i<iter;i++){
        Generate_surface(Point_Cloud,T_cpu,Sphere_Cloud_cpu,T_Sphere_cpu,&n_triangles_real_data_cpu);
    }
	finishCPU = clock();
	printf("numero de triangulos: %d\n",n_triangles_real_data_cpu);
    printf("CPU: %fms\n", ((double)(finishCPU - startCPU))*1000/(double)CLOCKS_PER_SEC/iter);
	
    /*****************************************************/
    /********************       GPU     ******************/
    /*****************************************************/
    unsigned int *T_gpu,*T_Sphere_gpu,n_triangles_real_data_gpu;
    T_gpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    T_Sphere_gpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    double *Sphere_Cloud_gpu;
    Sphere_Cloud_gpu = (double*)malloc(n_total_points * 3 *sizeof(double));
    

    double* Sphere_Cloud_dev;
    unsigned int *OneDonutMesh_dev,*T_temp_dev,*count_array_dev,*index_offset_array_dev;
    cudaMalloc((void**)(&Sphere_Cloud_dev), sizeof(double) * n_total_points * 3);
    cudaMalloc((void**)(&OneDonutMesh_dev), sizeof(unsigned int) * OneDonutFill_triangles * 3);
    cudaMalloc((void**)(&T_temp_dev), sizeof(unsigned int) * n_triangles_perDonut* (n_donuts-1) * 3);
    cudaMalloc((void**)(&count_array_dev), sizeof(unsigned int) * 1024);
    cudaMalloc((void**)(&index_offset_array_dev), sizeof(unsigned int) * 1024);

    clock_t startGPU;
	clock_t finishGPU;
    printf ("/********************CPU y GPU*********************/:\n");
	startGPU = clock();
    for (int i=0;i<iter;i++){
        Generate_surfaceGPU(Point_Cloud,T_gpu,Sphere_Cloud_gpu,T_Sphere_gpu,&n_triangles_real_data_gpu,
                            Sphere_Cloud_dev,OneDonutMesh_dev,T_temp_dev,count_array_dev,index_offset_array_dev);
    }
	finishGPU = clock();
	printf("numero de triangulos: %d\n",n_triangles_real_data_gpu);
    printf("GPU: %fms\n", ((double)(finishGPU - startGPU))*1000/(double)CLOCKS_PER_SEC/iter);
    cudaFree(Sphere_Cloud_dev);
    cudaFree(OneDonutMesh_dev);
    cudaFree(T_temp_dev);
    cudaFree(count_array_dev);
    cudaFree(index_offset_array_dev); 
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
    free (Point_Cloud);
    free (Sphere_Cloud_cpu);
    free (T_cpu);
    free (T_Sphere_cpu);
    return 0;
#elif TESTING == 1
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
    printf ("/********************  CPU  *********************/:\n");
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
    for (int z=0;z<iter;z++){
        //One_Donut_Fill(Point_Cloud,T);
        unsigned int v0,v1,v2;
        //Definimos los vértices
        //Realizamos la malla triangular para la Donut referencial
        for (unsigned int j = 0; j < n_AZBLK; j++){
            for (unsigned int k = 0; k < n_beams-1; k++)
            {   
                v0=j*n_beams+k;
                v2=v0+1;
                v1=(v0+n_beams+1)&mask;
                T[j*(n_beams-1)*3*2+k*6]=v0;
                T[j*(n_beams-1)*3*2+k*6+1]=v1;
                T[j*(n_beams-1)*3*2+k*6+2]=v2;
                v2=v1;
                v1=v2-1;
                T[j*(n_beams-1)*3*2+k*6+3]=v0;
                T[j*(n_beams-1)*3*2+k*6+4]=v1;
                T[j*(n_beams-1)*3*2+k*6+5]=v2;
            }
        }
    }
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("ODF_pt1 time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    
    start=clock();
    for (int z=0;z<iter;z++){
        //En base a la malla referencial hallamos las demás superficies
        double xp,yp,zp;
        unsigned int temp_vex,count=0;
        for (unsigned int i = 1; i < n_donuts; i++){
            //Analizamos cada vertice del tríangulo
            for (unsigned int j = 0; j < n_triangles_perDonut; j++){
                //Analizamos el punto del vertice v0
                temp_vex=(T[j*3]+i*n_points_perDonut);
                xp=Point_Cloud[temp_vex*3+0];
                yp=Point_Cloud[temp_vex*3+1];
                zp=Point_Cloud[temp_vex*3+2];
                if ((xp!=0)||(yp!=0)||(zp!=0)){
                    //analizamos el punto del vertice v1
                    temp_vex=(T[j*3+1]+i*n_points_perDonut);
                    xp=Point_Cloud[temp_vex*3+0];
                    yp=Point_Cloud[temp_vex*3+1];
                    zp=Point_Cloud[temp_vex*3+2];
                    if ((xp!=0)||(yp!=0)||(zp!=0)){
                        //analizamos el punto del vertice v2
                        temp_vex=(T[j*3+2]+i*n_points_perDonut);
                        xp=Point_Cloud[temp_vex*3+0];
                        yp=Point_Cloud[temp_vex*3+1];
                        zp=Point_Cloud[temp_vex*3+2];
                        if ((xp!=0)||(yp!=0)||(zp!=0)){
                            //Si todo lo anterior se cumple, guardamos el triángulo
                            T[n_triangles_perDonut*3+count*3+2]=temp_vex;
                            T[n_triangles_perDonut*3+count*3+1]=(T[j*3+1]+i*n_points_perDonut);
                            T[n_triangles_perDonut*3+count*3]=(T[j*3]+i*n_points_perDonut);
                            count++;
                        }
                    }
                }
            }
        }
    }
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("ODF_pt2 time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
    printf("total time to compare: %fms\n",timeCPU);
    
    start=clock();
    for (int z=0;z<iter;z++)
    TwoandTri_Donut_Fill(Point_Cloud,&T[OneDonutFill_triangles*3],&T[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);    
    stop=clock();
    timeCPU+=((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter;
    printf("LastFill time: %fms\n", ((double)(stop - start))*1000.0/(double)CLOCKS_PER_SEC/iter);
        
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

    printf ("/*********************  GPU  ********************/\n");
    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++){
        //GenerateSphereCUDA<<<4, 256 >>> (Point_Cloud_dev);
        cudaGenerateAZBLK<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
        cudaGenerateDonut<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
        cudaGenerateSphere<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
    }
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    printf("GS time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;
    
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

    unsigned int *OneMesh_dev;
    unsigned int *flag_array_dev,*index_offset_array_dev;
    cudaMalloc((void**)(&OneMesh_dev), sizeof(unsigned int) * n_triangles_perDonut * 3);
    cudaMalloc((void**)(&flag_array_dev), sizeof(unsigned int) * n_triangles_perDonut*(n_donuts-1));
    cudaMalloc((void**)(&index_offset_array_dev), sizeof(unsigned int) * 1024);
    unsigned int *T_gpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));

    unsigned int *flag_array=(unsigned int*)malloc(n_triangles_perDonut*(n_donuts-1) *sizeof(unsigned int));
    //3rd step. First part of the ODF
    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++){
        int temp_blok=64;
        cudaODF_part1<<<n_AZBLK*(n_beams-1)/temp_blok, temp_blok >>> (Point_Cloud_dev,OneMesh_dev);
        temp_blok=128;
        cudaODF_part2<<<n_triangles_perDonut*(n_donuts-1)/temp_blok, temp_blok >>> (Point_Cloud_dev,OneMesh_dev,flag_array_dev);
        cudaDeviceSynchronize();
        cudaerr=cudaMemcpy(T_gpu,OneMesh_dev, sizeof(unsigned int) *n_triangles_perDonut *3, cudaMemcpyDeviceToHost);
        if (cudaerr != 0)	printf("ERROR copying to OneMesh, index %d. CudaMalloc value=%i\n\r",z,cudaerr);
        cudaerr=cudaMemcpy(flag_array,flag_array_dev, sizeof(unsigned int) *n_triangles_perDonut*(n_donuts-1) , cudaMemcpyDeviceToHost);
        if (cudaerr != 0)	printf("ERROR copying to flag_array, index %d. CudaMalloc value=%i\n\r",z,cudaerr);
        
        cpuODF_part3(Point_Cloud_gpu,T_gpu,flag_array);
        
        ////4th step. eScan GPU
        //eScanGPU<<<1, 1024/2 >>> (index_offset_array_dev,count_array_dev);
        ////5th step. last part ODF
        //ODF_part2<<<numBlocks, threadsPerBlock >>> (OneMesh_dev,OneMesh_temp_dev,count_array_dev,index_offset_array_dev);
        //cudaerr=cudaMemcpy(T_gpu,OneMesh_dev, sizeof(unsigned int) *OneDonutFill_triangles * 3, cudaMemcpyDeviceToHost);
    }
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    printf("ODF time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;

    TwoandTri_Donut_Fill(Point_Cloud_gpu,&T_gpu[OneDonutFill_triangles*3],&T_gpu[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T_gpu[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);    
    printf("total time to compare: %fms\n",totalGPU);

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
    cudaFree(flag_array_dev);
    cudaFree(index_offset_array_dev);
    free(Point_Cloud_gpu);
    free(T_gpu); 
    return 0;
#else
    #define iter 100
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
    for (int z = 0; z < iter; z++){
        cudaGenerateAZBLK<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
        cudaDeviceSynchronize();
        cudaGenerateDonut<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
        cudaDeviceSynchronize();
        cudaGenerateSphere<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
        cudaDeviceSynchronize();
    }
    cudaerr=cudaMemcpy(Point_Cloud_gpu, Point_Cloud_dev, sizeof(double) *n_total_points * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to Point_Cloud_gpu. CudaMalloc value=%i\n\r",cudaerr);
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    printf("GS time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;

    FILE* archivo;
    archivo = fopen("../testfiles/CUDASphere_cloud.csv", "w+");
    fprintf(archivo, "X, Y, Z\n");
    for (unsigned int i=0; i < n_total_points; i++) {
        fprintf(archivo,"%.4f, %.4f, %.4f\n", Point_Cloud_gpu[i*3+0], Point_Cloud_gpu[i * 3 + 1], Point_Cloud_gpu[i * 3 + 2]);
    }
    fclose(archivo);
    return;
#endif
}
