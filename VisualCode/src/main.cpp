#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "ouster_reconstruction.h"

int main(){
     /*Allocate memory*/
    double *Point_Cloud,*Sphere_Cloud;
    Point_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    Sphere_Cloud=(double*)malloc(n_total_points * 3 *sizeof(double));
    FILE* archivo;
    
 #define FOR_TESTING 0
 #if FOR_TESTING==1
    unsigned int *T_ODF, *T_TwoDF,*T_TriDF,*T_MidDF;
    T_ODF=(unsigned int*)malloc(OneDonutFill_triangles*3*sizeof(unsigned int));
    T_TwoDF=(unsigned int*)malloc(TwoDonutFill_triangles*3*sizeof(unsigned int));
    T_TriDF=(unsigned int*)malloc(TriDonutFill_triangles*3*sizeof(unsigned int));
    T_MidDF=(unsigned int*)malloc(MidDonutFill_triangles*3*sizeof(unsigned int));
    Generate_sphere(Point_Cloud);
    Supress_redundant_data(Point_Cloud);
    One_Donut_Fill(Point_Cloud,T_ODF);
    TwoandTri_Donut_Fill(Point_Cloud,T_TwoDF,T_TriDF,T_MidDF);
    /*Escribimos la data obtenida en un archivo csv*/
    /*Creamos el csv de la esfera sin traslape*/
    archivo = fopen("files/Sphere_cloud.csv", "w+");
    fprintf(archivo, "X, Y, Z\n");
    for (unsigned int i=0; i < n_total_points; i++) {
        fprintf(archivo,"%.4f, %.4f, %.4f\n", Point_Cloud[i*3+0], Point_Cloud[i * 3 + 1], Point_Cloud[i * 3 + 2]);
    }
    fclose(archivo);
    /*Creamos el csv del mesh del OneDonutfill*/
    archivo = fopen("files/One_donut_fill.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < OneDonutFill_triangles; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_ODF[i*3+0], T_ODF[i * 3 + 1], T_ODF[i * 3 + 2]);
    }
    fclose(archivo);
    /*Creamos el csv del mesh del TwoDonutfill*/
    archivo = fopen("files/Two_donut_fill.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < TwoDonutFill_triangles; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_TwoDF[i*3+0], T_TwoDF[i * 3 + 1], T_TwoDF[i * 3 + 2]);
    }
    fclose(archivo);
    /*Creamos el csv del mesh del TriDonutfill*/
    archivo = fopen("files/Tri_donut_fill.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < TriDonutFill_triangles; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_TriDF[i*3+0], T_TriDF[i * 3 + 1], T_TriDF[i * 3 + 2]);
    }
    fclose(archivo);
    /*Creamos el csv del mesh del MidDonutfill*/
    archivo = fopen("files/Mid_donut_fill.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < MidDonutFill_triangles; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_MidDF[i*3+0], T_MidDF[i * 3 + 1], T_MidDF[i * 3 + 2]);
    }
    fclose(archivo);
    free(T_ODF);
    free(T_TwoDF);
    free(T_TriDF);
    free(T_MidDF);
#else
    /*****************************************************/
    /********************       CPU     ******************/
    /*****************************************************/
    //Leemos del csv los datos reales
    archivo = fopen("files/MinaData.csv", "r");
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
    //Obtenemos la reconstruccion
    unsigned int *T,*T_Sphere,n_triangles_real_data;
    T=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    T_Sphere=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    
    //evaluamos tiempo
    clock_t startCPU;
	clock_t finishCPU;
	printf("\nTiempo de ejecucion en CPU:\n");
	startCPU = clock();
    for (int i=0;i<500;i++){
        Generate_surface(Point_Cloud,T,Sphere_Cloud,T_Sphere,&n_triangles_real_data);
    }
	finishCPU = clock();
	printf("numero de triangulos: %d\n",n_triangles_real_data);
    printf("CPU: %fms\n", ((double)(finishCPU - startCPU))/(double)CLOCKS_PER_SEC);
	//Creamos el archivo csv para guardar el resultado
    archivo = fopen("files/MinaTriangleMesh.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < n_triangles_real_data; i++) {
        fprintf(archivo,"%d, %d, %d\n", T[i*3+0], T[i * 3 + 1], T[i * 3 + 2]);
    }
    fclose(archivo);
    //------------------------------------------
	//----------Generate the DXF file-----------
	//------------------------------------------
	//Open the DXF file
    archivo = fopen("files/MinaSurface.dxf", "w");
    //assert(archivo);
    //header
    fprintf(archivo, "0\nSECTION\n2\nENTITIES\n0\n");
	float x0,y0,z0,x1,y1,z1,x2,y2,z2;
	for (int i = 0; i < n_triangles_real_data; i++)
    {
        // get the coordinates of each point from the triangle
        x0 = Point_Cloud[T[i* 3+0]*3 + 0];
        y0 = Point_Cloud[T[i* 3+0]*3 + 1];
        z0 = Point_Cloud[T[i* 3+0]*3 + 2];
        
        x1 = Point_Cloud[T[i* 3+1]*3 + 0];
        y1 = Point_Cloud[T[i* 3+1]*3 + 1];
        z1 = Point_Cloud[T[i* 3+1]*3 + 2];
        
        x2 = Point_Cloud[T[i* 3+2]*3 + 0];
        y2 = Point_Cloud[T[i* 3+2]*3 + 1];
        z2 = Point_Cloud[T[i* 3+2]*3 + 2];
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
    //Finalmente, liberamos el resto de memoria
    free(Point_Cloud);
    free(Sphere_Cloud);
    free(T);   
    free(T_Sphere);   
#endif 
    return 0;
}
