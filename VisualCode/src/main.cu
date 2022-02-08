#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "math.h"
#include "time.h"
#ifndef MPI//definimos el valor de PI
#define MPI 3.14159265358979323846
#define MPI_2 1.57079632679489661923
#endif
#define D180_MPI 0.017453293 //: Degrees180/pi
/*----------------------------------------------------------------------------*/
/**
 * Definimos los parámetros del LIDAR
 * Se puede usar typedef para juntar los parámetros del LIDAR, o mejor no.
 */
/*Cantidad de rayos por azimuth*/
#define n_beams 16
/*Cantidad de azimuths por Donut*/
#define n_AZBLK 1024
/*Los ángulos por defecto de fabrica del LIDAR*/
double beam_altitude_angles[n_beams]= {15.379*D180_MPI,13.236*D180_MPI,11.128*D180_MPI,9.03*D180_MPI,6.941*D180_MPI,4.878*D180_MPI,2.788*D180_MPI,0.705*D180_MPI,-1.454*D180_MPI,-3.448*D180_MPI,-5.518*D180_MPI,-7.601*D180_MPI,-9.697*D180_MPI,-11.789*D180_MPI,-13.914*D180_MPI,-16.062*D180_MPI};
double beam_azimuth_angles[n_beams] = { -1.24*D180_MPI, -1.2145*D180_MPI, -1.1889*D180_MPI, -1.1634*D180_MPI, -1.1379*D180_MPI, -1.1123*D180_MPI, -1.0868*D180_MPI, -1.0613*D180_MPI, -1.0357*D180_MPI, -1.0102*D180_MPI, -0.98467*D180_MPI, -0.95913*D180_MPI, -0.9336*D180_MPI, -0.90807*D180_MPI, -0.88253*D180_MPI, -0.857*D180_MPI };
/*Cantidad puntos de la nube*/
#define n_points_perDonut (unsigned int)(n_AZBLK*n_beams)
/*Angulo entre azimuths*/
#define angle_between_azimuths (-2*MPI/n_AZBLK)
/*Ángulo de rotacion del motor*/
#define rot_angle (-33.53706667 * MPI / 180)
/*Creamos la matrix de rotación del eje del motor (eje Z)*/
double rot_motor_matrix[9]={ cos(rot_angle),-sin(rot_angle),0,sin(rot_angle),cos(rot_angle) ,0,0,0,1 };
/*Cantidad de Donuts en función del ángulo de rotación*/
#define n_donuts 6//(unsigned int) ceil(-MPI/rot_angle)
#define n_total_points (unsigned int)(n_donuts*n_points_perDonut)
/*Cantidad de triángulos*/
#define OneDonutFill_triangles          108226
#define TwoDonutFill_triangles          6883
#define TriDonutFill_triangles          360
#define MidDonutFill_triangles          2793
#define n_total_triangles (OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles+MidDonutFill_triangles)

//CUDA libraries
#define threadsPerBlock 1024
#define numBlocks (n_AZBLK/1024)
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#ifndef __CUDACC__  
	#define __CUDACC__
	#include <device_functions.h>
#endif

/*----------------------------------------------------------------------------*/
/**
 * Funciones de conversión
 */
void rad2deg(double *value) {
    *value = (*value) * 180 / MPI;
}
void deg2rad(double *value) {
    value[0] = value[0] * MPI / 180;
}
void mult_matrix(double* A,unsigned int m,unsigned int n,double*B,unsigned int l,double*C) {
    for (unsigned int i = 0; i < l;i++) {
        for (unsigned int j = 0; j < m; j++) {
            C[j*l+i]=0;
            for (unsigned int k = 0; k < n; k++) {
                C[j * l +i] += A[j*n+k]*B[k*l + i];
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
/** 
 * Funciones de matrices de rotación
 */
void rot_x_axis(double* XYZ_points,double angle ){
    double temp[3], rot_matrix[9] = { 1,0,0,0,cos(angle),-sin(angle),0,sin(angle),cos(angle) };
    memcpy(temp, XYZ_points, 3*sizeof(double));
    mult_matrix(rot_matrix, 3, 3, temp, 1, XYZ_points);
}
void rot_y_axis(double* XYZ_points, double angle) {
    double temp[3], rot_matrix[9] = { cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle) };
    memcpy(temp, XYZ_points, 3 * sizeof(double));
    mult_matrix(rot_matrix, 3, 3, temp, 1, XYZ_points);
}
void rot_z_axis(double* XYZ_points, double angle) {
    double temp[3], rot_matrix[9] = { cos(angle),-sin(angle),0,sin(angle),cos(angle) ,0,0,0,1};
    memcpy(temp, XYZ_points, 3 * sizeof(double));
    mult_matrix(rot_matrix, 3, 3, temp, 1, XYZ_points);
}
/*----------------------------------------------------------------------------*/
/**
 * \brief Generate synthetic sphere. La generación se divide en 3 estapas. La 
 * 1ra etapa consta en la generación de un azimut referencial; para ello, se 
 * define que el azimuth referencial será aquel que se genera en la vertical
 * inferior, luego, ubicamos cada punto en el plano XZ y luego le realizamos la
 * rotación en el eje X debido al desfase de rayos en cada azimut. En la 2da 
 * etapa se realiza el barrido del azimut referencial para obtener la Donut 
 * refernecial, y en la 3era etapa se realiza la rotación a la Donut referencial
 * para obtener la esfera modelada.
 * 
 * \param Point_Cloud es el puntero que tendrá los puntos de la nube de puntos
 * de la esfera sintética modelada
 * 
 * \return None
 */
#define Radius_sphere 1.0
void Generate_sphere(double* Point_Cloud) {
    /**
    * Nuestro sistema de referencia será : Eje Z será el eje de giro del motor.
    *                                      Eje X será el eje de la Donut referencial
    */
    double R= Radius_sphere;
    /*Generamos el azimuth refencial y los azimuts de cada sector*/
    for (int i = 0; i < n_beams; i++) {
        //ubicamos el punto del azimut referencial en el plano XZ
        Point_Cloud[3 * i + 0] = R* cos(beam_altitude_angles[i] - MPI_2); //x
        Point_Cloud[3 * i + 1] = 0;                                         //y
        Point_Cloud[3 * i + 2] = R * sin(beam_altitude_angles[i] - MPI_2);//z
        //Realizamos la rotacion del punto con respecto al eje x debido al desfase
        rot_x_axis(&Point_Cloud[3*i],beam_azimuth_angles[i]);
        /*Creamos los azimuts que inician en cada sector*/
        /*mirror points from quarter Donut*/
        Point_Cloud[(i + n_AZBLK / 4 * n_beams)*3+0] = Point_Cloud[3 * i + 0];
        Point_Cloud[(i + n_AZBLK / 4 * n_beams)*3+1] = Point_Cloud[3 * i + 2];
        Point_Cloud[(i + n_AZBLK / 4 * n_beams)*3+2] = -Point_Cloud[3 * i + 1];
        /*mirror points from midle Donut*/
        Point_Cloud[(i + n_AZBLK / 2 * n_beams)*3+0] = Point_Cloud[3 * i + 0];
        Point_Cloud[(i + n_AZBLK / 2 * n_beams)*3+1] = -Point_Cloud[3 * i + 1];
        Point_Cloud[(i + n_AZBLK / 2 * n_beams)*3+2] = -Point_Cloud[3 * i + 2];
        /*mirror points from 3 quater Donut*/
        Point_Cloud[(i + n_AZBLK * 3 / 4 * n_beams)*3 +0] = Point_Cloud[3 * i + 0];
        Point_Cloud[(i + n_AZBLK * 3 / 4 * n_beams)*3 +1] = -Point_Cloud[3 * i + 2];
        Point_Cloud[(i + n_AZBLK * 3 / 4 * n_beams)*3 +2] = Point_Cloud[3 * i + 1];

    }
    /*Definimos la matrix de rotación para los azimuth*/
    double rot_matrix[9] = { 1,0,0,0,cos(angle_between_azimuths),-sin(angle_between_azimuths),0,sin(angle_between_azimuths),cos(angle_between_azimuths) };
    /*Procedemos a realizar el barrido de cada sector para obtener la donnut referencial*/
    double XYZ[3],temp[9];
    for (int i = 1; i < n_AZBLK / 4; i++) {
        for (int j = 0; j < n_beams;j++) {
            //Calculate previous point
            XYZ[0] = Point_Cloud[((i - 1) * n_beams + j)*3 + 0];
            XYZ[1] = Point_Cloud[((i - 1) * n_beams + j)*3 + 1];
            XYZ[2] = Point_Cloud[((i - 1) * n_beams + j)*3 + 2];
            //rotate that point
            mult_matrix (rot_matrix,3,3,XYZ,1,temp);
            //Set the new azimuth
            Point_Cloud[(i * n_beams + j)*3 + 0] = temp[0];
            Point_Cloud[(i * n_beams + j)*3 + 1] = temp[1];
            Point_Cloud[(i * n_beams + j)*3 + 2] = temp[2];
            //mirror from quarter Donunt
            Point_Cloud[((i + n_AZBLK / 4) * n_beams + j)*3 + 0] = temp[0];
            Point_Cloud[((i + n_AZBLK / 4) * n_beams + j)*3 + 1] = temp[2];
            Point_Cloud[((i + n_AZBLK / 4) * n_beams + j)*3 + 2] = -temp[1];
            //mirror points from midle Donut
            Point_Cloud[((i + n_AZBLK / 2) * n_beams + j)*3+0] = temp[0];
            Point_Cloud[((i + n_AZBLK / 2) * n_beams + j)*3+1] = -temp[1];
            Point_Cloud[((i + n_AZBLK / 2) * n_beams + j)*3+2] = -temp[2];
            //mirror points from midle Donut
            Point_Cloud[((i + n_AZBLK * 3 / 4) * n_beams + j)*3 + 0] = temp[0];
            Point_Cloud[((i + n_AZBLK * 3 / 4) * n_beams + j)*3 + 1] = -temp[2];
            Point_Cloud[((i + n_AZBLK * 3 / 4) * n_beams + j)*3 + 2] = temp[1];
        }
    }
    //Rotamos la Donut referencial
    for (unsigned int i = 1; i < n_donuts; i++) {
        //multiplicamos a todos los n_point_perdonut anteriores con la matriz de rotación
        for (unsigned int j = 0; j < n_points_perDonut; j++) {
            mult_matrix(rot_motor_matrix, 3, 3, &Point_Cloud[((i - 1) * n_points_perDonut + j) * 3], 1, temp);
            Point_Cloud[(i * n_points_perDonut + j) * 3 + 0] = temp[0];
            Point_Cloud[(i * n_points_perDonut + j) * 3 + 1] = temp[1];
            Point_Cloud[(i * n_points_perDonut + j) * 3 + 2] = temp[2];
;        }
    }
}
/*----------------------------------------------------------------------------*/
/**
 * \brief Ecuación de la recta. Dado los parametros de una recta, se devuelve el
 * valor del punto y que corresponde al punto x
 * 
 * \param m.Es la pendiente de la recta
 * \param x.Es el punto x al cual se evalua para obtener el punto y correspondiente
 * \param xb.Punto xo donde pasa la recta
 * \param yb.Punto yo donde pasa la recta
  * 
 * \return y. resultado de realizar la ecuación
 */
double eq_line(double m,double x,double xb,double yb) {
    double y= m*(x-xb)+yb;
    return y;
}
/*----------------------------------------------------------------------------*/
/**
 * \brief Supress redundant data. En esta función se realiza la supresión de los 
 * puntos. Primero se debe definir los límites de la donut referencial y con 
 * respecto a las siguientes donuts se analiza el límite de la donut anterior. 
 * Definidos los límites, se eliminan aquello puntos que se encuentren en zonas 
 * no permitidas 
 * 
 * \param Point_Cloud es el puntero que tendrá los puntos de la nube de puntos ha 
 * ser reducida 
 * 
 * \return None
 */
void Supress_redundant_data(double* Point_Cloud){
    //Declare tomporary variables
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
    for (unsigned int j = 0; j < n_points_perDonut; j++){
        //Calculamos si la coordenada del punto x
        x=Point_Cloud[i*n_points_perDonut*3+j*3];
        //Analizamos si se encuentra en la zona de supresión
        if (X_Ln<=x&&x<=X_L1){
            Point_Cloud[i*n_points_perDonut*3+j*3]=0;
            Point_Cloud[i*n_points_perDonut*3+j*3+1]=0;
            Point_Cloud[i*n_points_perDonut*3+j*3+2]=0;
            //points_deleted=points_deleted+1;
        }
    }
    //Supress redundant for the rest of the Donuts
    //Creamos variable booleana para saber la zona del plano en donde se encuentra el punto
    bool left_side;
    for (unsigned int i = 2; i < n_donuts; i++){
        for (unsigned int j = 0; j < n_points_perDonut; j++){
            //Hallamos las coordenadas del punto a analizar
            x=Point_Cloud[i*n_points_perDonut*3+j*3];
            y=Point_Cloud[i*n_points_perDonut*3+j*3+1];
            //Evaluamos si se encuentra en la zona referencial
            if (X_Ln<=x){
                if(x<=X_L1){
                    Point_Cloud[i*n_points_perDonut*3+j*3]=0;
                    Point_Cloud[i*n_points_perDonut*3+j*3+1]=0;
                    Point_Cloud[i*n_points_perDonut*3+j*3+2]=0;
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
            y_temp=eq_line(L[(i-2)*5+4],x,L[(i-2)*5+left_side*2],L[(i-2)*5+left_side*2+1]);
            //Le colocamos un signo negativo, o no, para poder realizar un único condicional para ambos casos
            y_temp=y_temp*(1.0-2*left_side);
            y=y*(1.0-2*left_side);
            //Evaluamos la condición de supresión
            if (y>=y_temp){
                //Eliminamos los puntos
                Point_Cloud[i*n_points_perDonut*3+j*3]=0;
                Point_Cloud[i*n_points_perDonut*3+j*3+1]=0;
                Point_Cloud[i*n_points_perDonut*3+j*3+2]=0;
                //points_deleted=points_deleted+1;
            }
        }   
    }
}
/*----------------------------------------------------------------------------*/
/**
 * \brief One-Donut-Fill. Dado una nube de puntos sin traslape y siguiendo el 
 * patrón de medición definido por el sistema Ouster-motor, se obtiene como 
 * salida la malla triangular de las superficies que pertenecen a una única Donut.
 * Primero se realiza la triangulación de la Donut referencial, luego en base a 
 * esta triangulación se obtiene la de las demás Donuts, pero verificando si los
 * vértices son distintos de cero
 * 
 * \param Point_Cloud es el puntero que tendrá los puntos de la nube de puntos  
 * 
 * \param T es el puntero donde se almacenará los vértices de los triángulos
 * 
 * \return None
 */
#define mask (n_points_perDonut-1)
void One_Donut_Fill(double* Point_Cloud,unsigned int* T){
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
    //En base a la malla referencial hallamos las demás superficies
    double xp,yp,zp;
    unsigned int temp_vex,count=0,n_triangles_perDonut=n_AZBLK*2*(n_beams-1);
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
    //printf("numero de triangulos: %d\n",count+n_triangles_perDonut);    
}
/*----------------------------------------------------------------------------*/
/**
 * \brief get_tripivot. Dado el punto del vértice de la Donut, se hallan los otros
 * vértices del triángulo pívot. 
 * 
 * \param v0_pointer.  puntero del vértice v0
 * 
 * \param v1_pointer. puntero del vértice v1
 * 
 * \param point. Puntero que contiene las coordenadas del vertice v2 del tripivot
 * 
 * \param sector. Indica en que sector estamos de la Donut
 * 
 * \param i. Indica para que Donut se está hallando el tripivot
 * 
 * \param k_beam. Nos dice si nos enfocamos con el beam=0 o el beam=15
 * 
 * \return None
 */
void get_tripivot(unsigned int *vmin_pointer,unsigned int *vmax_pointer,double *point,unsigned int sector,unsigned int i,unsigned int k_beam){
    //Definimos las variables
    unsigned int offset,k_azimuth,vmin,vmax;
    double y_data,z_data,theta,rot_theta,alfa;
    double rot_point[3];
    //Se debe tener en cuenta los rango de acos y asin
    //  acos: [0~pi]
    //  atan: [-pi/2~pi/2]
    rot_point[0]=point[0];
    rot_point[1]=point[1];
    rot_point[2]=point[2];
    //Con el valor de x_data, podemos hallar el ángulo del punto con respecto a 
    //la Donut referencial, teniendo como eje de giro el eje z
    theta=-acos(point[0]/Radius_sphere);
    //Asimismo, debemos hallar el ángulo con respecto al al Donut previa, para ello
    //realizamos un antigiro para que la Donut previa "aparente" ser la referencial y
    //usar las funciones ya establecidas
    rot_z_axis(rot_point, -rot_angle*(i-1));
    rot_theta=-acos(rot_point[0]/Radius_sphere);
    if (theta*(1.0-2.0*((sector>>1)&0x1)) >= rot_theta*(1.0-2*((sector>>1)&0x1)))
    {
        //Los vértices pertenecen a la Donut referencias
        y_data=point[1];
        z_data=point[2];
        offset=0;
    }else{
        //los vértices pertenecerán a la previa Donut
        y_data=rot_point[1];
        z_data=rot_point[2];
        theta=rot_theta;
        offset=n_AZBLK*n_beams*(i-1);
    }
    //Hallamos el ángulo alfa, que es el angulo del azimuth
    alfa=-MPI_2 +atan(z_data/y_data);
    if (y_data>0) alfa=-MPI+alfa;
    //Añadimos -2pi al alfa para no tener probelmas con el bitwise and,
    //esto al final no perjudica ya que se hace el masking de bits
    //solo que tener en cuenta que la división entre -2pi/ang_bet_azit
    //da un total de 1024.
    alfa=alfa-2*MPI;
    //Con lo anterior, nos hemos asegurado que el alfa sea siempre negativo
    //Calculamos el azimuth que corresponde al alfa
    k_azimuth=(unsigned int)((alfa-beam_azimuth_angles[k_beam])/angle_between_azimuths);
    vmin=(k_azimuth)*n_beams+k_beam;
    //Enmascaramos
    vmin=vmin&mask;
    vmax=(vmin+n_beams)&mask;
    //Realizamos offset y sentido
    vmin+=offset;
    vmax+=offset;
    *vmin_pointer=vmin;
    *vmax_pointer=vmax;
                
}
/*----------------------------------------------------------------------------*/
/**
 * \brief side2sideFill. Dado una nube de puntos sin traslape y siguiendo 
 * el patrón de medición definido por el sistema Ouster-motor, se obtiene como 
 * 
 * \param parameter-name description
 * 
 * \return None
 */
void side2sideFill( unsigned int vL0_init,unsigned int vL1_init,
                    unsigned int vL0_fin,unsigned int vL1_fin,
                    int pasoL0,int pasoL1,
                    unsigned int* T,unsigned int* n_triangles_pointer)
{
    unsigned int vmax,vmin,v0,v1,v2,v1_fin,v_temp;//,v2_fin;
    int freepointsL0,freepointsL1,pasov1,pasov2,arista_mismo_vex;
    //Realizamos algunos ajustes para poder obterner el valor magnitud de
    //los puntos libres tanto para la izquierda y derecha
    if (pasoL0<0){
        vmax=vL0_init&mask;
        vmin=vL0_fin&mask;
    }else{
        vmax=vL0_fin&mask;
        vmin=vL0_init&mask;
    }
    freepointsL0=((vmax+(mask-vmin)+1)&mask)/abs(pasoL0);
    if (pasoL1<0){
        vmax=vL1_init&mask;
        vmin=vL1_fin&mask;
    }else{
        vmax=vL1_fin&mask;
        vmin=vL1_init&mask;
    }
    freepointsL1=((vmax+(mask-vmin)+1)&mask)/abs(pasoL1);
    //Con las siguientes formulas podemos hallar con que Donut 
    //estamos trabajando y obtener el offset adecuado
    unsigned int offsetL0=n_points_perDonut*floor(double((vL0_init-1))/n_points_perDonut);
    unsigned int offsetL1=n_points_perDonut*floor(double((vL1_init-1))/n_points_perDonut);
    unsigned int offset;
    //considero que "v2 es de L0" y "v1 es de L1"
    v2=vL0_init;
    v1=vL1_init;
    //v2_fin=vL0_fin;
    v1_fin=vL1_fin;
    pasov2=pasoL0;
    pasov1=pasoL1;
    //la siguiente variable es para el caso que el vertice en comun este en L2
    bool volteamos=false;
    //Analizamos si tenemos distintos puntos
    if (freepointsL1!=freepointsL0){
        //Realizamos un triangulo con mismo vertice
        arista_mismo_vex=freepointsL1-freepointsL0;
        offset=offsetL1;
        if (arista_mismo_vex<0){
            //En caso L1 tenga mas puntos el v_comun estará en L2
            v2=vL1_init;
            v1=vL0_init;
            pasov1=pasoL0;
            offset=offsetL0;
            arista_mismo_vex=arista_mismo_vex*-1;
            volteamos=true;
        }
        for (unsigned int j=0;j<arista_mismo_vex;j++){
            v0=v1;
            //debemos realizar el offset adecuado segun la Donut con la que
            //trabajamos. Primero le sumamos el mask y luego el paso, esto para obtener
            //la concatenación en la misma Donut. Luego le sumamos el offset de la Donut
            //le sumamos la unidad para que luego de hacer el bitmasking obtengamos el mismo
            //numero
            v1=(((v1+mask+1)+pasov1)&mask)+offset;
            T[n_triangles_pointer[0]*3+0]=v0;
            T[n_triangles_pointer[0]*3+1]=v1;
            T[n_triangles_pointer[0]*3+2]=v2;
            n_triangles_pointer[0]++;
        }
        if (volteamos){
            //en caso habiamos volteado, volvemos al caso inicial
            v_temp=v2;
            v2=v1;
            v1=v_temp;
            pasov1=pasoL1;
        }
    }
    //Realizamos el llenado "alineado"
    while (v1!=v1_fin){
        v0=v1;
        //Con esta formula podemos obtener la concatenación de Donuts
        v1=(((v2+mask+1)+pasov2)&mask)+offsetL0;
        T[n_triangles_pointer[0]*3+0]=v0;
        T[n_triangles_pointer[0]*3+1]=v1;
        T[n_triangles_pointer[0]*3+2]=v2;
        n_triangles_pointer[0]++;
        v2=v1;
        //Esta formula permite avanzar por medio del enmascaramiento sin
        //necesidad de usar condicionales (equivalente al operador modulo)
        v1=(((v0+mask+1)+pasov1)&mask)+offsetL1;
        T[n_triangles_pointer[0]*3+0]=v0;
        T[n_triangles_pointer[0]*3+1]=v1;
        T[n_triangles_pointer[0]*3+2]=v2;
        n_triangles_pointer[0]++;
    }
}

/*----------------------------------------------------------------------------*/
/**
 * \brief Two-and-Tri-Donut-Fill. Dado una nube de puntos sin traslape y siguiendo 
 * el patrón de medición definido por el sistema Ouster-motor, se obtiene como 
 * salida la malla triangular de las superficies que pertenecen a dos y tres Donut.
 * Para realizar esto, debemos seguir una serie de pasos:
 * 
 * Primero, se debe hallar los triángulos pivots, pero para ello se deben hallar
 * los triángulos de cada sector. Es por eso que primero establecemos unos valores
 * de índice inicial el cual estarán ubicado en el "centro de la esfera, estos 
 * valores de indice serán aumentados/reducidos para hallar el vértice no nulo que
 * formará parte del triángulo pivot.
 * 
 * Luego,
 * 
 * \param Point_Cloud es el puntero que tendrá los puntos de la nube de puntos  
 * 
 * \param T es el puntero donde se almacenará los vértices de los triángulos
 * 
 * \return None
 */
void TwoandTri_Donut_Fill(double* Point_Cloud,unsigned int* TwoDF,unsigned int* TriDF,unsigned int* MidDF){
    int paso;
    unsigned int n_tripivot,init_index,k_beam;
    unsigned int n_Twotriangles=0,n_Tritriangles=0,n_Midtriangles=0;
    double x_point,y_point,z_point;
    unsigned int Tripivot[16*3];
    //Creamos el arreglo que contendrá a los cuatro triángulos de cada sector que limitan la zona del medio
    unsigned int Tripivot_middle[4*3];
    //Para la última Donut se deben almacenar otros 4 triángulos
    unsigned int Tripivot_middle_particular[4*3];
    //last_tripivots almacena los últimos triángulos pivots de lso 4 sectores
    unsigned int last_Tripivots[4*3]={0,0,0,0,0,0,0,0,0,0,0,0};
    bool vex_found;
    //variables para clasificar la zona de llenado
    char escalera,tipo;
    //Variables de paso
    int pasoL0,pasoL1;
    //Variables para operar con los vértices
    unsigned int v0,v1,v2,v_temp,v_corner,v_corner_fin,v_init,v_fin;
    //Variables para operar con los vertices en el TriDonutFill
    unsigned int v0_mid,v1_mid,v2_mid,T_mid[3],v0_lim,v1_lim,index;//v2_lim,new_v_init;
    int paso_mid,freepoints_mid;
    unsigned int a,b,c;
    //Variables para el middlefill
    unsigned int vrigth_init,vrigth_fin,vleft_init,vleft_fin;
    int pasoRigth,pasoLeft;
    //A partir de la segunda Donut generamos las superficies
    for (int i = 1; i < n_donuts; i++){
        for (unsigned int sector = 0; sector < 4; sector++){
            n_tripivot=0;
            /**************************************************************************************/
            /****************************   Creamos los tripivots   *******************************/
            /**************************************************************************************/
            for (unsigned int j = 0; j < n_beams; j++){
                //Debemos hallar el vertice límite distinto de cero, para ello usamos un 
                //init_index que nos permitirá evaluar cada punto hasta hallar el que es 
                //distinto de cero. Para hacer esto, el init_index debe iniciar en la parte
                //inferior o superior de la esfera, es decir, donde los puntos son igual a 
                //cero. init_index será (n_points-n_beams)+j o (n_points/2-n_beams)+j, luego
                //se le añade el offset correspondiente a cada Donut.
                //Segun el sector, podemos definir si debería iniciar en la parte superior o
                //inferior:
                init_index=(n_AZBLK>>1)*(2-(((sector>>1)&0x1)^(sector&0x1)))-1;
                init_index=init_index*n_beams+j;
                //Definimos si el paso será positivo o negativo, esto dependerá del sector.
                paso=n_beams*(1.0-2*(sector&0x1));
                //Creamos bandera para saber si se halló un vértice para crear el tripivot
                vex_found=false;
                for (unsigned int count = 0; count < (n_AZBLK>>2); count++){
                    //Entramos en un bucle donde evaluamos los puntos de la Donut hasta hallar
                    //el punto que sea distinto de cero
                    x_point=Point_Cloud[i*n_points_perDonut*3+init_index*3+0];  
                    y_point=Point_Cloud[i*n_points_perDonut*3+init_index*3+1];
                    z_point=Point_Cloud[i*n_points_perDonut*3+init_index*3+2];
                    //Evaluamos condición
                    if ((x_point!=0)||(y_point!=0)||(z_point!=0)){
                        //Hallamos vertice distinto de cero
                        vex_found=true;
                        //Agregamos el offset de la Donut correspondiente
                        v2=init_index+i*n_points_perDonut;
                        break;
                    }
                    //Realizamos el paso correspondiente a init_index
                    init_index=(init_index+paso)&mask;
                }
                //Evaluamos el caso que no se consiguió el vértice
                if (!vex_found){
                    //Evaluamos para los siguientes beams
                    continue;
                }
                /*Procedemos a hallar los demás vértices del tripivot*/
                //Dependiendo del sector estamos más cerca del beam=0 o el beam=15
                k_beam=(n_beams-1)*(1.0-((sector>>1)&0x1));
                ///Hallamos los otros vértices
                get_tripivot(&v0,&v1,&Point_Cloud[v2*3],sector,i,k_beam);
                /*Debemos establecer el orden de los vértices, es decir horario o antihorario*/
                //Esto depende del sector
                if ((sector>>1)&0x1){
                    //definimos el sentido horario Ya que para los sectors 3 y 4
                    //El sentido de los vértices es distinto a los de los
                    //primeros sector. Entonces para seguir la jerarquía de los
                    //sentidos, cambiamos aqui
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    //%%% Esto podría ser o o importante                        %
                    //%%% Capaz, se puede definir un sentido para un lado y otro%
                    //%%% para los otros sector ()sector3 y sector4)            %
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    v_temp=v0;
                    v0=v1;
                    v1=v_temp;
                }
                /*Debemos verificar que los vértices no sean cero */
                //Analizamos el punto del vertice v0
                x_point=Point_Cloud[v0*3+0];
                y_point=Point_Cloud[v0*3+1];
                z_point=Point_Cloud[v0*3+2];
                if ((x_point==0)&&(y_point==0)&&(z_point==0)){
                    //Guardamos el valor de v1.//v_temp=v1;
                    get_tripivot(&v0,&v_temp,&Point_Cloud[v1*3],sector,i-1,k_beam);
                    if ((sector>>1)&0x1){
                        //Deseo el nuevo v0
                        v0=v_temp;
                    }
                    //Devolvemos el valor de v1.//v1=v_temp;
                }
                //Analizamos el punto del vertice v1
                x_point=Point_Cloud[v1*3+0];
                y_point=Point_Cloud[v1*3+1];
                z_point=Point_Cloud[v1*3+2];
                if ((x_point==0)&&(y_point==0)&&(z_point==0)){
                    get_tripivot(&v_temp,&v1,&Point_Cloud[v0*3],sector,i-1,k_beam);
                    if ((sector>>1)&0x1){
                        //Deseo el nuevo v1
                        v1=v_temp;
                    }
                }
                /*Debemos verificar concurrencia de vértices en los triángulos pívots con
                los últimos triangulos pivots de la Donut anterior con el mismo sector*/
                if (!(((sector>>1)&0x1)^(sector&0x1))){
                    //sector 1 y 4
                    if (v1==last_Tripivots[sector*3+1]){
                        v1=last_Tripivots[sector*3+2];
                    }
                }else{
                    //sector 2 y 3
                    if (v0==last_Tripivots[sector*3+0]){
                        v0=last_Tripivots[sector*3+2]; 
                    }
                }
                Tripivot[n_tripivot*3+0]=v0;
                Tripivot[n_tripivot*3+1]=v1;
                Tripivot[n_tripivot*3+2]=v2;
                n_tripivot++;
            }
            //Guardamos los triángulos hallados 
            for (unsigned int j= 0; j < n_tripivot; j++){
                TwoDF[n_Twotriangles*3+0]=Tripivot[j*3+0];
                TwoDF[n_Twotriangles*3+1]=Tripivot[j*3+1];
                TwoDF[n_Twotriangles*3+2]=Tripivot[j*3+2];
                n_Twotriangles++;
            }
            //Hallamos los 4 triangulos que limitan la zona del medio de cada sector,
            //para ello, guardamos el triángulo más cercano al medio
            Tripivot_middle[sector*3+0]=Tripivot[((n_tripivot-1)*(sector>>1))*3+0];
            Tripivot_middle[sector*3+1]=Tripivot[((n_tripivot-1)*(sector>>1))*3+1];
            Tripivot_middle[sector*3+2]=Tripivot[((n_tripivot-1)*(sector>>1))*3+2];
            if ((i==(n_donuts-1))&&(n_tripivot>1)){
                //Este caso ocurre en la última Donut, también hallar los 4 triangulos
                Tripivot_middle_particular[sector*3+0]=Tripivot[((n_tripivot-1)*(0x1^(sector>>1)))*3+0];
                Tripivot_middle_particular[sector*3+1]=Tripivot[((n_tripivot-1)*(0x1^(sector>>1)))*3+1];
                Tripivot_middle_particular[sector*3+2]=Tripivot[((n_tripivot-1)*(0x1^(sector>>1)))*3+2];
            }
            /**************************************FIN TRIPIVOT************************************/
            //Ahora con los triángulos pivot hallados, procedemos a realizar el llenado de las zonas.
            /**************************************************************************************/
            /**************************************************************************************/
            /**********                                                                ************/
            /**********                 TRIDONUTFILL && TWODONUTFILL                   ************/
            /**********                                                                ************/
            /**************************************************************************************/
            /**************************************************************************************/
            for (unsigned int j = 0; j < n_tripivot-1; j++){
                //Debemos clasifiicar la zona a llenar, por eso debemos obtener el punto siguiente 
                //del vertice v2
                v_temp=Tripivot[j*3+2];                
                //Particularmente, en los sector 0 y 2 hallamos el punto siguiente pero en los otros
                //sectores 1 y 3, hallamos el 17mo punto consecuente. Al analizar si este punto es
                //nulo o no, podemos definir si es una zona tipo escalera o rampa(ver pagina 19 de 
                //la presentación)
                x_point=Point_Cloud[(v_temp+1+n_beams*((sector)&0x1))*3+0];
                y_point=Point_Cloud[(v_temp+1+n_beams*((sector)&0x1))*3+1];
                z_point=Point_Cloud[(v_temp+1+n_beams*((sector)&0x1))*3+2];
                //Evaluamos la condición y definimos la forma de la zona
                if ((x_point==0)||(y_point==0)||(z_point==0)){
                    //El punto es cero, ha sido suprimido
                    escalera=sector&0x1;
                }else{
                    escalera=!(sector&0x1);
                }
                //Definimos el tipo
                tipo=0x1^((sector&0x1)^escalera);
                //Luego de definir la forma y el tipo. Podemos definir el valor de pasoL1 para llegar al v_fin
                pasoL0=n_beams*(1.0-2*(escalera^tipo));
                //Ahora debemos definir el vertice que estará en la esquina y su lugar puede estar en el beam 
                //actual o en el beam siguiente o en el beam del siguiente/anterior azimut
                v_corner=Tripivot[(j+tipo)*3+2]+(1.0-2*tipo)+pasoL0*(0x1^escalera);
                //definimos los límites
                v_init=Tripivot[(j+tipo)*3+1-escalera];
                v_fin=Tripivot[(j+(0x1^tipo))*3+1-(0x1^escalera)];
                v_corner_fin=Tripivot[(j+(0x1^tipo))*3+2];
                pasoL1=(1.0-2*(((sector>>1)&0x1)^tipo))*pasoL0;
                //Para verificar si es de TWO o TRI llenado, se tiene que analizar la condición de
                //que los vértices pertenezcan a una misma donut. Algo particular que podemos 
                //destacar es que si es TriDonutFill, necesariamente uno de los vértices pertenece
                //a la Donuut Referencial. 
                //Para saber si NO es TriDonutFill basta verificar el MSB del vértice, en este caso
                //si n_point_perDonut=(2>>14), habría que analizar el bit 14 y ver que este sea distinto
                //de cero, lo que es lo mismo decir que el vértice no pertene a la Donut referencial
                //Por lo tanto, para verificar si es TwoDonutFill basta cumplir alguna de estas 
                //condiciones:
                //  -Si al hacer el bitshift, ambos vertices tienen un valor distinto de cero
                //  -Si al hacer el bitshift, ambos vertices tienenun  valor igual a cero.
                //Lo anterior se puede hacer con bitwise and: vex&(~mask), para ambos casos
                if ((((~mask)&v_init)==(v_fin&(~mask)))||(((v_fin&(~mask))!=0)&&((v_init&(~mask))!=0))){
                    /* TwoDonutFill*/
                    //Creamos el primer triangulo
                    v0=Tripivot[(j+tipo)*3+2];
                    v1=v_init;
                    v2=v_corner;
                    //Guardamos el triángulo
                    TwoDF[n_Twotriangles*3+0]=v0;
                    TwoDF[n_Twotriangles*3+1]=v1;
                    TwoDF[n_Twotriangles*3+2]=v2;
                    n_Twotriangles++;
                    //Creamos la superficie formada por los dos triángulos pivot
                    side2sideFill(v_corner,v_init,v_corner_fin,v_fin,pasoL0,pasoL1,TwoDF,&n_Twotriangles);
                }else{
                    /* TriDonutFill*/
                    //Definimos un nuevo triángulo, que será el del medio
                    //primero tenemos que partir de un vertice para hallar al
                    //v2_mid. Este vertice de partida será siempre del T_actual
                    //para los sector 1 y 2 y para los sector3y4 será el T_next
                    //notar que se está hallando el tripivot de la Donut
                    //anterior (i-1). Pero para "angulos particulares" el tripivot
                    //pertenece a la siguiente Donut!!! (con suerte, este no es el caso)
                    v2_mid=Tripivot[(j+((sector>>1)&0x1))*3+0];//---->es indiferente si es el vertice 0 o 1, están en la misma recta(?)
                    paso_mid=pasoL1*(1.0-2.0*(tipo^((sector>>1)&0x1)));
                    //no hace falta hacer mask ya que siempre es con la Donut
                    //anterior y los Trifill ocurren a partir de la donut 2, es
                    //decir i>=3. Además, el tema de concatenación ocurre más en 
                    //la donut referencial
                    x_point=Point_Cloud[(v2_mid+paso_mid)*3+0];
                    y_point=Point_Cloud[(v2_mid+paso_mid)*3+1];
                    z_point=Point_Cloud[(v2_mid+paso_mid)*3+2];
                    freepoints_mid=0;
                    while ((x_point!=0)||(y_point!=0)||(z_point!=0)){
                        v2_mid+=paso_mid;
                        freepoints_mid++;//realizado el conteo de cuantos puntos libres
                        x_point=Point_Cloud[(v2_mid+paso_mid)*3+0];
                        y_point=Point_Cloud[(v2_mid+paso_mid)*3+1];
                        z_point=Point_Cloud[(v2_mid+paso_mid)*3+2];    
                    }
                    //Hallamos el triangulo pivot Repetimos el código de líneas arriba
                    //------------------TRIPIVOT MIDDLE----------------//
                    k_beam=(n_beams-1)*(1.0-((sector>>1)&0x1));
                    //realizamos la operación
                    get_tripivot(&v0_mid,&v1_mid,&Point_Cloud[v2_mid*3],sector,i-1,k_beam);
                    if ((sector>>1)&0x1){
                        //definimos el sentido horario Ya que para los sectors 3 y 4
                        //El sentido de los vértices es distinto a los de los
                        //primeros sector. Entonces para seguir la jerarquía de los
                        //sentidos, cambiamos aqui
                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        //%%% Esto podría ser o o importante                        %
                        //%%% Capaz, se puede definir un sentido para un lado y otro%
                        //%%% para los otros sector ()sector3 y sector4)                  %
                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        v_temp=v0_mid;
                        v0_mid=v1_mid;
                        v1_mid=v_temp;
                    }
                    T_mid[0]=v0_mid;
                    T_mid[1]=v1_mid;
                    T_mid[2]=v2_mid;
                    /*------------Fin TRIPIVOT MIDDLE-------------*/
                    /*----------Primer llenado----------*/
                    //Definimos los v_lim Para hacer el primer llenado
                    a=(sector>>1)&0x1;
                    b=sector&0x1;
                    c=escalera;
                    //v0_lim será el vertice limite para el v_corner
                    v0_lim=Tripivot[(j+(tipo^0x1))*3+2];
                    //v1_lim será el vertice limite para el v_init
                    index= 2 - (a^(0x1^(b^c))) - (0x1^((a^b)|c));
                    v1_lim=T_mid[index];
                    //procedemos a realizar el primer llenado
                    //Creamos el primer triangulo
                    v0=Tripivot[(j+tipo)*3+2];
                    v1=v_init;
                    v2=v_corner;
                    TriDF[n_Tritriangles*3+0]=v0;
                    TriDF[n_Tritriangles*3+1]=v1;
                    TriDF[n_Tritriangles*3+2]=v2;
                    n_Tritriangles++;
                    side2sideFill(v2,v1,v0_lim,v1_lim,pasoL0,pasoL1,TriDF,&n_Tritriangles);
                    /*----------Fin Primer llenado----------*/
                    /*----------Segundo llenado----------*/
                    index=((sector>>1)&0x1)^((sector&0x1)^escalera);
                    index= 2- index- (escalera&&(0x1^((sector&0x1)^((sector>>1)&0x1))));
                    //Creamos el triangulo de transicion
                    v0=v1_lim;
                    v1=T_mid[index];//Definimos el nuevo v_init
                    v2=v_corner_fin;
                    TriDF[n_Tritriangles*3+0]=v0;
                    TriDF[n_Tritriangles*3+1]=v1;
                    TriDF[n_Tritriangles*3+2]=v2;
                    n_Tritriangles++;
                    side2sideFill(v2,v1,v_corner_fin,v_fin,pasoL0,pasoL1,TriDF,&n_Tritriangles);
                    /*----------Fin Segundo llenado----------*/
                }
            }
            /*--------------------FIN FILL sector--------------------*/
            //la siguiente variable guardará los tripivots pasados para verificar que 
            //no exista concurrencia de tripivots con la siguiente Donut
            if (i!=n_donuts){
                //hacemos este condicional ya que para la última Donut no
                //necesitamos realizar esto
                last_Tripivots[sector*3+0]=Tripivot[(0x1^((sector>>1)&0x1))*(n_tripivot-1)*3+0];
                last_Tripivots[sector*3+1]=Tripivot[(0x1^((sector>>1)&0x1))*(n_tripivot-1)*3+1];
                last_Tripivots[sector*3+2]=Tripivot[(0x1^((sector>>1)&0x1))*(n_tripivot-1)*3+2];
            }            
        }
        /*Middle Fill*/
        //sector1 con sector2
        vleft_init=Tripivot_middle[2];
        vrigth_init=Tripivot_middle[1];
        vleft_fin=Tripivot_middle[3+2];
        vrigth_fin=Tripivot_middle[3+0];
        pasoLeft=n_beams;
        pasoRigth=pasoLeft;
        side2sideFill(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,MidDF,&n_Midtriangles);
        //sector3 con sector4
        vleft_init=Tripivot_middle[3*3+2];
        vrigth_init=Tripivot_middle[3*3+1];
        vleft_fin=Tripivot_middle[2*3+2];
        vrigth_fin=Tripivot_middle[2*3+0];
        pasoLeft=-n_beams;
        pasoRigth=pasoLeft;
        side2sideFill(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,MidDF,&n_Midtriangles);
        if ((i==(n_donuts-1))&&(n_tripivot>1)){
            vrigth_init=Tripivot_middle_particular[0];
            vleft_init=Tripivot_middle_particular[2];
            vrigth_fin=Tripivot_middle_particular[3+1];
            vleft_fin=Tripivot_middle_particular[3+2];
            pasoRigth=-n_beams;
            pasoLeft=-pasoRigth;
            side2sideFill(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,MidDF,&n_Midtriangles);
            vrigth_init=Tripivot_middle_particular[3*3+0];
            vleft_init=Tripivot_middle_particular[3*3+2];
            vrigth_fin=Tripivot_middle_particular[2*3+1];
            vleft_fin=Tripivot_middle_particular[2*3+2];
            pasoRigth=n_beams;
            pasoLeft=-pasoRigth;
            side2sideFill(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,MidDF,&n_Midtriangles);
        }
    }
}
/*----------------------------------------------------------------------------*/
/**
 * \brief Generate_surface. Genera la superficie
 * 
 * \param Point_Cloud es el puntero que tendrá los puntos de la nube de puntos 
 * de la data real  
 * 
 * \param T es el puntero donde se almacenará los vértices de los triángulos de 
 * la data real
 * 
 * \param pointer_n_triangles. Es un puntero que contendrá el valor de la cantidad 
 * de triángulos generados. Este valor definirá el límite de líneas para los 
 * archivos de reconstrucción
 * 
 * \return None
 */
void Generate_surface(double* Point_Cloud,unsigned int* T,unsigned int *pointer_n_triangles){
    double* Sphere_Cloud;
    Sphere_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    unsigned int *T_temp;
    T_temp=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    Generate_sphere(Sphere_Cloud);
    Supress_redundant_data(Sphere_Cloud);
    One_Donut_Fill(Sphere_Cloud,T_temp);
    TwoandTri_Donut_Fill(Sphere_Cloud,&T_temp[OneDonutFill_triangles*3],&T_temp[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T_temp[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);    
    /*Procedemos a chequear que los vertices son no nulos*/
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
/*----------------------------------------------------------------------------*/
/**
 * \brief multmatrix. function multiplica matrices para ser usado dentro del GPU 
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
/**
 * \brief GenerateSphereCUDA. genera esfera 
 */
__global__ void GenerateSphereCUDA(double* Point_Cloud){
    double beam_altitude_angles[n_beams]= {15.379*D180_MPI,13.236*D180_MPI,11.128*D180_MPI,9.03*D180_MPI,6.941*D180_MPI,4.878*D180_MPI,2.788*D180_MPI,0.705*D180_MPI,-1.454*D180_MPI,-3.448*D180_MPI,-5.518*D180_MPI,-7.601*D180_MPI,-9.697*D180_MPI,-11.789*D180_MPI,-13.914*D180_MPI,-16.062*D180_MPI};
    double beam_azimuth_angles[n_beams] = { -1.24*D180_MPI, -1.2145*D180_MPI, -1.1889*D180_MPI, -1.1634*D180_MPI, -1.1379*D180_MPI, -1.1123*D180_MPI, -1.0868*D180_MPI, -1.0613*D180_MPI, -1.0357*D180_MPI, -1.0102*D180_MPI, -0.98467*D180_MPI, -0.95913*D180_MPI, -0.9336*D180_MPI, -0.90807*D180_MPI, -0.88253*D180_MPI, -0.857*D180_MPI };
    //Obtenemos el ID del thread
    int thid = threadIdx.x + blockIdx.x * blockDim.x;
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
        double rot_matrix[9] = { 1,0,0,0,cos(angle_between_azimuths*thid),-sin(angle_between_azimuths*thid),0,sin(angle_between_azimuths*thid),cos(angle_between_azimuths*thid) };
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
    double rot_motor_matrix[9]={ cos(rot_angle),-sin(rot_angle),0,sin(rot_angle),cos(rot_angle) ,0,0,0,1 };
    for (unsigned int i = 1; i < n_donuts; i++) {
        //multiplicamos a todos los n_point_perdonut anteriores con la matriz de rotación
        for(int j=0;j<n_beams;j++){
            mult_matrix_dev(rot_motor_matrix, 3, 3, &Point_Cloud[((i - 1) * n_points_perDonut + thid*n_beams+j) * 3], 1, temp);
            Point_Cloud[(i * n_points_perDonut + thid*n_beams+j) * 3 + 0] = temp[0];
            Point_Cloud[(i * n_points_perDonut + thid*n_beams+j) * 3 + 1] = temp[1];
            Point_Cloud[(i * n_points_perDonut + thid*n_beams+j) * 3 + 2] = temp[2];        
        }
    }
    __syncthreads();
}

/**
 * \brief SupressOverlapCUDA. funcion en paralelo 
 */

__device__ double eq_line_dev(double m,double x,double xb,double yb) {
    double y= m*(x-xb)+yb;
    return y;
}
/**
 * \brief SupressOverlapCUDA. funcion en paralelo 
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
    __syncthreads();
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
 * \brief OneDonutFillCUDA. First Mesh 
 */
__global__ void OneDonutFillCUDA(double* Point_Cloud,unsigned int* T,unsigned int* index_offset){
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
    __syncthreads();
    //En base a la malla referencial hallamos las demás superficies
    double xp,yp,zp;
    unsigned int count=0,offset,temp_vex,n_triangles_perThread=2*(n_beams-1);
    //creamos un Mesh array temporary
    unsigned int *T_per_thread;
    T_per_thread=(unsigned int*)malloc(5*n_triangles_perThread* 3 *sizeof(unsigned int));
    for (unsigned int i = 1; i < n_donuts; i++){
        //Analizamos cada vertice del tríangulo
        for (unsigned int j = 0; j < n_triangles_perThread; j++){
             //Analizamos el punto del vertice v0
            offset=(thid*n_triangles_perThread+j)*3;
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
                        T_per_thread[count*3+2]=temp_vex;
                        T_per_thread[count*3+1]=(T[offset+1]+i*n_points_perDonut);
                        T_per_thread[count*3]=(T[offset]+i*n_points_perDonut);
                        count++;
                    }
                }
            }
        }
    }
    index_offset[thid]=count;
    __syncthreads();
    unsigned int n_triangles_perDonut=n_AZBLK*2*(n_beams-1);
    //Ha este punto, cada hilo contiene una cierta cantidad de n triangulos que han de ser colocadas en el array original
    if(thid>0){
        offset=0;
        for(int z=0;z<thid;z++){
            offset+=index_offset[z];
        }
        offset*=3;
        //copy triangles 
        for(int z=0;z<count*3;z++){
            T[n_triangles_perDonut*3+offset+z]=T_per_thread[z];
        }
    }
    free(T_per_thread);
    __syncthreads();
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
    __syncthreads();
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

void Generate_surfaceGPU(double* Point_Cloud,unsigned int* T,unsigned int *pointer_n_triangles){
    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      CUDA        ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    cudaError_t cudaerr;
    double* Sphere_Cloud_dev;
    cudaMalloc((void**)(&Sphere_Cloud_dev), sizeof(double) * n_total_points * 3);
    //1st step. GenerateSphere
    GenerateSphereCUDA<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    cudaDeviceSynchronize();
    //2nd step. Overlap removing
    SupressOverlapCUDA<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev);
    double* Sphere_Cloud;
    Sphere_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    cudaerr=cudaMemcpy(Sphere_Cloud, Sphere_Cloud_dev, sizeof(double) *n_total_points * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to SphereCloud. CudaMalloc value=%i\n\r",cudaerr);
    unsigned int *OneDonutMesh_dev,*T_temp_dev,*count_array_dev,n_triangles_perDonut=n_AZBLK*2*(n_beams-1);
    cudaMalloc((void**)(&OneDonutMesh_dev), sizeof(unsigned int) * OneDonutFill_triangles * 3);
    cudaMalloc((void**)(&T_temp_dev), sizeof(unsigned int) * n_triangles_perDonut* (n_donuts-1) * 3);
    cudaMalloc((void**)(&count_array_dev), sizeof(unsigned int) * threadsPerBlock);
    //3rd step. First part of the ODF
    ODF_part1<<<numBlocks, threadsPerBlock >>> (Sphere_Cloud_dev,OneDonutMesh_dev,T_temp_dev,count_array_dev);
    unsigned int *index_offset_array_dev;
    cudaMalloc((void**)(&index_offset_array_dev), sizeof(unsigned int) * 1024);
    //4th step. eScan GPU
    cudaDeviceSynchronize();
    eScanGPU<<<numBlocks, threadsPerBlock/2 >>> (index_offset_array_dev,count_array_dev);
    //5th step. last part ODF
    cudaDeviceSynchronize();
    ODF_part2<<<numBlocks, threadsPerBlock >>> (OneDonutMesh_dev,T_temp_dev,count_array_dev,index_offset_array_dev);
    cudaDeviceSynchronize();
    unsigned int *T_temp;
    T_temp=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));
    cudaerr=cudaMemcpy(T_temp,OneDonutMesh_dev, sizeof(unsigned int) *OneDonutFill_triangles * 3, cudaMemcpyDeviceToHost);
    if (cudaerr != 0)	printf("ERROR copying to T_temp. CudaMalloc value=%i\n\r",cudaerr);
	/////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      CPU        ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //6th step. last Fill
    TwoandTri_Donut_Fill(Sphere_Cloud,&T_temp[OneDonutFill_triangles*3],&T_temp[(OneDonutFill_triangles+TwoDonutFill_triangles)*3],&T_temp[(OneDonutFill_triangles+TwoDonutFill_triangles+TriDonutFill_triangles)*3]);
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
    cudaFree(Sphere_Cloud_dev);
    cudaFree(OneDonutMesh_dev);
    cudaFree(T_temp_dev);
    cudaFree(count_array_dev);
    cudaFree(index_offset_array_dev);
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
    printf("CPU: %fms\n", ((double)(finishCPU - startCPU))/(double)CLOCKS_PER_SEC);
	
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
    printf("GPU: %fms\n", ((double)(finishGPU - startGPU))/(double)CLOCKS_PER_SEC);
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
    #define iter 500.0
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
    GenerateSphereCUDA<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    printf("GS time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;
    cudaDeviceSynchronize();

    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++)
    SupressOverlapCUDA<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev);
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    cudaerr=cudaMemcpy(Point_Cloud_gpu, Point_Cloud_dev, sizeof(double) *n_total_points * 3, cudaMemcpyDeviceToHost);
    printf("OR time:  %fms\n\r", timeGPU / iter);
    totalGPU+=timeGPU/iter;

    unsigned int *T_gpu,*OneMesh_dev,*OneMesh_temp_dev,n_triangles_perDonut=n_AZBLK*2*(n_beams-1);
    unsigned int *count_array_dev,*index_offset_array_dev;
    
    cudaMalloc((void**)(&OneMesh_dev), sizeof(unsigned int) * n_total_triangles * 3);
    cudaMalloc((void**)(&OneMesh_temp_dev), sizeof(unsigned int) * n_triangles_perDonut* (n_donuts-1) * 3);
    cudaMalloc((void**)(&count_array_dev), sizeof(unsigned int) * threadsPerBlock);
    cudaMalloc((void**)(&index_offset_array_dev), sizeof(unsigned int) * 1024);
    T_gpu=(unsigned int*)malloc(n_total_triangles * 3 *sizeof(unsigned int));

    //3rd step. First part of the ODF
    cudaDeviceSynchronize();
    cudaEventRecord(start_gpu);
    for (int z = 0; z < iter; z++){
        ODF_part1<<<numBlocks, threadsPerBlock >>> (Point_Cloud_dev,OneMesh_dev,OneMesh_temp_dev,count_array_dev);
        cudaDeviceSynchronize();
        //4th step. eScan GPU
        eScanGPU<<<numBlocks, threadsPerBlock/2 >>> (index_offset_array_dev,count_array_dev);
        cudaDeviceSynchronize();
        //5th step. last part ODF
        ODF_part2<<<numBlocks, threadsPerBlock >>> (OneMesh_dev,OneMesh_temp_dev,count_array_dev,index_offset_array_dev);
        cudaDeviceSynchronize();
    
    }
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);
	cudaEventElapsedTime(&timeGPU, start_gpu, stop_gpu);
    cudaerr=cudaMemcpy(T_gpu,OneMesh_dev, sizeof(unsigned int) *OneDonutFill_triangles * 3, cudaMemcpyDeviceToHost);
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