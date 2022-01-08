#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "math.h"
#include "time.h"
#ifndef MPI//definimos el valor de PI
#define MPI 3.14159265358979323846
#define MPI_2 1.57079632679489661923
#define D180_MPI 0.017453293
#endif
//#include "mkl.h"
/*----------------------------------------------------------------------------*/
/**
 * Definimos los parámetros del LIDAR
 * Se puede usar typedef para juntar los parámetros del LIDAR
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
#define n_donuts (unsigned int) ceil(-MPI/rot_angle)
#define n_total_points (unsigned int)(n_donuts*n_points_perDonut)
/*Cantidad de triángulos*/
#define OneDonutFill_triangles          108226
#define TwoandTriDonutFill_triangles    10036

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
    double x,y,y_temp,x1,y1,xn,yn,m;
    //Set vertical limits
    double X_L1=Point_Cloud[0],X_Ln=Point_Cloud[15*3+0];
    //Array which contains lineś parameters
    double L[5*(n_donuts-2)];
    //Declare parameters from Donut_1 to Donut_n-2
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
    //Supress redundant data only for Donut 1
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
            y_temp=y_temp*(1-2*left_side);
            y=y*(1-2*left_side);
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
    printf("numero de triangulos: %d\n",count+n_triangles_perDonut);    
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
void get_tripivot(unsigned int *v0_pointer,unsigned int *v1_pointer,double *point,unsigned int sector,unsigned int i,unsigned int k_beam){
    //Definimos las variables
    unsigned int offset,k_azimuth,v0_temp,v1_temp,v_temp;
    double y_data,z_data,theta,rot_theta,alfa;
    double *rot_point;
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
    rot_theta=-acos(rot_point[1]/Radius_sphere);
    if (theta*(1-2*((sector>>1)&0x1)) >= rot_theta*(1-2*((sector>>1)&0x1)))
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
    alfa=alfa-2*pi;
    //Con lo anterior, nos hemos asegurado que el alfa sea siempre negativo
    //Calculamos el azimuth que corresponde al alfa
    k_azimuth=(unsigned int)((alfa-beam_azimuth_angles[k_beam])/angle_between_azimuths);
    v0_temp=(k_azimuth)*n_beams+k_beam;
    //Enmascaramos
    v0_temp=v0_temp&mask;
    v1_temp=(v0_temp+n_beams)&mask;
    //Realizamos offset y sentido
    v0_temp+=offset;
    v1_temp+=offset;
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
        v_temp=v0_temp;
        v0_temp=v1_temp;
        v1_temp=v_temp;
    }
    *v0_pointer=v0_temp;
    *v1_pointer=v1_temp;
                
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
#define NO_VEX_FOUND -1 
void TwoandTri_Donut_Fill(double* Point_Cloud,unsigned int* T){
    int paso;
    unsigned int n_tripivot,init_index,v0,v1,v2,k_beam;
    double x_point,y_point,z_point;
    unsigned int *last_Tripivots;
    bool vex_found;
    //A partir de la segunda Donut generamos las superficies
    for (int i = 0; i < n_donuts; i++){
        for (unsigned int sector = 0; sector < 4; sector++){
            n_tripivot=0;
            /*Creamos los tripivots*/
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
                paso=n_beams*(1-2*(sector&0x1));
                //Creamos bandera para saber si se halló un vértice para crear el tripivot
                vex_found=false;
                for (unsigned int count = 0; count < n_AZBLK; count++){
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
                        //Realizamos un conteo de los tripivot
                        n_tripivot++;
                        break;
                    }
                    //Realizamos el paso correspondiente a init_index
                    init_index=(init_index+paso)&mask;
                }
                //Evaluamos el caso que no se consiguió el vértice
                if (!vex_found){
                    //Evaluamos para los siguientes beams
                    continue
                }
                /*Procedemos a hallar los demás vértices del tripivot*/
                //Dependiendo del sector estamos más cerca del beam=0 o el beam=15
                k_beam=(n_beams-1)*(1-(sector>>1)&0x1);
                ///Hallamos los otros vértices
                get_tripivot(&v0,&v1,&Point_Cloud[v2],sector,i,k_beam);
                //continuamos aqui:




                
            }
            
        }
        
        
    }
    
}

/*----------------------------------------------------------------------------*/
int main()
{
    /**
    * Nuestro sistema de referencia será : Eje Z será el eje de giro del motor.
    *                                      Eje X será el eje de la Donut referencial
    */
    /*Allocate memory*/
    double* Point_Cloud;
    Point_Cloud = (double*)malloc(n_total_points * 3 *sizeof(double));
    unsigned int *T_ODF, *T_TTDF;
    T_ODF=(unsigned int*)malloc(OneDonutFill_triangles*3*sizeof(unsigned int));
    T_TTDF=(unsigned int*)malloc(TwoandTriDonutFill_triangles*3*sizeof(unsigned int));
    Generate_sphere(Point_Cloud);
    Supress_redundant_data(Point_Cloud);
    One_Donut_Fill(Point_Cloud,T_ODF);
    TwoandTri_Donut_Fill(Point_Cloud,T_TTDF);
    /*Escribimos la data obtenida en un archivo csv*/

    FILE* archivo;
    archivo = fopen("Sphere_cloud.csv", "w+");
    fprintf(archivo, "X, Y, Z\n");
    for (unsigned int i=0; i < n_total_points; i++) {
        fprintf(archivo,"%.4f, %.4f, %.4f\n", Point_Cloud[i*3+0], Point_Cloud[i * 3 + 1], Point_Cloud[i * 3 + 2]);
    }
    fclose(archivo);

    archivo = fopen("One_donut_fill.csv", "w+");
    fprintf(archivo, "V1, V2, V3\n");
    for (unsigned int i=0; i < n_AZBLK*2*(n_beams-1)*n_donuts; i++) {
        fprintf(archivo,"%d, %d, %d\n", T_ODF[i*3+0], T_ODF[i * 3 + 1], T_ODF[i * 3 + 2]);
    }
    fclose(archivo);


    //------------------------------//
    /*Testing function*/
    double A[12] = { 1,2,3,4,5,6,7,8,9,10,11,12 };
    double B[12] = { 10,7,4,1,11,8,5,2,12,9,6,3};
    double C[16];
    mult_matrix(A, 1, 3, B, 4, C);
    for (int z=0; z < 1; z++) {
        for (int w=0; w < 4; w++) {
            printf("%.3f\t", C[z*4+w]);
        }
        printf("\n");
    }
    //-----------------------------// 
    /*Testing rot axis*/
    double temp[3] = { 2,4,5 };
    printf("Antes:\n%.3f\t%.3f\t%.3f\n", temp[0], temp[1], temp[2]);
    rot_x_axis(temp,MPI);
    printf("Despues:\n%.3f\t%.3f\t%.3f\n", temp[0], temp[1], temp[2]);
    //----------------------------//
    /*Testing n_donuts value*/
    printf("Hola mundo %d %f\n", n_donuts,beam_altitude_angles[2]);
    /* for (int i=0;i<n_beams; i++) {
        printf("%f\n", beam_altitude_angles[i]);
    }*/
    return 0;
}
