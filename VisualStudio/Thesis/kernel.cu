﻿#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "math.h"
#include "time.h"
#ifndef MPI//definimos el valor de PI
#define MPI 3.14159265358979323846
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
double beam_altitude_angles[n_beams]= {15.379,13.236,11.128,9.03,6.941,4.878,2.788,0.705,-1.454,-3.448,-5.518,-7.601,-9.697,-11.789,-13.914,-16.062};
double beam_azimuth_angles[n_beams] = { -1.24, -1.2145, -1.1889, -1.1634, -1.1379, -1.1123, -1.0868, -1.0613, -1.0357, -1.0102, -0.98467, -0.95913, -0.9336, -0.90807, -0.88253, -0.857 };
/*Cantidad puntos de la nube*/
#define n_points_perDonut (n_AZBLK*n_beams)
/*Angulo entre azimuths*/
#define angle_between_azimuths (-2*MPI/n_AZBLK)
/*Ángulo de rotacion del motor*/
#define rot_angle (-33.53706667 * MPI / 180)
/*Cantidad de Donuts en función del ángulo de rotación*/
unsigned char n_donuts = (unsigned char)ceil(-MPI/rot_angle);
#define n_total_points (n_points_perDonut*n_donuts)
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
    double x,y,z,R= Radius_sphere;
    for (int i = 0; i < n_beams; i++) {
        /*ubicamos el punto del azimuth*/
        /*x = r * cos(beam_altitude_angles[i] - mpi / 2);
        y = 0;
        z = r * sin(beam_altitude_angles[i] - mpi / 2);




        temp = rot_x_axis(beam_azimuth_angles[i]) * [x, y, z]';
            point_cloud(j, 1) = temp(1);
        point_cloud(j, 2) = temp(2);
        point_cloud(j, 3) = temp(3);*/
    }


}
/*----------------------------------------------------------------------------*/

int main()
{  
    /** 
    * Nuestro sistema de referencia será : Eje Z será el eje de giro del motor.
    *                                      Eje X será el eje de la Donut referencial
    */
    /*Creamos la matrix de rotación del eje del motor (eje Z)*/
    double rot_matrix[9] = { cos(rot_angle),-sin(rot_angle),0,sin(rot_angle),cos(rot_angle) ,0,0,0,1 };
    /*Allocate memory*/
    double* Point_Cloud;
    Point_Cloud= (double*)malloc(n_total_points*sizeof(double));











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
