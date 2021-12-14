#include <stdio.h>
#include <stdint.h>
#include "math.h"
#ifndef MPI//definimos el valor de PI
#define MPI 3.14159265358979323846
#endif
//--------------------------------------------//
//-----Definimos los parámetros del LIDAR-----//
//--------------------------------------------//
//Cantidad de rayos por azimuth
#define n_beams 16
//Cantidad de azimuths por Donut
#define n_AZBLK 1024
//Los ángulos por defecto de fabrica del LIDAR
double beam_altitude_angles[n_beams]= {15.379,13.236,11.128,9.03,6.941,4.878,2.788,0.705,-1.454,-3.448,-5.518,-7.601,-9.697,-11.789,-13.914,-16.062};
double beam_azimuth_angles[n_beams] = { -1.24, -1.2145, -1.1889, -1.1634, -1.1379, -1.1123, -1.0868, -1.0613, -1.0357, -1.0102, -0.98467, -0.95913, -0.9336, -0.90807, -0.88253, -0.857 };
//Cantidad puntos de la nube
#define n_points (n_AZBLK*n_beams)
//Angulo entre azimuths
#define angle_between_azimuths (-2*MPI/n_AZBLK)
//Ángulo de rotacion del motor
#define rot_angle (-33.53706667 * MPI / 180);
//Cantidad de Donuts en función del ángulo de rotación
unsigned char n_donuts = (unsigned char)(-180 + beam_altitude_angles[0] - beam_altitude_angles[15])*MPI/180/rot_angle;

//--------------------------------------------//
//-----     Funciones de conversión      -----//
//--------------------------------------------//
void rad2deg(double *value) {
    *value = (*value) * 180 / MPI;
}
void deg2rad(double *value) {
    value[0] = value[0] * MPI / 180;
}

//--------------------------------------------//
//----- Funciones de matrices de rotación-----//
//--------------------------------------------//
void rot_matrix() {
}

void rot_x_axis() {
}

void rot_y_axis() {
}

void rot_z_axis() {
}
int main()
{  
     
    printf("Hola mundo %f %f\n",MPI,beam_altitude_angles[2]);
    return 0;
}
