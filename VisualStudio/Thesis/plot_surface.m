close all; clear;clc
load("csv_data.mat")
T=readmatrix("One_donut_fill.csv")+1;
nube=readmatrix("Sphere_cloud.csv");
isequal(nube,Point_Cloud)
trimesh(T,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
temp=[OneDonutFill;TwoDonutFill;TriDonutFill;TriMiddleFill];
%trimesh(temp,nube(:,1),nube(:,2),nube(:,3))
