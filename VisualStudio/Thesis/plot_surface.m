close all; clear;clc
load("csv_data.mat")
T=[];
temp=readmatrix("One_donut_fill.csv")+1;
T=[T;temp];
temp=readmatrix("Two_donut_fill.csv")+1;
T=[T;temp];
temp=readmatrix("Tri_donut_fill.csv")+1;
T=[T;temp];
temp=readmatrix("Mid_donut_fill.csv")+1;
T=[T;temp];
Point_Cloud=readmatrix("Sphere_cloud.csv");
trimesh(T,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
%%
figure
Point_Cloud=readmatrix("dataMina.csv");
%T=readmatrix("Triangle_mesh.csv")+1;
trimesh(T,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
