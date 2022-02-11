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
T=readmatrix("Triangle_mesh.csv")+1;
trimesh(T,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
%%
% figure
% Point_Cloud=readmatrix("esfera.csv");
% T=readmatrix("T.csv")+1;
% trimesh(T,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
% figure
% scatter3(Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
%%
clc;close
T=[];
temp=readmatrix("CUDAOneMesh.csv")+1;
T=[T;temp];
Point_Cloud=readmatrix("CUDASphere_cloud.csv");
trimesh(T,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
%%
scatter3(Point_Cloud(1:16,1),Point_Cloud(1:16,2),Point_Cloud(1:16,3))
hold
for i=1024*5+2:1024*6
scatter3(Point_Cloud((i-1)*16+1:i*16,1),Point_Cloud((i-1)*16+1:i*16,2),Point_Cloud((i-1)*16+1:i*16,3))
pause(0.00001);
end
%%
temp=readmatrix("CUDAOneMesh.csv");
temp2=readmatrix("One_donut_fill.csv");

for i=1:length(temp)
    if ~isequal(temp(i,:),temp2(i,:))
        display(i);
        display(temp(i,:));
        display(temp2(i,:));
        lost_values(i,:)=temp(i,:);
        break;
    end
end
%%
temp3=readmatrix("CUDASphere_cloud.csv");
%%
thid=0;
for i=1:16*3
index_nBeams(i)=fix(thid/3);
index_XYZ(i)=fix(thid-index_nBeams(i)*3);
thid=thid+1;
end
index_nBeams
index_XYZ