clc;clear; close all
n_faces=16*6-5;
[x,y,z]=sphere(n_faces);
scatter3(x(:,1),y(:,1),z(:,1)) 
hold
for i=2:n_faces
    scatter3(x(:,i),y(:,i),z(:,i)) 
end