clear;clc
cloud=readmatrix("dataXYZ1.csv");
pointsPerDonut=1024*16;
Dona1=cloud(1:pointsPerDonut,:);
Dona2=cloud(pointsPerDonut+1:pointsPerDonut*2,:);
Dona3=cloud(pointsPerDonut*2+1:pointsPerDonut*3,:);
Dona4=cloud(pointsPerDonut*3+1:pointsPerDonut*4,:);
Dona5=cloud(pointsPerDonut*4+1:pointsPerDonut*5,:);
Dona6=cloud(pointsPerDonut*5+1:pointsPerDonut*6,:);
writematrix(Dona1,"Donas/Don1.csv")
writematrix(Dona2,"Donas/Don2.csv")
writematrix(Dona3,"Donas/Don3.csv")
writematrix(Dona4,"Donas/Don4.csv")
writematrix(Dona5,"Donas/Don5.csv")
writematrix(Dona6,"Donas/Don6.csv")
%noblecilla.ae@pucp.edu.pe