%% STEP 1.Generate overlapped structure
close all;clear;clc

[temp,~,~]=xlsread('dataXYZ1.csv');

n_elevation_points=16;
n_AZBLK=1024;
n_rotation=6;
points_per_donut=n_AZBLK*n_elevation_points;
points_data=cell(1,n_rotation);
for i=1:n_rotation
    points_data{i}=temp(points_per_donut*(i-1)+1:points_per_donut*i,:);
end
figure
hold
for i=1:n_rotation
    temp=points_data{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off
%% STEP 2.Convert structure into a sphere
points_data_spheric=points_data;
for i=1:n_rotation
    for j=1:points_per_donut
        temp=points_data{i}(j,:);
        mag_temp=norm(temp);
        if mag_temp~=0
            points_data_spheric{i}(j,:)=temp/mag_temp;
        end
    end
end

figure(2)
hold
for i=1:n_rotation
    temp=points_data_spheric{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

%% STEP 3.supress "redudant" data
alfa1=16.0620*pi/180;
alfa2=-15.3790*pi/180;
point1=[0.265203,0.964189];
point2=[0.754467,0.656332];
rot_angle=-acos(sum(point1.*point2)/norm(point1)/norm(point2));
%rot_angle=33.6°;

% syms rot_angle
% rot_matrix=[cos(rot_angle) -sin(rot_angle)  0;
%             sin(rot_angle) cos(rot_angle)   0;
%             0               0               1];


p1=[0,0,0];
p2=p1;
l1_param=[0,0,0];
l2_param=l1_param;
offset=pi/2;

%Let's set the boundary points of the donuts and every pair of equation lines
for i=1:n_rotation
    p1(i,:)=[cos(rot_angle*(i-1)+alfa1+offset),sin(rot_angle*(i-1)+alfa1+offset),0];
    p2(i,:)=[cos(rot_angle*(i-1)+alfa2+offset),sin(rot_angle*(i-1)+alfa2+offset),0];
    if i~=1
        l1_param(i,:)=[tan(rot_angle*(i-1)+offset),p1(i,1),p1(i,2)];
        l2_param(i,:)=[tan(rot_angle*(i-1)+offset),p2(i,1),p2(i,2)];
    end
end

%Verify relevant data
new_points_data_spheric=points_data_spheric;
%Donut2
i=2;
for j=1:length(new_points_data_spheric{i})
    x_data=new_points_data_spheric{i}(j,1);
    y_data=new_points_data_spheric{i}(j,2);
    z_data=new_points_data_spheric{i}(j,3);
    if x_data>p1(1,1) && x_data<p2(1,1)
        new_points_data_spheric{i}(j,:)=[0,0,0];
    end
end
%Donut 3 to 6
for i=3:n_rotation
    for j=1:length(new_points_data_spheric{i})
        x_data=new_points_data_spheric{i}(j,1);
        y_data=new_points_data_spheric{i}(j,2);
        z_data=new_points_data_spheric{i}(j,3);
        if x_data<p1(1,1) 
            y_temp=line_equation(x_data,l1_param(i-1,:));
            if y_data<y_temp
                new_points_data_spheric{i}(j,:)=[0,0,0];
            end
        elseif x_data>p2(1,1)
            y_temp=line_equation(x_data,l2_param(i-1,:));
            if y_data>y_temp
                new_points_data_spheric{i}(j,:)=[0,0,0];
            end
        else
            new_points_data_spheric{i}(j,:)=[0,0,0];
        end
    end
end
%%
figure(3)
hold
%Data a Carlos->temp=[];
for i=1:n_rotation
    temp=[new_points_data_spheric{i}];
    %Data a Carlos->temp=[temp;new_points_data_spheric{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3));
end
hold off

%% STEP 4.Generate structure without overlapped areas
new_points_data=points_data;

total_count=0;
for i=1:n_rotation
    count=0;
    for j=1:points_per_donut
        if isequal(new_points_data_spheric{i}(j,:),[0,0,0])
            new_points_data{i}(j,:)=[0,0,0];
            count=count+1;
        end
    end 
    fprintf("De la dona "+i+" se eliminaron "+count+" puntos. ("+count/163.84+"%%)\n");
    total_count=total_count+count;
end
fprintf("En total se eliminaron "+total_count+" puntos. ("+total_count/163.84/6+"%%)\n");
   
figure (4)
hold
%Data a Carlos->temp=[];
for i=1:n_rotation
    temp=[new_points_data{i}];
%   Data a Carlos->temp=[temp; new_points_data{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

%% Analyzing data (from step3)
temp=zeros(length(points_data{1}),2);
figure(5)
for i=1:n_rotation
    subplot(2,3,i)
    for j=1:16:length(points_data{i})
        for k=0:15
            if isequal(new_points_data{i}(j+k,:),[0 0 0])
                temp(j+k,:)=[0 0];
            else
                temp(j+k,:)=[(j-1)/16+1,16-k];
            end
        end
    end
    scatter(temp(:,1),temp(:,2))
    title('Donut '+string(i))
end

%% Triangulate (from step3)
figure(6);
hold
%Data a Carlos-> T=[];
%                k_puntos=16384;
for i=1:n_rotation
    temp=new_points_data_spheric{i};
    T=triangulate_v1(temp,n_elevation_points,n_AZBLK);
    %Data a Carlos->T=[T;triangulate_v1(temp,n_elevation_points,n_AZBLK)+(i-1)*k_puntos];
    temp=new_points_data{i};
    trimesh(T,temp(:,1),temp(:,2),temp(:,3))
end

%% Triangulate_v2
temp=[];
index=1;
for i=1:n_rotation
    for j=1:length(new_points_data{i})
        if ~isequal(new_points_data{i}(j,:),[0 0 0])
            temp(index,:)=new_points_data{i}(j,:);
            index=index+1;
        end
    end
end
P = delaunay(temp);
figure(7)
trisurf(P,temp(:,1),temp(:,2),temp(:,3),'FaceColor','cyan')
%% Making Surface Plots From Scatter Data
% How do you turn a collection of XYZ triplets into a surface plot? This is
% the most frequently asked 3D plotting question that I got when I was in
% Tech Support.
%% Load the data
% load seamount
% who -file seamount
%[temp,~,~]=xlsread('dataXYZ1.csv');
x=temp(:,1);
y=temp(:,2);
z=temp(:,3);
%%
% The problem is that the data is made up of individual (x,y,z)
% measurements. It isn't laid out on a rectilinear grid, which is what the
% SURF command expects. A simple plot command isn't very useful.
figure (8)
plot3(x,y,z,'.-')
%% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.
tri = delaunay(x,y);
plot(x,y,'.')
%%
% How many triangles are there?
[r,c] = size(tri);
disp(r)
%% Plot it with TRISURF
h = trisurf(tri, x, y, z);
axis vis3d
%% Clean it up
axis off
l = light('Position',[-50 -15 29])
set(gca,'CameraPosition',[208 -50 7687])
lighting phong
shading interp
colorbar EastOutside