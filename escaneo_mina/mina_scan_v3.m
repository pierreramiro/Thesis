close all;clear;clc
%% Set LIDAR's parameters
%Angle elevation range;
    AZBLK_angle=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062]*pi/180;
    n_rays=length(AZBLK_angle);
%Offset angle of AZBLK's points
    %AZBLK_offset=regexprep(num2str(linspace(-1.24,-0.857,16)),'\s+',',')
    AZBLK_offset=[-1.24,-1.2145,-1.1889,-1.1634,-1.1379,-1.1123,-1.0868,-1.0613,-1.0357,-1.0102,-0.98467,-0.95913,-0.9336,-0.90807,-0.88253,-0.857]*pi/180;
%Encoder range: [0,90111] ticks
    total_ticks=90112;
%Mode 1023 azimuth. Every azimuth increments in 88 ticks
    %1 tick = 360/90112= 0.003995°.     Aprox -> 0.004°
    %88 tick = 360*88/90112 = 0.3516°
    n_AZBLK=1024;
    %It should be noted that, we obtain 1024*16 = 16384 points in every LIDAR scan
    n_points=n_AZBLK*n_rays;
    ticks_between_azimuths=total_ticks/n_AZBLK;
    angle_between_azimuths=-2*pi*(ticks_between_azimuths/total_ticks);
                         %=-2*pi/n_AZBLK
%Angle rotation,n° donuts,angular velocity,time per scan.
    n_donuts=6;
    rot_angle=-33.53706667*pi/180;
    rot_matrix=[cos(rot_angle) -sin(rot_angle)  0;
                sin(rot_angle) cos(rot_angle)   0;
                0               0               1];
%pre-allocate spheric point clouds
    temp=zeros(n_AZBLK*n_rays,3);
    spheric_point_cloud=cell(1,n_donuts);
    for i=1:n_donuts
        spheric_point_cloud{i}=temp;
    end
%% Generate synthetic sphere
R=1;%radius sphere
for i=1:n_donuts
    for j=1:n_AZBLK
        for k=1:n_rays
            x=R*sin(AZBLK_angle(k));
            z=-R*cos(AZBLK_angle(k))*cos(angle_between_azimuths*(j-1)+AZBLK_offset(k));
            y=R*cos(AZBLK_angle(k))*sin(angle_between_azimuths*(j-1)+AZBLK_offset(k));
            temp((j-1)*n_rays+k,:)=(rot_matrix^(i-1)*[x;y;z])';
        end
    end
    spheric_point_cloud{i}=temp;    
end
figure (1)
hold
for i=1:n_donuts
    temp=spheric_point_cloud{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off
title('Nube de puntos de 6 donas, con traslape')
%% Supress "redundant" data
%set parameters
alfa1=AZBLK_angle(1);
alfa2=AZBLK_angle(end);

p1=zeros(n_donuts,3);
p2=p1;
l1_param=[0,0,0];
l2_param=l1_param;
offset=-pi/2;
%Set boundary points of the donuts and every pair of equation lines
for i=1:n_donuts
    angle_temp=rot_angle*(i-1)+alfa1+offset;
    p1(i,:)=[cos(angle_temp),sin(angle_temp),0];
    angle_temp=rot_angle*(i-1)+alfa2+offset;
    p2(i,:)=[cos(angle_temp),sin(angle_temp),0];
    if i~=1
        l1_param(i,:)=[tan(rot_angle*(i-1)+offset),p1(i,1),p1(i,2)];
        l2_param(i,:)=[tan(rot_angle*(i-1)+offset),p2(i,1),p2(i,2)];
    end
end
%supress data
spheric_point_cloud_processed=spheric_point_cloud;

%Donut2
i=2;
for j=1:length(spheric_point_cloud_processed{i})
    x_data=spheric_point_cloud_processed{i}(j,1);
    y_data=spheric_point_cloud_processed{i}(j,2);
    z_data=spheric_point_cloud_processed{i}(j,3);
    if x_data<=p1(1,1) && x_data>=p2(1,1)
        spheric_point_cloud_processed{i}(j,:)=[0,0,0];
    end
end
%Donut 3 to 6
for i=3:n_donuts
    for j=1:length(spheric_point_cloud_processed{i})
        x_data=spheric_point_cloud_processed{i}(j,1);
        y_data=spheric_point_cloud_processed{i}(j,2);
        z_data=spheric_point_cloud_processed{i}(j,3);
        if x_data>p1(1,1) 
            y_temp=line_equation(x_data,l1_param(i-1,:));
            if y_data>=y_temp
                spheric_point_cloud_processed{i}(j,:)=[0,0,0];
            end
        elseif x_data<p2(1,1)
            y_temp=line_equation(x_data,l2_param(i-1,:));
            if y_data<=y_temp
                spheric_point_cloud_processed{i}(j,:)=[0,0,0];
            end
        else
            spheric_point_cloud_processed{i}(j,:)=[0,0,0];
        end
    end
end

spheric_point_cloud_without_overlapped=zeros(n_points*n_donuts,3);
figure(2)
hold
for i=1:n_donuts
    temp=[spheric_point_cloud_processed{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3));
    spheric_point_cloud_without_overlapped((i-1)*n_points+1:i*n_points,:)=temp;
end
hold off
title('Nube de puntos de 6 donas, sin traslape')
%% Triangulate with MyRobustCrust
[Triangle_mesh,~]=MyRobustCrust(spheric_point_cloud_without_overlapped);
temp=spheric_point_cloud_without_overlapped;
figure (3)
trimesh(Triangle_mesh,temp(:,1),temp(:,2),temp(:,3));
title('Superficie tridimensional del modelo esférico')

%% Reconstruction
[point_cloud_without_overlapped,~,~]=xlsread('dataXYZ1.csv');
Triangle_mesh_mina=Triangle_mesh;
for i=length(Triangle_mesh_mina):-1:1
    v=Triangle_mesh_mina(i,1);
    if isequal(point_cloud_without_overlapped(v,:),[0 0 0])
        Triangle_mesh_mina(i,:)=[];
        continue
    end
    v=Triangle_mesh_mina(i,2);
    if isequal(point_cloud_without_overlapped(v,:),[0 0 0])
        Triangle_mesh_mina(i,:)=[];
        continue
    end
    v=Triangle_mesh_mina(i,3);
    if isequal(point_cloud_without_overlapped(v,:),[0 0 0])
        Triangle_mesh_mina(i,:)=[];
        continue
    end    
end

figure (4)
temp=spheric_point_cloud_without_overlapped;
trimesh(Triangle_mesh_mina,temp(:,1),temp(:,2),temp(:,3));
title('Superficie tridimensional del modelo esférico')

%%
plot_mesh(point_cloud_without_overlapped,Triangle_mesh_mina);

figure(6)
temp=point_cloud_without_overlapped;
trimesh(Triangle_mesh_mina,temp(:,1),temp(:,2),temp(:,3));
%% Aditional. Triangulation with Delaunay
T=delaunay(spheric_point_cloud_without_overlapped);
figure(7)
temp=spheric_point_cloud_without_overlapped;
trimesh(T,temp(:,1),temp(:,2),temp(:,3));