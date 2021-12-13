%% STEP1. Generate overlapped structure
close all;clear;clc
[temp,~,~]=xlsread('dataXYZ1.csv');

n_rays=16;
n_AZBLK=1024;
n_donuts=6;
n_points=n_AZBLK*n_rays;
point_cloud=cell(1,n_donuts);
for i=1:n_donuts
    point_cloud{i}=temp(n_points*(i-1)+1:n_points*i,:);
end
figure (1)
hold
for i=1:n_donuts
    temp=point_cloud{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
title('Nube de puntos sin procesar')
%% STEP 2.Convert structure into a sphere
point_cloud_spheric=point_cloud;
for i=1:n_donuts
    for j=1:n_points
        temp=point_cloud{i}(j,:);
        if ~isequal(temp,[0 0 0])
            point_cloud_spheric{i}(j,:)=temp/norm(temp);
        end
    end
end
figure (2)
hold
for i=1:n_donuts
    temp=point_cloud_spheric{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
% %% See how kth-azimuth is plotted step by step per donut
% actual_number_fig=2;
% for i=1:n_donut
%     figure (i+1+actual_number_fig)
%     hold
%     temp=point_cloud_spheric{i};
%     scatter3(temp(:,1),temp(:,2),temp(:,3))
%     offset=1024/16;
%     count=0;
%     for azimuth=offset:offset:n_AZBLK*1/2
%         index=16*azimuth-15;
%         temp=point_cloud_spheric{i}((index:index+15),:);
%         scatter3(temp(:,1),temp(:,2),temp(:,3),150,'filled')
%         count=count+1;
%     end
%     title("Nube de puntos esférica traslapada "+num2str(i))
%     legend
% end
%% STEP3. supress "redundant" data
%set parameters
rot_angle=-33.53706667*pi/180;

alfa1=15.379*pi/180;
alfa2=-16.062*pi/180;

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
new_point_cloud_spheric=point_cloud_spheric;

%Donut2
i=2;
del_points=[0 0 0 0 0 0];
for j=1:length(new_point_cloud_spheric{i})
    x_data=new_point_cloud_spheric{i}(j,1);
    y_data=new_point_cloud_spheric{i}(j,2);
    z_data=new_point_cloud_spheric{i}(j,3);
    if x_data<=p1(1,1) && x_data>=p2(1,1)
        new_point_cloud_spheric{i}(j,:)=[0,0,0];
        del_points(i)=del_points(i)+1;
    end
end
%Donut 3 to 6
for i=3:n_donuts
    for j=1:length(new_point_cloud_spheric{i})
        x_data=new_point_cloud_spheric{i}(j,1);
        y_data=new_point_cloud_spheric{i}(j,2);
        z_data=new_point_cloud_spheric{i}(j,3);
        if x_data>p1(1,1) 
            y_temp=line_equation(x_data,l1_param(i-1,:));
            if y_data>=y_temp
                new_point_cloud_spheric{i}(j,:)=[0,0,0];
                del_points(i)=del_points(i)+1;
            end
        elseif x_data<p2(1,1)
            y_temp=line_equation(x_data,l2_param(i-1,:));
            if y_data<=y_temp
                new_point_cloud_spheric{i}(j,:)=[0,0,0];
                del_points(i)=del_points(i)+1;
            end
        else
            new_point_cloud_spheric{i}(j,:)=[0,0,0];
            del_points(i)=del_points(i)+1;
        end
    end
end
%Data a Carlos:
spheric_point_cloud_without_overlapped=zeros(n_points*n_donuts,3);

figure(3)
hold
for i=1:n_donuts
    temp=[new_point_cloud_spheric{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3));
    %Data a Carlos
    spheric_point_cloud_without_overlapped((i-1)*n_points+1:i*n_points,:)=temp;
end
hold off
%% STEP 4.Generate structure without overlapped areas
new_point_cloud=point_cloud;
total_count=0;
for i=1:n_donuts
    count=0;
    for j=1:n_points
        if isequal(new_point_cloud_spheric{i}(j,:),[0,0,0])
            new_point_cloud{i}(j,:)=[0,0,0];
            count=count+1;
        end
    end 
    fprintf("De la dona "+i+" se eliminaron "+del_points(i)+" puntos ("+del_points(i)/163.84+"%%)"+" y se perdieron "+(count-del_points(i))+" puntos ("+(count-del_points(i))/163.84+"%%)\n");
    total_del=total_count+del_points(i);
    total_count=total_count+count;
end
fprintf("En total se eliminaron "+total_del+" puntos. ("+total_count/163.84/6+"%%)\n");
fprintf("En total se perdieron "+(total_count-total_del)+" puntos. ("+(total_count-total_del)/163.84/6+"%%)\n");
fprintf("Quedaron "+(16384*6-total_count)+" puntos");

%Data a Carlos
points_cloud_without_overlapped=zeros(n_points*n_donuts,3);
figure (4)
hold
for i=1:n_donuts
    temp=[new_point_cloud{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3))
    %Data a Carlos
    points_cloud_without_overlapped((i-1)*n_points+1:i*n_points,:)=temp;
end
ptCloud=pointCloud(points_cloud_without_overlapped);
pcshow(ptCloud)
hold off
%% Delaunay 
%(Solo para funciona para el exterior)
T=delaunay(spheric_point_cloud_without_overlapped);
temp=points_cloud_without_overlapped;
figure(5)
trimesh(T,temp(:,1),temp(:,2),temp(:,3))
%     


%% Triangulate_v1

%Data a Carlos
Triangle_mesh=[];

figure(6);
hold
for i=1:1
    temp=new_point_cloud_spheric{i};
    T=triangulate_v1(temp,n_rays,n_AZBLK);
    temp=new_point_cloud{i};
    trimesh(T,temp(:,1),temp(:,2),temp(:,3))
    %Data a Carlos
    Triangle_mesh=[Triangle_mesh;T+(i-1)*n_points];    
end
hold off
% %% Analyzing data (from step3)
% temp=zeros(n_points,2);
% figure(6)
% for i=1:n_donuts
%     subplot(2,3,i)
%     for j=1:16:n_points
%         for k=0:15
%             if isequal(new_point_cloud{i}(j+k,:),[0 0 0])
%                 temp(j+k,:)=[0 0];
%             else
%                 temp(j+k,:)=[(j-1)/16+1,16-k];
%             end
%         end
%     end
%     scatter(temp(:,1),temp(:,2))
%     title('Donut '+string(i))
% end
