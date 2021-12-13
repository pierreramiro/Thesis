clear;close all;clc
%characteristic of LIDAR
    %Angle elevation range; [-16.677; 16,545]
        %In this particular device we have 16 laser beams/points
        %Also, every azimuth block has an angle offset (-1.24°)
        %AZBLK=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062]
    %About encode range: [0,90111]
    %In mode 1023 azimuth. Every azimuth increments in 88 ticks
        %1 tick = 360/90112= 0.003995°.     Aprox -> 4°
        %88 tick = 360*88/90112 = 0.3516°
    %In total, we obtain 1024*16 = 16384 points in every LIDAR scan
      
%Set Local variables:
AZBLK=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062];
azimuth_offset=-1.24*(pi/180); 
n_AZBLK=1024;
ticks_per_azimuth=90112/n_AZBLK;
ang_diff=2*pi/n_AZBLK;

%Interception with synthetic spheric room
    %sphere with radius=1
n_rotation=6;
temp=zeros(n_AZBLK*length(AZBLK),3);
points_data=cell(1,n_rotation);
for i=1:n_rotation %pre-allocate
    points_data{1,i}=temp;
end

%Rotation -> 31°
rot_angle=-31*pi/180;
rot_matrix=[cos(rot_angle) -sin(rot_angle)  0;
            sin(rot_angle) cos(rot_angle)   0;
            0               0               1];
R=1;%radius sphere
for i=1:n_rotation
    for j=1:n_AZBLK
        x=sin(ang_diff*(j-1)+azimuth_offset);
        z=-cos(ang_diff*(j-1)+azimuth_offset);
        for k=1:length(AZBLK)
            y=tan(AZBLK(k)*pi/180);
            new_xyz=rot_matrix^(i-1)*[x;y;z];
            %sphere
            [thetha,phi,~]=cart2sph(new_xyz(1),new_xyz(2),new_xyz(3));
            [temp_x,temp_y,temp_z]=sph2cart(thetha,phi,R);           
            temp(length(AZBLK)*(j-1)+k,:)=[temp_x,temp_y,temp_z];
        end
    end
    points_data{1,i}=temp;    
end

figure
hold
for i=1:6%n_rotation
    temp=points_data{1,i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off
%%
syms x 
elev_max=AZBLK(1)*(pi/180);%elevation angle
elev_min=AZBLK(end)*(pi/180);
y=R*sin(elev_max);
eqns=rot_matrix(1:2,1:2)^2*[x;y];
eqn=eqns(2)==R*sin(elev_min);
x=eval(solve(eqn,x));
y=R*sin(elev_min);
z=sqrt(R*R-(x*x+y*y));
up_point=[x,y,z];
down_point=[x,y,-z];

temp=surfaces_allowed_v1(points_data,up_point,down_point);
temp(temp(:,1)==0,:) = [] ;
n_points=length(temp);
for i=1:4
    for j=1:n_points
        temp=[temp;(rot_matrix^(i-1)*temp(j,:)')'];
    end
end
rot_180_matrix=[-1 0 0;
                0 -1 0;
                0  0 1];
temp=[temp;(rot_180_matrix*temp')'];

%nueva area
temp_particular=particular_surfaces_allowed_v1(points_data,up_point,down_point);
temp=[temp;temp_particular];

temp_last_surfaces=last_surfaces_v1(points_data);
temp=[temp;temp_last_surfaces];

figure
scatter3(temp(:,1),temp(:,2),temp(:,3),[],[1 1 0])
%scatter(temp(:,1),temp(:,2))