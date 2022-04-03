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
%Tiene que ser menor al la suma de angulos de elevacion, pero mayor a 30
%para cumplir con el barrido total
rot_angle=-31.4*pi/180;
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
for i=1:n_rotation
    temp=points_data{1,i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

%In total, we have (4+1+1)*2 areas without overlap
%If we look at the x-Y plane it seems like triangles
%There are 4*2 triangles, which surrounded the circle, with the same area
%And other two pair of triangles with different areas
%Understanding this, lets found every 3 points for each triangle
    %It should be noted that the 2 points of each triangle that surround the 
    %circle has the next restriction:
        %Considering the 1st point as the point that has more angle in-plane X-Y
        %The 1st point is the point related to the previous azimuth with a minimun
        %elevation
        %The 2nd point is the point related to the next azimuth with a maximum 
        %elevation
    %To obtain the 3rd point is necessary to found the eq line with the
    %rotation angle 
%In sum, lets obtain the X component from the data points, then determinate 
%in which area it belongs each point with the Y component. In adittion, we 
%can use the restriction with the Z component of 3rd point before evaluate 
%the area.

p1=[ 0 0 0];
p2=[ 0 0 0];
p3=[ 0 0 0];
syms x
for i=1:4
%Area1:4
    alfa1=AZBLK(end)*pi/180+rot_angle*(i-1);
    alfa2=AZBLK(1)*pi/180+rot_angle*(i+1);
    p1(i,:)=[cos(alfa1) sin(alfa1) 0];
    p2(i,:)=[cos(alfa2) sin(alfa2) 0];
    l1(i)=tan(rot_angle*(i-1))*(x-p1(i,1))+p1(i,2);
    l2(i)=tan(rot_angle*(i+1))*(x-p2(i,1))+p2(i,2);
    temp_x=eval(solve(l1(i)==l2(i),x));
    temp_y=subs(l1(i),x,temp_x);
    p3(i,:)=[temp_x temp_y,sqrt(R*R-temp_x*temp_x-temp_y*temp_y)];
end
%Area5
i=5;
alfa1=AZBLK(end)*pi/180+rot_angle*(i-1);
alfa2=(-180-AZBLK(end))*pi/180;%it changes
p1(i,:)=[cos(alfa1) sin(alfa1) 0];
p2(i,:)=[cos(alfa2) sin(alfa2) 0];
l1(i)=tan(rot_angle*(i-1))*(x-p1(i,1))+p1(i,2);
l2(i)=p2(i,2);%it changes
temp_x=eval(solve(l1(i)==l2(i),x));
temp_y=subs(l1(i),x,temp_x);
p3(i,:)=[temp_x temp_y,sqrt(R*R-temp_x*temp_x-temp_y*temp_y)];

%Area6
i=6;
alfa1=AZBLK(end)*pi/180+rot_angle*(i-1);
alfa2=(-180-AZBLK(end))*pi/180+rot_angle;%it changes
p1(i,:)=[cos(alfa1) sin(alfa1) 0];
p2(i,:)=[cos(alfa2) sin(alfa2) 0];
l1(i)=tan(rot_angle*(i-1))*(x-p1(i,1))+p1(i,2);
l2(i)=tan(rot_angle)*(x-p2(i,1))+p2(i,2);%it changes
temp_x=eval(solve(l1(i)==l2(i),x));
temp_y=subs(l1(i),x,temp_x);
p3(i,:)=[temp_x temp_y,sqrt(R*R-temp_x*temp_x-temp_y*temp_y)];
%notar que ang
%acos(dot(p3(5,1:2),p3(6,1:2))/(norm(p3(5,1:2))^2))*180/pi->rot_angle
%es decir, el angulo entre los ultimos puntos (en el plano X-Y)) es igual 
%al angulo de rotación.
%Area7to10
for i=7:10
    alfa1=(-180-AZBLK(1))*pi/180+(i-7)*rot_angle;
    alfa2=(-180-AZBLK(end))*pi/180+(i-5)*rot_angle;
    p1(i,:)=[cos(alfa1) sin(alfa1) 0];
    p2(i,:)=[cos(alfa2) sin(alfa2) 0];
    l1(i)=l1(i-6);
    l2(i)=l2(i-6);
    temp_x=eval(solve(l1(i)==l2(i),x));
    temp_y=subs(l1(i),x,temp_x);
    p3(i,:)=[temp_x temp_y,sqrt(R*R-temp_x*temp_x-temp_y*temp_y)];
end
%Area 11 to 12
for i=11:12
    alfa1=(-180-AZBLK(1))*pi/180+(i-7)*rot_angle;
    alfa2=(AZBLK(1))*pi/180+rot_angle*(i-11);
    p1(i,:)=[cos(alfa1) sin(alfa1) 0];
    p2(i,:)=[cos(alfa2) sin(alfa2) 0];
    l1(i)=l1(i-6);
    l2(i)=l2(i-6);
    temp_x=eval(solve(l1(i)==l2(i),x));
    temp_y=subs(l1(i),x,temp_x);
    p3(i,:)=[temp_x temp_y,sqrt(R*R-temp_x*temp_x-temp_y*temp_y)];
end

points_data_with_overlap=points_data;
points_data_without_overlap=points_data;
z_max=max(p3(:,3));
z_mid=min(p3(:,3));

for i=1:length(points_data)
   for j=1:length(points_data{1})
       x_point=points_data{i}(j,1);
       y_point=points_data{i}(j,2);
       z_point=points_data{i}(j,3);
       if abs(z_point)>z_max
           %overlapped points
           points_data_without_overlap{i}(j,:)=[0 0 0];
       elseif abs(z_point)>z_mid
           %Must be Area 1:4
           
       else
           %Must be Area5 or Area6
           if abs(x_point)<abs(p3(5,1))
               points_data_without_overlap{i}(j,:)=[0 0 0];
           elseif abs(x_point)<abs(p3(6,1))
               %Area 5  
               y1_temp=eval(subs(l1(5),x,x_point));
               y2_temp=eval(subs(l2(5),x,x_point));
               %need to correct. There are 12 areas, not 6
           else
               %Area 6 o 5
           end
       end
           
   end
end
