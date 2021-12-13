clear;close all;clc
%characteristic of LIDAR
    %Angle elevation range; [-16.677; 16,545]
        %In this particular device we have 16 laser beams/points
        %Also, every azimuth block has an angle offset (-1.24°)
        %AZBLK=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062]
    %About encode range: [0,90111]
    %In mode 1023 azimuth. Every azimuth increments in 88 ticks
        %1 tick = 360/90112= 0.003995°.     Aprox -> 0.004°
        %88 tick = 360*88/90112 = 0.3516°
    %In total, we obtain 1024*16 = 16384 points in every LIDAR scan
      
%Set Local variables:
AZBLK=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062];
azimuth_offset=-1.24*(pi/180); %without considering offset of every beam
n_AZBLK=1024;
ticks_per_azimuth=90112/n_AZBLK;
ang_diff=2*pi/n_AZBLK;
n_rotation=6;


%Interception with synthetic spheric room
%pre-allocation
temp=zeros(n_AZBLK*length(AZBLK),3);
points_data=cell(1,n_rotation);
for i=1:n_rotation
    points_data{i}=temp;
end

%Rotation -> 31° 
%Tiene que ser menor al la suma de angulos de elevacion, pero mayor a 30
%para cumplir con el barrido total
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
    points_data{i}=temp;    
end
figure (1)
hold
for i=1:n_rotation
    temp=points_data{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

%%
p1=[0,0,0];
p2=p1;
l1_param=[0,0,0];
l2_param=l1_param;
alfa1=AZBLK(1)*pi/180;
alfa2=AZBLK(end)*pi/180;
    
%Let's set the boundary points of the donuts and every pair of equation lines
for i=1:n_rotation
    p1(i,:)=[cos(rot_angle*(i-1)+alfa1),sin(rot_angle*(i-1)+alfa1),0];
    p2(i,:)=[cos(rot_angle*(i-1)+alfa2),sin(rot_angle*(i-1)+alfa2),0];
    l1_param(i,:)=[tan(rot_angle*(i-1)),p1(i,1),p1(i,2)];
    l2_param(i,:)=[tan(rot_angle*(i-1)),p2(i,1),p2(i,2)];
end

%Verify relevant data
new_points_data=points_data;
%Donut2
i=2;
for j=1:length(new_points_data{i})
    x_data=new_points_data{i}(j,1);
    y_data=new_points_data{i}(j,2);
    z_data=new_points_data{i}(j,3);
    if y_data<p1(1,2) && y_data>p2(1,2)
        new_points_data{i}(j,:)=[0,0,0];
    end
end
%Donut 3 and 4
for i=3:4
    for j=1:length(new_points_data{i})
        x_data=new_points_data{i}(j,1);
        y_data=new_points_data{i}(j,2);
        z_data=new_points_data{i}(j,3);
        if y_data>p1(1,2) 
            y_temp=line_equation(x_data,l1_param(i-1,:));
            if y_data<y_temp
                new_points_data{i}(j,:)=[0,0,0];
            end
        elseif y_data<p2(1,2)
            y_temp=line_equation(x_data,l2_param(i-1,:));
            if y_data>y_temp
                new_points_data{i}(j,:)=[0,0,0];
            end
        else
            new_points_data{i}(j,:)=[0,0,0];
        end
    end
end
%Donut 5 and 6
for i=5:6
    for j=1:length(new_points_data{i})
        x_data=new_points_data{i}(j,1);
        y_data=new_points_data{i}(j,2);
        z_data=new_points_data{i}(j,3);
        if y_data>p1(1,2) 
            y_temp=line_equation(x_data,l1_param(i-1,:));
            if y_data>y_temp%this part is different
                new_points_data{i}(j,:)=[0,0,0];
            end
        elseif y_data<p2(1,2)
            y_temp=line_equation(x_data,l2_param(i-1,:));
            if y_data<y_temp%this part is different
                new_points_data{i}(j,:)=[0,0,0];
            end
        else
            new_points_data{i}(j,:)=[0,0,0];
        end
    end
end

figure (2)
hold
%temp=[];
for i=1:n_rotation
    temp=[new_points_data{i}];
    %Para generar data a Carlos->temp=[temp;new_points_data{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

%% Analysis data
temp=zeros(length(points_data{1}),2);
figure (3)
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
%% Triangulate
figure(4);
hold
%k_puntos=16384;
%T=[];
for i=1:n_rotation
    temp=new_points_data{i};
    T=[triangulate_v1(temp,length(AZBLK),n_AZBLK)];
    %T=[T;triangulate_v1(temp,length(AZBLK),n_AZBLK)+(i-1)*k_puntos];
    trimesh(T,temp(:,1),temp(:,2),temp(:,3))
end
