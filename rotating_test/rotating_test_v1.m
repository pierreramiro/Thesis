clear;clc;close all
n_fig=0;
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
    angle_between_azimuths=2*pi*(ticks_between_azimuths/total_ticks);
                         %=2*pi/n_AZBLK
%Angle rotation,n° donuts,angular velocity,time per scan.
    time_per_scan=0.100;%100 ms
    n_donuts=6;
    rot_angle=-33.53706667*pi/180;
    w_lidar=rot_angle/time_per_scan; %%%%%%%%<---------set Angular velocity
    rot_matrix=[cos(rot_angle) -sin(rot_angle)  0;
                sin(rot_angle) cos(rot_angle)   0;
                0               0               1];
    delta_t_ray2ray=time_per_scan/n_points;
    alfa=delta_t_ray2ray * w_lidar;
    rot_due2w_r2r=[cos(alfa) -sin(alfa)  0;
                    sin(alfa) cos(alfa)   0;
                    0         0          1];
    delta_t_ray2AZBLK=delta_t_ray2ray;
    alfa=delta_t_ray2AZBLK * w_lidar;
    rot_due2w_r2A=[cos(alfa) -sin(alfa)  0;
                    sin(alfa) cos(alfa)   0;
                    0         0          1];
%pre-allocate point clouds
    temp=zeros(n_AZBLK*n_rays,3);
    point_cloud=cell(1,n_donuts);
    for i=1:n_donuts
        point_cloud{i}=temp;
    end
%% Generate synthetic sphere
R=1;%radius sphere
start_offset=0;
temp=zeros(n_rays,3);
count=0;
rot_due2w_acum=eye(3);
for i=1:n_donuts
    for j=1:n_AZBLK
        for k=1:n_rays-1
            x=R*-sin(AZBLK_angle(k));           
            y=(R*cos(AZBLK_angle(k))) * -sin(AZBLK_offset(k)-(j-1)*angle_between_azimuths);
            z=(R*cos(AZBLK_angle(k))) * -cos(AZBLK_offset(k)-(j-1)*angle_between_azimuths);
                                 %El signo menos de arriba _^_ es relevante.
                                 %Se debe a que el scan es sentido horario,
                                 %y se habia definido esta variable con
                                 %signo positivo.
            %let's consider the rotation due to the motor's "w"
            temp(k,:)=(rot_due2w_acum*[x;y;z])';               
            rot_due2w_acum=rot_due2w_r2r * rot_due2w_acum; 
            count=count+1;    
        end
        k=k+1;
        x=R*-sin(AZBLK_angle(k));           
        y=(R*cos(AZBLK_angle(k))) * -sin(AZBLK_offset(k)-(j-1)*angle_between_azimuths);
        z=(R*cos(AZBLK_angle(k))) * -cos(AZBLK_offset(k)-(j-1)*angle_between_azimuths);
        temp(k,:)=(rot_due2w_acum*[x;y;z])';
        %save the 16 rays of the j-th azimuth
        point_cloud{i}((j-1)*k+1:j*k,:)=temp;
        rot_due2w_acum=rot_due2w_r2A * rot_due2w_acum;
        count=count+1;
    end
    if w_lidar==0
        point_cloud{i}=(rot_matrix^(i-1)*point_cloud{i}')';
    end
end
n_fig=n_fig+1;
figure(n_fig)
hold
for i=1:n_donuts
    temp=point_cloud{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
title('Esfera completa')
%%
n_fig=n_fig+1;
figure(n_fig)
temp=point_cloud{1};
scatter3(temp(:,1),temp(:,2),temp(:,3))
title('Dona 1')
%%
n_fig=n_fig+1;
figure(n_fig)
hold
for i=1:2
    temp=point_cloud{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
title('Dona 1 y Dona 2')
%%
n_fig=n_fig+1;
figure(n_fig)
hold
temp=point_cloud{1};
scatter3(temp(:,1),temp(:,2),temp(:,3))
temp=point_cloud{6};
scatter3(temp(:,1),temp(:,2),temp(:,3),'c')
title('Dona 1 y Dona 6')
%% Set parameter of curves
i=1;
elev_angle=AZBLK_angle(1);
c1=sin(elev_angle);
c2=cos(elev_angle);

syms alfa1 alfa2 phi1
eq_1=-c1*cos(alfa1)+c2*sin(phi1)*sin(alfa1)==-c1*cos(rot_angle-alfa1)-c2*sin(phi1)*sin(rot_angle-alfa1);
eq_2=c1*sin(alfa1)+c2*sin(phi1)*cos(alfa1)==c1*sin(rot_angle-alfa1)-c2*sin(phi1)*cos(rot_angle-alfa1);

temp=vpasolve([eq_1;eq_2]);
alfa=eval(temp.alfa1);
phi=eval(temp.phi1);

sen_phi=sin(phi);
sen_alfa=sin(alfa);
cos_alfa=cos(alfa);
cos_phi=cos(phi);
c3=cos(rot_angle*(i-1));
c4=sin(rot_angle*(i-1));

x=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c3 + (c1*sen_alfa+c2*sen_phi*cos_alfa)*c4 );
y=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c4 - (c1*sen_alfa+c2*sen_phi*cos_alfa)*c3 );    
z=-R*c2*cos_phi;
cross_point(i,1:3)=[x,y,z];

%queda tarea hallar ymax y ymin

%% 
n_fig=n_fig+1;
figure(n_fig)
R=1;
rot_angle=-0.5853;
phi=-2*pi/0.1*(0:0.00001:0.9);
elev_angle=AZBLK_angle(end);
alfa=-rot_angle/2/pi*phi;
i=1;
sen_phi=sin(phi);
sen_alfa=sin(alfa);
cos_alfa=cos(alfa);
c1=sin(elev_angle);
c2=cos(elev_angle);
c3=cos(rot_angle*(i-1));
c4=sin(rot_angle*(i-1));

x=R*( (-c1*cos_alfa+c2*sen_phi.*sen_alfa)*c3 + (c1*sen_alfa+c2*sen_phi.*cos_alfa)*c4 );
y=R*( (-c1*cos_alfa+c2*sen_phi.*sen_alfa)*c4 - (c1*sen_alfa+c2*sen_phi.*cos_alfa)*c3 );    
z=-R*c2.*cos(phi);
plot3(x,y,z)

% elev_angle=AZBLK_angle(1);
% c1=sin(elev_angle);
% c2=cos(elev_angle);
% x=R*( (-c1*cos_alfa+c2*sen_phi.*sen_alfa)*c3 + (c1*sen_alfa+c2*sen_phi.*cos_alfa)*c4 );
% y=R*( (-c1*cos_alfa+c2*sen_phi.*sen_alfa)*c4 - (c1*sen_alfa+c2*sen_phi.*cos_alfa)*c3 );    
% z=-R*c2.*cos(phi);
% hold
% plot3(x,y,z)
% scatter3(cross_point(1),cross_point(2),cross_point(3))
% n_fig=n_fig+1;
% figure(n_fig)
% plot(x)
% n_fig=n_fig+1;
% figure(n_fig)
% plot(y)
% n_fig=n_fig+1;
% figure(n_fig)
% plot(z)
%% Proceed to supress redundant data
new_point_cloud=point_cloud;
%Donut 2
i=2;
count=0;
elev_angle=AZBLK_angle(1);
c1=sin(elev_angle);
c2=cos(elev_angle);
c3=cos(rot_angle*(i-2));
c4=sin(rot_angle*(i-2));
x_start=point_cloud{i}(16,1);
x_cross=cross_point(1,1);
for j=1:length(new_point_cloud{i})
    x_data=new_point_cloud{i}(j,1);
    y_data=new_point_cloud{i}(j,2);
    z_data=new_point_cloud{i}(j,3);
    if x_data>0
       if z_data>0
       
       else
           if y_data>0
              if x_data<x_start
                  new_point_cloud{i}(j,:)=[0,0,0];
              else
                  [x_temp,y_temp]=xy_curve_eq_v2(AZBLK_angle(end),i-1,rot_angle,z_data);
                  if x_data<x_temp
                      new_point_cloud{i}(j,:)=[0,0,0];              
                  end
              end
           end
       end
    else 
        if z_data>0
        
        else
            if x_data>x_cross
                new_point_cloud{i}(j,:)=[0,0,0];    
            else
                [x_temp,y_temp]=xy_curve_eq_v3(AZBLK_angle(1),i-1,rot_angle,z_data);
                if x_data>x_temp
                    new_point_cloud{i}(j,:)=[0,0,0];    
                end
            end
        end
        
    end
        
    
%     if inside_curve_eq_v1(rot_angle,x_data,y_data,z_data,cross_point,c1,c2,c3,c4);
%        new_point_cloud{i}(j,:)=[0,0,0];
%        count=count+1;
%     end
end
n_fig=n_fig+1;
figure(n_fig)
hold
for i=1:2
    temp=new_point_cloud{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
%% Triangle Mesh Generation
%%
% n_fig=n_fig+1;
% figure(n_fig)
% zp=point_cloud{2}(:,3);
% yp=point_cloud{2}(:,2);
% count_x=0;
% count_y=0;
% for i=1:length(zp)
% [x(i),y(i)]=xy_curve_eq(AZBLK_angle(end),1,rot_angle,yp(i),zp(i));
% if imag(x(i))
%     count_x=count_x+1;
% end
% if imag(y(i))
%     count_y=count_y+1;
% end
% end
% scatter3(x,y,zp)
% n_fig=n_fig+1;
% figure(n_fig)
% plot(x')
% n_fig=n_fig+1;
% figure(n_fig)
% plot(y')
% n_fig=n_fig+1;
% figure(n_fig)
% plot(zp)