clear;clc;close all
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
%pre-allocate point clouds
    temp=zeros(n_AZBLK*n_rays,3);
    point_cloud=cell(1,n_donuts);
    for i=1:n_donuts
        point_cloud{i}=temp;
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
    point_cloud{i}=temp;    
end
figure (1)
hold
for i=1:n_donuts
    temp=point_cloud{i};
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

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
new_point_cloud=point_cloud;

%Donut2
i=2;
for j=1:length(new_point_cloud{i})
    x_data=new_point_cloud{i}(j,1);
    y_data=new_point_cloud{i}(j,2);
    z_data=new_point_cloud{i}(j,3);
    if x_data<=p1(1,1) && x_data>=p2(1,1)
        new_point_cloud{i}(j,:)=[0,0,0];
    end
end
%Donut 3 to 6
for i=3:n_donuts
    for j=1:length(new_point_cloud{i})
        x_data=new_point_cloud{i}(j,1);
        y_data=new_point_cloud{i}(j,2);
        z_data=new_point_cloud{i}(j,3);
        if x_data>p1(1,1) 
            y_temp=line_equation(x_data,l1_param(i-1,:));
            if y_data>=y_temp
                new_point_cloud{i}(j,:)=[0,0,0];
            end
        elseif x_data<p2(1,1)
            y_temp=line_equation(x_data,l2_param(i-1,:));
            if y_data<=y_temp
                new_point_cloud{i}(j,:)=[0,0,0];
            end
        else
            new_point_cloud{i}(j,:)=[0,0,0];
        end
    end
end
%Data a Carlos:
point_cloud_without_overlapped=zeros(n_points*n_donuts,3);

figure(2)
hold
for i=1:n_donuts
    temp=[new_point_cloud{i}];
    scatter3(temp(:,1),temp(:,2),temp(:,3));
    %Data a Carlos
    point_cloud_without_overlapped((i-1)*n_points+1:i*n_points,:)=temp;
end
hold off

%% Triangulate

%Data a Carlos
Triangle_mesh=[];

figure(3);
hold
for i=1:n_donuts
    temp=new_point_cloud{i};
    T=triangulate_v1(temp,n_rays,n_AZBLK);
    %trimesh(T,temp(:,1),temp(:,2),temp(:,3))
    %Data a Carlos
    Triangle_mesh=[Triangle_mesh;T+(i-1)*n_points];    
end
hold off
close(3)
%% Triangle mesh of the holes
minimal_offset=mod(AZBLK_offset(end),angle_between_azimuths);
%Donut 2
i=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                  ZONA INFERIOR                %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp=new_point_cloud{i};%cargamos los datos de la Donut2 de la esfera sin traslape
%Creamos los 16 triangulos de unión
T=[];
for p=16:-1:1
    for j=p:16:n_points
        if ~isequal(temp(j,:),[0 0 0])
            break
        end
    end
    %Hallamos el angulo del azimuth con el punto hallado previamente
    x=temp(j,1);
    y=temp(j,2);
    z=temp(j,3);
    angulo_circulo=-acos(x/R);
    alfa=asin(-y/R/sin(angulo_circulo));
    decimales=mod(alfa,angle_between_azimuths);
    %Procedemos a buscar los Azimuths de la Donut 1 que contienen a los dos 
    %vertices que se conectarán con el vertice de la dona 2
    k_azimuth=floor(alfa/angle_between_azimuths)-1;
    if k_azimuth<=0
        k_azimuth=k_azimuth+n_AZBLK;
    end
    %arriba definimos el minimal_offset
    if decimales>minimal_offset %recordar que ambos valores son negativos
        %Caso 1
        if k_azimuth==1
            v1=(i-2)*n_points+n_points;
            v2=(i-2)*n_points+n_rays;
        else
            v1=(i-2)*n_points+(k_azimuth-1)*n_rays;
            v2=(i-2)*n_points+k_azimuth*n_rays;
        end
        %Probar el alfa con:
        %sin(-point_cloud{1}(v1-32,2)/sin(AZBLK_angle(end)-pi/2))
    else
        %Caso 2
        if k_azimuth==n_AZBLK
            v1=(i-2)*n_points+n_points;
            v2=(i-2)*n_points+n_rays;
        else
            v1=(i-2)*n_points+k_azimuth*n_rays;
            v2=(i-2)*n_points+(k_azimuth+1)*n_rays;
        end
    end
    v3=(i-1)*n_points+j;
    %añadimos uno de los 16 triangulos a la malla triangular
    T(17-p,:)=[v1 v2 v3];
end
%Procedemos a unir los puntos libres que se encuentran entre los triangulos
%Denominaremos dos tipos de puntos libres, aquellos que pertenecen a la Donut 1
%como free_points_L1 y aquellos de la Donut 2 como free_points_L2
%Dado que son 16 triangulos, habrán 15 zonas donde realizaremos la unión

freepoints=zeros(15,2);   %este vector nos muestra cuantos puntos libres
                            %habrán en la L1 (columna 1) y L2 (columna 2),
                            %de esta forma, verificamos que los puntos
                            %libres de L2 sean igual o mayor!
for j=1:n_rays-1
    free_points_L1=(T(j+1,3)+1-T(j,3))/16;
    v1=T(j+1,1);
    v2=T(j,2);
    if v2>v1
        %debemos concatenar la dona, es decir, unir la parte final del
        %azimuth 1024 con el primer azimuth!
        free_points_L2=(v1-(v2-n_points))/16;
    else
        free_points_L2=(v1-v2)/16;
    end
    aristas_mismo_vertice=free_points_L2-free_points_L1;
    %creamos el primer triangulo con mismo vertice
    v1=T(j+1,1) ;
    v2=T(j+1,3);
    v3=T(j+1,3)+1; % el v3 es el vertice en comun
    T=[T;v1,v2,v3];
    % ahora creamos los triangulos que comparten el mismo vertice v3
    for k=1:aristas_mismo_vertice
        v2=T(end,1);
        v1=v2-16;
        if v1<=0
            v1=v1+n_points*(i-1);
        end
        T=[T;v1 v2 v3];
    end
    %finalmente, nos quedan "zonas rectangulares" por completar
    for k=1:free_points_L1
        v1=T(end,1)-16;
        v2=T(end,1);
        v3=T(end,3);
        if v1<=0
            v1=v1+n_points*(i-1);
        end
        T=[T;v1 v2 v3];
        v2=v3;
        v3=v2-16;
        T=[T;v1 v2 v3];
    end
    freepoints(j,:)=[free_points_L1,free_points_L2];
end
% figure
% temp=point_cloud_without_overlapped;
% trimesh([Triangle_mesh;T],temp(:,1),temp(:,2),temp(:,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                    ZONA MEDIA                 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ahora, realizamos la unión entre la zona media de la Donut 2 y la Donut 1
%para ello, usamos el vértice 2 y 3 del triangulo 16
v2=T(16,2);
v3=T(16,3);
while(~isequal(point_cloud_without_overlapped(v3+16,:),[0 0 0]))
    v1=v2;
    v2=v3+16;
    T=[T;v1 v2 v3];
    v2=v1+16;
    v3=v3+16;
    T=[T;v1 v2 v3];
end
Triangle_mesh=[Triangle_mesh;T];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                  ZONA SUPERIOR                %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Ahora continuamos creando los triangulos pivtos de la zona superior, en
%este caso, ya tenemos unos de los puntos de estos triangulos: el v3. Para
%este punto realizamos un analisis particular para verificar si es o no el
%adecuado, ya que este triángulo viene de la unión previa. 
%   En suma, se observa por la gráfica que el triangulo no es el adecuado y es 
%necesario hallar el nuevo tirangulo pivot.
%   El cálculo a realizar es analogamente al realizado para los triangulos 
%de la zona inferior
T=[];
for p=1:16
    x=point_cloud_without_overlapped(v3,1);
    y=point_cloud_without_overlapped(v3,2);
    z=point_cloud_without_overlapped(v3,3);
    angulo_circulo=-acos(x/R);%angulo de giro en el eje x
    alfa=asin(-y/R/sin(angulo_circulo));
    %dado que el asin tiene un rango de valores [-pi/2,pi/2] debemos hacer
    %la correción del alfa verdadero, ya que conocemos que este debe ser
    %menor a -pi/2
    alfa=-alfa-pi;
    decimales=mod(alfa,angle_between_azimuths);
    k_azimuth=floor(alfa/angle_between_azimuths)-1;%notar que para este caso los
                                                   %valores de k_azimuth son
                                                   %positivos, es decir, que no
                                                   %es necesario concatenar la
                                                   %dona.
    if decimales>minimal_offset
        %Caso1    
        v1=(i-2)*n_points+(k_azimuth-1)*n_rays;
        v2=(i-2)*n_points+k_azimuth*n_rays;
    else
        %Caso 2
        v1=(i-2)*n_points+k_azimuth*n_rays;
        v2=(i-2)*n_points+(k_azimuth+1)*n_rays;
    end
    %añadimos uno de los 16 triangulos a la malla triangular
    T(p,:)=[v1 v2 v3];
    if p==16
        continue
    else
        for j=v3+1:16:n_points*i
            if isequal(point_cloud_without_overlapped(j,:),[0 0 0])
                v3=j-16;
                break
            end
        end
    end
end

%Procedemos a rellenar las zonas rectangulares
% freepoints=zeros(15,2);   %este vector nos muestra cuantos puntos libres
                            %habrán en la L1 (columna 1) y L2 (columna 2),
                            %de esta forma, verificamos que los puntos
                            %libres de L2 sean igual o mayor!
                            %De lo contrario, se debe analizar la nueva
                            %unión de puntos

for j=1:n_rays-1
    free_points_L1=(T(j+1,3)-(T(j,3)+17))/16; % el valor de 17 se 
    %debe a, que se tiene que aumentar 1 ray y 1 azimuth, dando asi la suma
    %de (16+1)=17
    v1=T(j+1,1);
    v2=T(j,2);
    free_points_L2=(v1-v2)/16;
    freepoints(j,:)=[free_points_L1,free_points_L2];
    aristas_mismo_vertice=free_points_L2-free_points_L1;
    %creamos el primer triangulo con mismo vertice, este triangulo se crea
    %independientemente existan aristas repetidas
    v1=T(j,3) ;
    v2=T(j,2);
    v3=v1+17; % el v3 es el vertice en comun
    T=[T;v1,v2,v3];
    % ahora creamos los triangulos que comparten el mismo vertice v3
    for k=1:aristas_mismo_vertice
        v1=T(end,2);
        v2=v2+16;
        T=[T;v1 v2 v3];
    end
    %finalmente, nos quedan "zonas rectangulares" por completar
    for k=1:free_points_L1
        v1=v2;
        v2=v3+16;
        T=[T;v1 v2 v3];
        v2=v1+16;
        v3=v3+16;
        T=[T;v1 v2 v3];
    end
end

%%
%Verificamos:
figure (4)
temp=point_cloud_without_overlapped;   
trimesh([Triangle_mesh;T],temp(:,1),temp(:,2),temp(:,3))
title('primer intento de union')
%%%%%%%%%%%%


%% Analyzing data (from step3)
% temp=zeros(n_points,2);
% for i=1:n_donuts
%     figure(3+i)
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