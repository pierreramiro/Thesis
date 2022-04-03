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
for i=1:n_donuts
    temp=new_point_cloud{i};
    T=triangulate_v1(temp,n_rays,n_AZBLK);
    %trimesh(T,temp(:,1),temp(:,2),temp(:,3))
    %Data a Carlos
    Triangle_mesh=[Triangle_mesh;T+(i-1)*n_points];    
end
%% Hole 1. Tri-Mesh
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
T_tripivot=[];
k_ray=16;
k_offset_azimuth=floor(AZBLK_offset(end)/angle_between_azimuths);
minimal_offset=mod(AZBLK_offset(end),angle_between_azimuths);
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
    [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
    v3=(i-1)*n_points+j;
    %añadimos uno de los 16 triangulos a la malla triangular
    T_tripivot(17-p,:)=[v1 v2 v3];
end
%Procedemos a unir los puntos libres que se encuentran entre los triangulos
%Denominaremos dos tipos de puntos libres, aquellos que pertenecen a la Donut 1
%como free_points_L1 y aquellos de la Donut 2 como free_points_L2
%Dado que son 16 triangulos, habrán 15 zonas donde realizaremos la unión

freepoints=zeros(15,2);   %este vector nos muestra cuantos puntos libres
                            %habrán en la L1 (columna 1) y L2 (columna 2),
                            %de esta forma, verificamos que los puntos
                            %libres de L2 sean igual o mayor!
T=[];
for j=1:n_rays-1
    free_points_L1=(T_tripivot(j+1,3)+1-T_tripivot(j,3))/16;
    v1=T_tripivot(j+1,1);
    v2=T_tripivot(j,2);
    if v2>v1
        %debemos concatenar la dona, es decir, unir la parte final del
        %azimuth 1024 con el primer azimuth!
        free_points_L2=(v1-(v2-n_points))/16;
    else
        free_points_L2=(v1-v2)/16;
    end
    aristas_mismo_vertice=free_points_L2-free_points_L1;
    %creamos el primer triangulo con mismo vertice
    v1=T_tripivot(j+1,1) ;
    v2=T_tripivot(j+1,3);
    v3=T_tripivot(j+1,3)+1; % el v3 es el vertice en comun
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
v2=T_tripivot(16,2);
v3=T_tripivot(16,3);
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
tripivot_offset=16;
for p=1:16
    x=point_cloud_without_overlapped(v3,1);
    y=point_cloud_without_overlapped(v3,2);
    z=point_cloud_without_overlapped(v3,3);
    [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
    %añadimos uno de los 16 triangulos a la malla triangular
    T_tripivot(p+tripivot_offset,:)=[v1 v2 v3];
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

T=[];
for j=1:n_rays-1
    free_points_L1=(T_tripivot(j+1+tripivot_offset,3)-(T_tripivot(j+tripivot_offset,3)+17))/16; % el valor de 17 se 
    %debe a, que se tiene que aumentar 1 ray y 1 azimuth, dando asi la suma
    %de (16+1)=17
    v1=T_tripivot(j+1+tripivot_offset,1);
    v2=T_tripivot(j+tripivot_offset,2);
    free_points_L2=(v1-v2)/16;
    freepoints(j,:)=[free_points_L1,free_points_L2];
    aristas_mismo_vertice=free_points_L2-free_points_L1;
    %creamos el primer triangulo con mismo vertice, este triangulo se crea
    %independientemente existan aristas repetidas
    v1=T_tripivot(j+tripivot_offset,3) ;
    v2=T_tripivot(j+tripivot_offset,2);
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
Triangle_mesh=[Triangle_mesh;T];
%% Donut 3
i=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                  ZONA INFERIOR                %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=[];
k_ray=16;
%construimos lso triángulos pivot de la zona inferior
for p=1:16
    j=n_rays+1-p;
    while isequal(new_point_cloud{i}(j,:),[0 0 0])
        j=j+n_rays;
    end
    temp=rot_matrix^-(i-2)*new_point_cloud{i}(j,:)';
    x=temp(1);
    y=temp(2);
    z=temp(3);
    [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
    if isequal(new_point_cloud{i-1}(v1,:),[0 0 0])
        temp=new_point_cloud{i}(j,:);
        x=temp(1);
        y=temp(2);
        z=temp(3);
        [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
    else
        v1=(i-2)*n_points+v1;
        v2=(i-2)*n_points+v2;
    end
    v3=(i-1)*n_points+j;
    T=[T;v1,v2,v3];
end
T_tripivot=[T_tripivot;T];

%procedemos a llenar la zona inferior
    %freepoints=zeros(15,2);
    %Se verificó que hay mas ptos libres en L2 que L1
%IDEA: añadir un parametro a T, que indique en que dona está
special_zone=0;
for j=1:15
    free_points_L1=(T(j+1,3)+1-T(j,3))/16;
    v1=T(j+1,1);
    v2=T(j,2);
    if (v1-n_points*(i-2))*(v2-n_points*(i-2))<0
        %analizamos si la zona de llenado contiene triángulso pivots de
        %distintas Donuts
        special_zone=special_zone+1;
        continue
    else
        if v2>v1
        %Caso particular de la zona inferior que debemos concatenar la Donut
            v2=v2-n_points;
        end
        free_points_L2=(v1-v2)/16;
    end
        %freepoints(j,:)=[free_points_L1,free_points_L2];
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
    %freepoints
Triangle_mesh=[Triangle_mesh;T(17:end,:)]; %Los primeros 16 valores de T 
                                           %contienen a los vertices de los
                                           %triangulos pivot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                    ZONA MEDIA                 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ahora, realizamos la unión entre la zona media de la Donut 2 y la Donut 1
%para ello, usamos el vértice 2 y 3 del triangulo 16
v2=T_tripivot(end,2);
v3=T_tripivot(end,3);
T=[];
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
j=Triangle_mesh(end,3);
T=[];
for p=1:n_rays
    %obtenmos el punto extremo que será el v3 del triangulo pivot y
    %referenciamos a un nuevo sistema de referencia
    temp=rot_matrix^-(i-2)*point_cloud_without_overlapped(j,:)';
    x=temp(1);
    y=temp(2);
    z=temp(3);
    %hallamos el vertice según el azimuth
    [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
    if isequal(new_point_cloud{i-1}(v1,:),[0 0 0])
        %en caso la intersección no es con la Donut_k-1, entonces
        %analizamos la intersección con la Donut_ref
        temp=point_cloud_without_overlapped(j,:);
        %volvemos a obtener el punto sin giro del eje z
        x=temp(1);
        y=temp(2);
        z=temp(3);
        [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
    else
        v1=(i-2)*n_points+v1;
        v2=(i-2)*n_points+v2;
    end
    v3=j;
    T=[T;v1 v2 v3];
    if p==n_rays
        %no hago el cálculo del nuevo v3 para la ultima iteración
        continue
    else
        j=j+1;
        if isequal(point_cloud_without_overlapped(j+n_rays,:),[0 0 0])
            %Verificamos si el siguiente v3 del triangulo pivot se encuentra en
            %un azimut menor al actual
            while(isequal(point_cloud_without_overlapped(j,:),[0 0 0]))
            j=j-n_rays;
            end
        else
            while(~isequal(point_cloud_without_overlapped(j,:),[0 0 0]))
                j=j+n_rays;
            end
                j=j-n_rays;
        end
    end
end

%procedemos a llenar la zona inferior
    %freepoints=zeros(15,2);
    %Se verificó que hay mas ptos libres en L2 que L1
%IDEA: añadir un parametro a T, que indique en que dona está
special_zone=0;
for j=1:15
    free_points_L1=(T(j+1,3)-(T(j,3)+17))/16; % el valor de 17 se 
    %debe a, que se tiene que aumentar 1 ray y 1 azimuth, dando asi la suma
    %de (16+1)=17
    v1=T(j+1,1);
    v2=T(j,2);
    if (v1-n_points*(i-2))*(v2-n_points*(i-2))<0
        special_zone=special_zone+1;
        continue
    else
%         if v2>v1 
%             %esto es para el caso de donas concatenadas
%             v2=v2-n_points;
%         end
        free_points_L2=(v1-v2)/16;
    end
        %freepoints(j,:)=[free_points_L1,free_points_L2];
    aristas_mismo_vertice=free_points_L2-free_points_L1;
    %creamos el primer triangulo con mismo vertice, este triangulo se crea
    %independientemente existan aristas repetidas
    v1=T(j,3);
    v2=T(j,2);
    v3=v1+n_rays+1; % el v3 es el vertice en comun
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
%freepoints
T_tripivot=[T_tripivot;T(1:16,:)];
Triangle_mesh=[Triangle_mesh;T(17:end,:)];
%% Donut 4
i=4;
%Deseo conocer, de esta Donut, si existen filas completamente eliminadas y
%los puntos en donde están los limites de cada fila
boundary_points=get_boundary_points(new_point_cloud{i},n_rays,n_AZBLK);
T=[];
for j=2:3
    %Tenemos dos slot de triangulos, los superiores e inferiores. Aunque
    %tambien hay los del lado opuesto (Son, en realidad, 4 slots)
    for p=1:n_rays
        %obtenemos el v3 que nos ayudará a obtenes los v1 y v2 de c/triangulo
        v3=boundary_points{p}(j,1)+(i-1)*n_points;
        %hallamos las coordenadas del vertice
        point=boundary_points{p}(j,2:4);
        %le realizamos la rotacion para usar la función get_azimuth_vex()
        temp=rot_matrix^-(i-2)*point';
        x=temp(1);
        y=temp(2);
        z=temp(3);
        [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
        %debemos verificar y los vertices v1 y v2 obtenidos de la Donut
        %anterior son distintos de cero, caso contrario debemos volver a
        %hallar los vertices v1 y v2 pero con la Donut referencial
        if isequal(new_point_cloud{i-1}(v1,:),[0 0 0]) || isequal(new_point_cloud{i-1}(v2,:),[0 0 0])
            %hallamos los nuevos v1 y v2 de la Donut referencial
            temp=point;
            x=temp(1);
            y=temp(2);
            z=temp(3);
            [v1,v2]=get_azimuth_vex(x,y,z,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset);
        else
            %realizamos el offset a los vertices obtenidos
            v1=(i-2)*n_points+v1;
            v2=(i-2)*n_points+v2;
        end
        %guardamos el triangulo
        T=[T;v1,v2,v3];
    end
end
for slot=1:2
    for j=1:n_rays-1
    %establecemos variables para identificar los vertices de los dos
    %triangulos que contienen el espacio a llenar
    v3_act=T(j+n_rays*(slot-1),3);
    v2_act=T(j+n_rays*(slot-1),2);
    v1_act=T(j+n_rays*(slot-1),1);
    v3_next=T(j+1+n_rays*(slot-1),3);
    v2_next=T(j+1+n_rays*(slot-1),2);
    v1_next=T(j+1+n_rays*(slot-1),1);
    %definimos los parámetros según el modo y tipo
    if isequal(point_cloud_without_overlapped(v3_act-n_rays,:),[0 0 0])
        %Escalera_tipo_1 o Rampa_tipo_2
        v_init=v2_next;
        v_fin=v1_act;
        if v3_act>=v3_next
            %MODO ESCALERA TIPO 1
            %Obtenemos los puntos libres de L1
            common_vex=v3_act+1;
            free_points_L1=(common_vex-v3_next)/16;
            %creamos el triangulo con mismo vertice
            v1=v_fin;
            v2=v3_act;
            v3=common_vex;
            %esta variable servirá para cambiar el sentido del llenado
            %ya que el v_init sería el v_fin, coon respecto al common_vex
            escalera=boolean(1);
            tipo=1;
        else
            %MODO RAMPA TIPO 2
            %Obtenemos los puntos libres de L1
            common_vex=v3_next-17;
            free_points_L1=(common_vex-v3_act)/16;
            %creamos el triangulo con mismo vertice
            v1=v3_next;
            v2=v_init;
            v3=common_vex;
            escalera=boolean(0);
            tipo=2;
        end
    else
        %Escalera_tipo_2 o Rampa_tipo_1
        v_init=v2_act;
        v_fin=v1_next;
        if v3_act>=v3_next
            %MODO ESCALERA TIPO 2
            %Obtenemos los puntos libres de L1
            common_vex=v3_next-1;
            free_points_L1=(v3_act-common_vex)/16;
            %creamos el triangulo con mismo vertice
            v1=v_fin;
            v2=v3_next;
            v3=common_vex;
            escalera=boolean(1);
            tipo=2;
        else
            %MODO RAMPA TIPO 1
            %Obtenemos los puntos libres de L1
            common_vex=v3_act+17;
            free_points_L1=(v3_next-common_vex)/16;
            %creamos el triangulo con mismo vertice
            v1=v3_act;
            v2=v_init;
            v3=common_vex;
            escalera=boolean(0);
            tipo=1;
        end
    end
    %La variable escalera será usada como (-1)escalera
    %Luego de tener el modo y tipo, debemos ver si los vertices 
    %correspondientes a sus triangulos están en distintas Donuts, sino, se 
    %debe hacer un distinto llenado
    if (v_init-n_points)*(v_fin-n_points)>=0
        if v_init>v_fin
            free_points_L2=(v_fin-(v_init-n_points))/16;
        else
            free_points_L2=(v_fin-v_init)/16;
        end
        %Tambien debemos verificar que free_point_L2 sea mayor que los de L1
        if free_points_L2 < free_points_L1
            error("Los puntos en L2 son menos que los de L1, hacer otro llenado")
        end
        %Ahora hacemos el llenado que es el que más se repite
        %El triangulo con el common_vex ya fue creado al definir el modo y tipo
        T=[T;v1,v2,v3];
        aristas_mismo_vertice=free_points_L2-free_points_L1;
        if escalera
            %si es escalera el llenado es restando
            %A continuación, el llenado de los triangulos con misma arista
            for k=1:aristas_mismo_vertice
                v2=v1;
                v1=v1-n_rays;
                T=[T;v1,v2,v3];
            end
            %Ahora, el llenado que une los puntos libres de L1 y L2
            for k=1:free_points_L1
                %Primer triangulo
                v2=v1;
                v1=v1-n_rays;
                T=[T;v1,v2,v3];
                %Segundo triangulo y dependiendo del tipo, se suma o resta
                v2=v3;
                if tipo==1
                    v3=v3-n_rays;
                else
                    v3=v3+n_rays;
                end
                T=[T;v1,v2,v3];
            end
        else
            %El llenado del modo RAMPA
            %A continuación, el llenado con el common_vex
            for k=1:aristas_mismo_vertice
                v1=v2;
                v2=v1+n_rays;
                T=[T;v1,v2,v3];
            end
            %Procedemos, al llenado de los triángulo rectángulos
            for k=1:free_points_L1
                %Primer triangulo
                v1=v2;
                if tipo==1
                    v2=v3+n_rays;
                else
                    v2=v3-n_rays;
                end
                T=[T;v1,v2,v3];
                %Segundo triangulo
                v3=v2;
                v2=v1+n_rays;
                T=[T;v1,v2,v3];
            end
        end
    else
        %llenado especial, triangulos de distintas Donuts
        free_points_L2=NaN;
    end
    %free_points_vector(j,:)=[free_points_L1 free_points_L2];
    end
end
%T_tripivot=[T_tripivot;T];
%%
%Verificamos:
figure (4)
temp=point_cloud_without_overlapped;   
trimesh([Triangle_mesh;T_tripivot;T],temp(:,1),temp(:,2),temp(:,3))
title('primer intento de union')
