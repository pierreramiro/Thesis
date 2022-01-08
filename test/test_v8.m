clear; close all;clc
Total_time_start=tic;
%numero de donuts
%angle donuts
%rango elevación
%n_azimuths
%% Set LIDAR's parameters
%Angle elevation range;
    beam_altitude_angles=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062]*pi/180;
    %beam_altitude_angles=flip([-20,-17.3333,-14.6667,-12,-9.33333,-6.66667,-4,-1.33333,1.33333,4,6.66667,9.33333,12,14.6667,17.3333,20]*pi/180);
    n_beams=length(beam_altitude_angles);
%Offset angle of AZBLK's points
    %beam_azimuth_angles=regexprep(num2str(linspace(-1.24,-0.857,16)),'\s+',',')
    beam_azimuth_angles=[-1.24,-1.2145,-1.1889,-1.1634,-1.1379,-1.1123,-1.0868,-1.0613,-1.0357,-1.0102,-0.98467,-0.95913,-0.9336,-0.90807,-0.88253,-0.857]*pi/180;
%Encoder range: [0,90111] ticks
    total_ticks=90112;
%Mode 1023 azimuth. Every azimuth increments in 88 ticks
    %1 tick = 360/90112= 0.003995°.     Aprox -> 0.004°
    %88 tick = 360*88/90112 = 0.3516°
    n_AZBLK=1024;
    %It should be noted that, we obtain 1024*16 = 16384 points in every LIDAR scan
    n_points=n_AZBLK*n_beams;
    ticks_between_azimuths=total_ticks/n_AZBLK;
    angle_between_azimuths=-2*pi/n_AZBLK;
%Angle rotation,n° donuts,angular velocity,time per scan.
    rot_angle=-33.53706667*pi/180;  %Rango operable [-135° ~ -15°]
                                    %Ojo, notar que en -15 estamos por
                                    %debajo del beam_altitude_angle
    n_donuts=1;
    while (n_donuts*rot_angle>-pi-beam_altitude_angles(end)+beam_altitude_angles(end))
        n_donuts=n_donuts+1;
    end
    rot_matrix=[cos(rot_angle) -sin(rot_angle)  0;
                sin(rot_angle) cos(rot_angle)   0;
                0               0               1];
    rot_x_axis=@(alfa) [1 0         0;
                        0 cos(alfa) -sin(alfa);
                        0 sin(alfa) cos(alfa)];
    rot_y_axis=@(alfa) [cos(alfa)   0   sin(alfa);
                0           1   0;
                -sin(alfa)  0   cos(alfa)];
    rot_z_axis=@(alfa) [cos(alfa) -sin(alfa)  0;
                sin(alfa) cos(alfa)   0;
                0           0     1];
    
%pre-allocate point clouds
    Point_Cloud=zeros(n_AZBLK*n_beams*n_donuts,3);
%% Generate synthetic sphere
Time_sphere_generation_start=tic;
R=1;%RADIO
%Generate referencial azimuth--->XZ
for j=1:n_beams
    x=R*cos(beam_altitude_angles(j)-pi/2);
    y=0;
    z=R*sin(beam_altitude_angles(j)-pi/2);
    temp=rot_x_axis(beam_azimuth_angles(j))*[x,y,z]';
    Point_Cloud(j,1)=temp(1);
    Point_Cloud(j,2)=temp(2);
    Point_Cloud(j,3)=temp(3);
end
%Define rotation matrix for the azimuths
temp_rot_mat=rot_x_axis(angle_between_azimuths);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Referencial Donut%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate mirror azimuths for every azimuth calculated
for j=1:n_beams
    %Set acutal points
    temp_x=Point_Cloud(j,1);
    temp_y=Point_Cloud(j,2);
    temp_z=Point_Cloud(j,3);
    %mirror points from quarter Donut
    Point_Cloud(j+n_AZBLK/4*n_beams,1)=temp_x;
    Point_Cloud(j+n_AZBLK/4*n_beams,2)=temp_z;
    Point_Cloud(j+n_AZBLK/4*n_beams,3)=-temp_y;
    %mirror points from midle Donut
    Point_Cloud(j+n_AZBLK/2*n_beams,1)=temp_x;
    Point_Cloud(j+n_AZBLK/2*n_beams,2)=-temp_y;
    Point_Cloud(j+n_AZBLK/2*n_beams,3)=-temp_z;
    %mirror points from midle Donut
    Point_Cloud(j+n_AZBLK*3/4*n_beams,1)=temp_x;
    Point_Cloud(j+n_AZBLK*3/4*n_beams,2)=-temp_z;
    Point_Cloud(j+n_AZBLK*3/4*n_beams,3)=temp_y;
end
for i=2:n_AZBLK/4
    for j=1:n_beams
        %Calculate previous point
        x=Point_Cloud((i-2)*n_beams+j,1);
        y=Point_Cloud((i-2)*n_beams+j,2);
        z=Point_Cloud((i-2)*n_beams+j,3);
        %rotate that point
        temp=temp_rot_mat*[x,y,z]';
        %Set the new azimuth
        Point_Cloud((i-1)*n_beams+j,1)=temp(1);
        Point_Cloud((i-1)*n_beams+j,2)=temp(2);
        Point_Cloud((i-1)*n_beams+j,3)=temp(3);
        %mirror from quarter Donunt
        Point_Cloud((i-1+n_AZBLK/4)*n_beams+j,1)=temp(1);
        Point_Cloud((i-1+n_AZBLK/4)*n_beams+j,2)=temp(3);
        Point_Cloud((i-1+n_AZBLK/4)*n_beams+j,3)=-temp(2);
        %mirror points from midle Donut
        Point_Cloud((i-1+n_AZBLK/2)*n_beams+j,1)=temp(1);
        Point_Cloud((i-1+n_AZBLK/2)*n_beams+j,2)=-temp(2);
        Point_Cloud((i-1+n_AZBLK/2)*n_beams+j,3)=-temp(3);
        %mirror points from midle Donut
        Point_Cloud((i-1+n_AZBLK*3/4)*n_beams+j,1)=temp(1);
        Point_Cloud((i-1+n_AZBLK*3/4)*n_beams+j,2)=-temp(3);
        Point_Cloud((i-1+n_AZBLK*3/4)*n_beams+j,3)=temp(2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotate Referencial Donut%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:n_donuts
    temp=rot_matrix^(i-1)*Point_Cloud(1:n_points,:)';
    temp=temp';
    Point_Cloud(n_points*(i-1)+1:n_points*i,:)=temp;
end
Time_sphere_generation_end=toc(Time_sphere_generation_start);
fprintf("El tiempo de ejecución para generar la esfera sintética fue de: %.4f segundos\n",Time_sphere_generation_end);
figure (1)
hold on
for i=1:n_donuts
    temp=Point_Cloud(n_points*(i-1)+1:n_points*i,:);
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off
%% Supress redundant data
Time_fusion_point_cloud_start=tic;
%Set boundary vertical limits
X_L1=Point_Cloud(1,1);
X_Ln_beam=Point_Cloud(n_beams,1);
%Set eq line
eq_line=@(x,xb,yb,m) m*(x-xb)+yb;
%Calculate parameters for eq line
L=zeros(n_donuts,3);
for i=2:n_donuts
    x1=Point_Cloud((i-1)*n_points+1,1);%x1
    y1=Point_Cloud((i-1)*n_points+1,2);%y1
    x1_next=Point_Cloud((i-1)*n_points+n_beams+1,1);%x1_next
    y1_next=Point_Cloud((i-1)*n_points+n_beams+1,2);%y1_next
    m=(y1_next-y1)/(x1_next-x1);

    L(i,1)=x1;
    L(i,2)=y1;
    L(i,3)=Point_Cloud((i-1)*n_points+n_beams,1);%xn_beam
    L(i,4)=Point_Cloud((i-1)*n_points+n_beams,2);%yn_beam
    L(i,5)=m;
end
%Supress redundant data for Donut 2
i=2;
Total_points_deleted=0;
points_deleted=0;
for j=1:n_points 
    x=Point_Cloud((i-1)*n_points+j,1);
    %%y=Point_Cloud((i-1)*n_points+j,2);
    %%z=Point_Cloud((i-1)*n_points+j,3);
    if X_Ln_beam<=x && x<=X_L1
        Point_Cloud((i-1)*n_points+j,:)=[0 0 0];
        points_deleted=points_deleted+1;
    end
end
Total_points_deleted=Total_points_deleted+points_deleted;
fprintf("\tSe eliminaron %d puntos de la Donut %d (%.2f%% eliminados)\n",points_deleted,i,points_deleted*100/n_points);
%Supress redundant data for all Donuts
for i=3:n_donuts
   points_deleted=0;
   for j=1:n_points
       x=Point_Cloud((i-1)*n_points+j,1);
       y=Point_Cloud((i-1)*n_points+j,2);
       %%z=Point_Cloud((i-1)*n_points+j,3);
       if X_Ln_beam <= x
           if x<=X_L1
               Point_Cloud((i-1)*n_points+j,:)=[0 0 0];
               points_deleted=points_deleted+1;
               continue;
           else
               left_side=boolean(0);
           end
       else
           left_side=boolean(1);%Está a la zquierda
       end
       y_temp=eq_line(x,L(i-1,1+left_side*2),L(i-1,2+left_side*2),L(i-1,5));
       y_temp=y_temp*(-1)^left_side;
       y=y*(-1)^left_side;
       if y>=y_temp
           Point_Cloud((i-1)*n_points+j,:)=[0 0 0];
           points_deleted=points_deleted+1;
       end
   end 
   Total_points_deleted=Total_points_deleted+points_deleted;
   fprintf("\tSe eliminaron %d puntos de la Donut %d (%.2f%% eliminados)\n",points_deleted,i,points_deleted*100/n_points);
end
fprintf("\tEn total, se eliminaron %d de un total de %d puntos (%.2f%% eliminados)\n",Total_points_deleted,n_points*n_donuts,Total_points_deleted*100/n_donuts/n_points);
%verificamos si la última Donut fue elminada totalmente
if points_deleted==n_points
    n_donuts=n_donuts-1;
end
Time_fusion_point_cloud_end=toc(Time_fusion_point_cloud_start);
fprintf("El tiempo de ejecución para la fusión de la nube de puntos es de: %.4f segundos\n\n",Time_fusion_point_cloud_end);
figure (2)
hold on
for i=1:n_donuts
    temp=Point_Cloud((i-1)*n_points+1:n_points*i,:);
    scatter3(temp(:,1),temp(:,2),temp(:,3))
end
hold off

%% Triangulate One-Donut Fill
Time_OneDonutFill_start=tic;
%%%%%ANALIZAR esto para bits
mask=(n_points-1);
%0x03FFF;%bit masking of 32 bits
%mas puntos

T=[];
%Generate reference triangulation
for j=1:n_AZBLK
    for g=1:n_beams-1
        v1=(j-1)*n_beams+g;
        v3=v1+1;
        v2=bitand(v1+n_beams,mask)+1;
        T=[T;v1,v2,v3];
        v3=v2;
        v2=v3-1;
        T=[T;v1,v2,v3];
    end
end
OneDonutFill=T;
n_triangles_per_donut=length(OneDonutFill);
for i=2:n_donuts
    count=0;
    T=[];
    for j=1:n_triangles_per_donut
        point1=Point_Cloud(OneDonutFill(j,1)+(i-1)*n_points,:);
        if ~isequal(point1,[0 0 0])
            point2=Point_Cloud(OneDonutFill(j,2)+(i-1)*n_points,:);
            if ~isequal(point2,[0 0 0])
                point3=Point_Cloud(OneDonutFill(j,3)+(i-1)*n_points,:);
                if ~isequal(point3,[0 0 0])
                    count=count+1;
                    T(count,:)=OneDonutFill(j,:)+(i-1)*n_points;
                end
            end
        end
    end
    OneDonutFill=[OneDonutFill;T];
end
Time_OneDonutFill_end=toc(Time_OneDonutFill_start);
fprintf("El tiempo de ejecución para el One-Donut Fill es de: %.4f segundos\n",Time_OneDonutFill_end);
%% Triangulate Two-Donut Fill AND Tri-Donut Fill
Time_TDonutFill_start=tic;
TwoDonutFill=[];
TriDonutFill=[];
TriMiddleFill=[];
T=[];
T_Temp=[];
count_temp=0;
last_Tripivots=zeros(4,3);
for i=2:n_donuts
    %Define sector of fill
    Tripivot_middle=[];
    Tripivot_middle_particular=[];
    %last_Tripivots=[];
    for sector=1:4
        %Transform sector to bits
        bits_sector=de2bi(sector-1,2);
        Tripivot=[];
        n_tripivot=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                                 %%%
        %%%             TRIPIVOT            %%%
        %%%                                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:n_beams
            %Found boundary points
            %init index será (n_points-n_beams)+j o (n_points/2-n_beams)+j
            init_index=n_AZBLK/2*(2-xor(bits_sector(1),bits_sector(2)))-2;
            init_index=bitand(init_index*n_beams+j-1,mask)+1;
            %^^^Note that is missing the offset of points (i-1)*n_points^^^
            %paso: sera el diferencial el cual sumará o restará al init_index
            paso=n_beams*(-1)^bits_sector(1);
            v3=0;%considerar algun tipo de bandera para saber si se 
            %halló el v3 o no. Esto interesa para el caso de la Donut 6
            for count=1:n_AZBLK/4
                %Entramos en bucle hasta encontrar los puntos límetes de
                %cada sector. En el condicional agregamos el offset para cada
                %Donut correspondiente
                eval_point=Point_Cloud(init_index+(i-1)*n_points,:);
                if ~isequal(eval_point,[0 0 0])
                    %se encontró el v3
                    v3=init_index+(i-1)*n_points;
                    %Realizamos un conteo de los triángulos
                    n_tripivot=n_tripivot+1;
                    break
                end
                init_index=bitand(init_index+paso-1,mask)+1;
            end
            if v3==0
                %fprintf("no consegui para "+i+"en sector "+sector+" y j igual a "+j +"\n");
                continue
            end
%%%%--------%Hallamos el triángulo pivot
            point=Point_Cloud(v3,:);
            %Calculamos el azimuth que corresponde al alfa
            %segun el sector, trabajamos con un cierto beam_elevate_angle
            k_beam=(n_beams-1)*(1-bits_sector(2))+1;
            %realizamos la operación
            [v1,v2]=get_tripivot(point,R,rot_matrix,sector,n_beams,n_AZBLK,i,...
            k_beam,beam_azimuth_angles(k_beam),angle_between_azimuths,mask);
            if bits_sector(2)
                %definimos el sentido horario Ya que para los sectors 3 y 4
                %El sentido de los vértices es distinto a los de los
                %primeros sector. Entonces para seguir la jerarquía de los
                %sentidos, cambiamos aqui
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Esto podría ser o o importante                        %
                %%% Capaz, se puede definir un sentido para un lado y otro%
                %%% para los otros sector ()sector3 y sector4)                  %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                temp=v1;
                v1=v2;
                v2=temp;
            end
            %Tripivot(n_tripivot,:)=[v1 v2 v3];
            if isequal(Point_Cloud(v1,:),[0 0 0])
                v_temp=v2;
                point=Point_Cloud(v_temp,:);
                [v1,v2]=get_tripivot(point,R,rot_matrix,sector,n_beams,n_AZBLK,i-1,...
                k_beam,beam_azimuth_angles(k_beam),angle_between_azimuths,mask);
                if bits_sector(2)
                    %Deseo el nuevo v1
                    v1=v2;
                end
                v2=v_temp;
                %Tripivot(n_tripivot,:)=[v1 v2 v3];
            end
            if isequal(Point_Cloud(v2,:),[0 0 0])
                v_temp=v1;
                point=Point_Cloud(v_temp,:);
                [v1,v2]=get_tripivot(point,R,rot_matrix,sector,n_beams,n_AZBLK,i-1,...
                k_beam,beam_azimuth_angles(k_beam),angle_between_azimuths,mask);
                if bits_sector(2)
                    %Deseo el nuevo v2
                    v2=v1;
                end
                v1=v_temp;
                %Tripivot(n_tripivot,:)=[v1 v2 v3];
            end
            %debemos verificar que no exista concurrencia de Tripivots
            if not(xor(bits_sector(1),bits_sector(2)))
                %sector 1 y 4
                if v2==last_Tripivots(sector,2)
                   v2=last_Tripivots(sector,3);
                   %Tripivot(n_tripivot,:)=[v1 v2 v3]; 
                end
            else
                %sector 2 y 3
                if v1==last_Tripivots(sector,1)
                   v1=last_Tripivots(sector,3);
                   %Tripivot(n_tripivot,:)=[v1 v2 v3]; 
                end
            end
            Tripivot(n_tripivot,:)=[v1 v2 v3];
        end
%%%%%%Fin codígo tripivot
        TwoDonutFill=[TwoDonutFill;Tripivot];
        %Hallamos los triangulos que limitan la zona del medio
        triangle_sector=Tripivot(1+(n_tripivot-1)*bits_sector(2),:);
        Tripivot_middle=[Tripivot_middle;triangle_sector];
        if i==n_donuts &&length(Tripivot)>1 
            triangle_sector=Tripivot(1+(n_tripivot-1)*not(bits_sector(2)),:);
            Tripivot_middle_particular=[Tripivot_middle_particular;triangle_sector];
        end
        %--------------FIN TRIPIVOT--------------%
        %Con los triangulos pivot obtenidos, procedemos a realizar el
        %llenado de las zonas que encierran los triangulos. Hay dos tipos
        %de zonas: los sector y los del medio
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                                     %%%
        %%%     TRIDONUT FILL && TWODONUTFILL   %%%
        %%%                                     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:n_tripivot-1
            v3_act=Tripivot(j,3);
%             if v3_act==0
%                 continue
%             end%esto ya no va
            %verificiamos que modo y tipo de zona será
            %analizamos el valor del punto siguiente
            temp_point=Point_Cloud(v3_act+1+n_beams*bits_sector(1),:);
            if isequal(temp_point,[0 0 0])
                %No existe punto
                escalera=bits_sector(1);
            else
                %si existe punto
                escalera=not(bits_sector(1));
            end
            tipo= not(xor(bits_sector(1),escalera));
            %Luego de definir el modo y el tipo. Podemos definir el pasoL1
            %para llegar al v3_next
            pasoL1=n_beams*(-1)^xor(escalera,tipo);
            %Ahora debemos definir el vertice que estará en la esquina y
            %su lugar puede estar en el beam actual o en el beam siguiente
            v_corner=Tripivot(j+tipo,3)+(-1)^tipo+pasoL1*not(escalera);
            %v3_next=Tripivot(j+1,3);            
            %definimos los límites
            v_init=Tripivot(j+tipo,2-escalera);
            %v_init=Tripivot(j+tipo, 1+xor(bits_sector(1),tipo) );
            v_fin=Tripivot(j+not(tipo),2-not(escalera));
            %v_fin=Tripivot(j+not(tipo),2-xor(bits_sector(1),tipo) );
            v_corner_fin=Tripivot(j+not(tipo),3);
            pasoL2=(-1)^xor(bits_sector(2),tipo)*pasoL1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%-----DEPRECATED-----%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %------>>Para verificar si es de doble o tri llenado, se resta n_points
            %------>>ya que si fuera TriDonutFill, una de las Donuts invoolucradas
            %------>>seria la referencial. De esta forma obtenemos un signo negativo
            %------>>if (v_init-n_points)*(v_fin-n_points)>=0 
            %%%%%%%%%%%%%%%%%%%%^^^^^^^-----DEPRECATED-----^^^^%%%%%%%%%%%%%%%%%%%%%%%
            %Se cambió la condición anterior con esta nueva que evalua los bits MSB en vez de usar valores negativos
            condicionUno=isequal(bitand(uint32(v_init),bitcmp(uint32(mask))),bitand(uint32(v_fin),bitcmp(uint32(mask))));
            condicionDos=~isequal(bitand(uint32(v_init),bitcmp(uint32(mask))),0)&&~isequal(bitand(uint32(v_fin),bitcmp(uint32(mask))),0);
            if condicionUno||condicionDos
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%  TWO-DONUT FILL %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Creamos el primer triangulo
                v1=Tripivot(j+tipo,3);
                v2=v_init;
                v3=v_corner;
                %%%%%%%%%%%%%%%%%%%
                temp=sidetosideFill_v3(v_corner,v_init,v_corner_fin,v_fin,pasoL1,pasoL2,n_points);
                TwoDonutFill=[TwoDonutFill;v1 v2 v3;temp];
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%  TRI-DONUT FILL %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Definimos un nuevo triángulo, que será el del medio
                %primero tenemos que partir de un vertice para hallar al
                %v3_mid. Este vertice de partida será siempre del T_actual
                %para los sector 1 y 2 y para los sector3y4 será el T_next
                %notar que se está hallando el tripivot de la Donut
                %anterior (i-1)!!!!
                v3_mid=Tripivot(j+bits_sector(2),1);%es indiferente si es el vertice 1 o 2
                paso_mid=pasoL2*(-1)^xor(tipo,bits_sector(2));
                %no hace falta hacer mask ya que siempre es con la Donut
                %anterior y los Trifill ocurren a partir de la donut 2, es
                %decir i>=3
                eval_point=Point_Cloud(v3_mid+paso_mid,:);
                freepoints_mid=0;
                while(~isequal(eval_point,[0 0 0]))
                    v3_mid=v3_mid+paso_mid;
                    freepoints_mid=freepoints_mid+1;%realizado el conteo de cuantos puntos libres
                    eval_point=Point_Cloud(v3_mid+paso_mid,:);
                end
                %Hallamos el triangulo pivot Repetimos e código de líneas
                %arriba
                %------------------TRIPIVOT MIDDLE----------------%
                k_beam=(n_beams-1)*(1-bits_sector(2))+1;
                %realizamos la operación
                point=Point_Cloud(v3_mid,:);
                [v1_mid,v2_mid]=get_tripivot(point,R,rot_matrix,sector,n_beams,n_AZBLK,i-1,...
                k_beam,beam_azimuth_angles(k_beam),angle_between_azimuths,mask);
                if bits_sector(2)
                    %definimos el sentido horario Ya que para los sectors 3 y 4
                    %El sentido de los vértices es distinto a los de los
                    %primeros sector. Entonces para seguir la jerarquía de los
                    %sentidos, cambiamos aqui
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Esto podría ser o o importante                        %
                    %%% Capaz, se puede definir un sentido para un lado y otro%
                    %%% para los otros sector ()sector3 y sector4)                  %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    temp=v1_mid;
                    v1_mid=v2_mid;
                    v2_mid=temp;
                end
                T_mid=[v1_mid,v2_mid,v3_mid];
                %------------Fin TRIPIVOT MIDDLE-------------%
                %----------Primer llenado----------%
                %Definimos los v_lim Para hacer el primer llenado
                a=bits_sector(2);
                b=bits_sector(1);
                c=escalera;
                %v_lim1 será el vertice limite para el v_corner
                v_lim1=Tripivot(j+not(tipo),3);
                %v_lim2 será el vertice limite para el v_init
                index=3-double(xor(a,not(xor(b,c))))-double(not(xor(a,b)|c));
                v_lim2=T_mid(index);
                %procedemos a realizar el primer llenado
                %Creamos el primer triangulo
                v1=Tripivot(j+tipo,3);
                v2=v_init;
                v3=v_corner;
                TriDonutFill=[TriDonutFill;v1 v2 v3];
                temp=sidetosideFill_v3(v3,v2,v_lim1,v_lim2,pasoL1,pasoL2,n_points);
                TriDonutFill=[TriDonutFill;temp];
                %----------Fin Primer llenado----------%
                %----------Segundo llenado----------%
                index=double(xor(bits_sector(2),xor(bits_sector(1),escalera)));
                index=3-index-double(escalera&&not(xor(bits_sector(1),bits_sector(2))));
                %Definimos el nuevo v_init
                new_v_init=T_mid(index);
                %Creamos el triangulo de transicion
                v1=v_lim2;
                v2=new_v_init;
                v3=v_corner_fin;
                temp=sidetosideFill_v3(v3,v2,v_corner_fin,v_fin,pasoL1,pasoL2,n_points);
                TriDonutFill=[TriDonutFill;v1 v2 v3;temp];
                %----------Fin Segundo llenado----------%
                %Analizamos cuantos casos de tripivot hay. En este caso hay 16
                %count_temp=count_temp+1;
                %display(freepoints_mid)
            end
        end   
        %--------------------FIN FILL sectorS---------------------%
        %la siguiente variable guardará los tripivots pasados para
        %verificar que no exista concurrencia de tripivots con la siguiente
        %Donut
        if i~=n_donuts
            %hacemos este condicional ya que para la última Donut no
            %necesitamos realizar esto
            last_Tripivots(sector,:)=Tripivot(1+double(not(bits_sector(2)))*j,:);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%           MIDDLE FILL         %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sector1 con sector2
    vleft_init=Tripivot_middle(1,3);
    vrigth_init=Tripivot_middle(1,2);
    vleft_fin=Tripivot_middle(2,3);
    vrigth_fin=Tripivot_middle(2,1);
    pasoLeft=n_beams;
    pasoRigth=pasoLeft;
    temp=sidetosideFill_v3(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,n_points);
    TriMiddleFill=[TriMiddleFill;temp];
    %sector3 con sector4
    vleft_init=Tripivot_middle(4,3);
    vrigth_init=Tripivot_middle(4,2);
    vleft_fin=Tripivot_middle(3,3);
    vrigth_fin=Tripivot_middle(3,1);
    pasoLeft=-n_beams;
    pasoRigth=pasoLeft;
    temp=sidetosideFill_v3(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,n_points);
    TriMiddleFill=[TriMiddleFill;temp];
    if i==n_donuts &&length(Tripivot)>1
        vrigth_init=Tripivot_middle_particular(1,1);
        vleft_init=Tripivot_middle_particular(1,3);
        vrigth_fin=Tripivot_middle_particular(2,2);
        vleft_fin=Tripivot_middle_particular(2,3);
        pasoRigth=-n_beams;
        pasoLeft=-pasoRigth;
        temp=sidetosideFill_v3(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,n_points);
        TriMiddleFill=[TriMiddleFill;temp];
        vrigth_init=Tripivot_middle_particular(4,1);
        vleft_init=Tripivot_middle_particular(4,3);
        vrigth_fin=Tripivot_middle_particular(3,2);
        vleft_fin=Tripivot_middle_particular(3,3);
        pasoRigth=n_beams;
        pasoLeft=-pasoRigth;
        temp=sidetosideFill_v3(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,n_points);
        TriMiddleFill=[TriMiddleFill;temp];
    end
    %---------------Fin MIDDLE FILL-----------------%
end
Time_TDonutFill_end=toc(Time_TDonutFill_start);
total_triangles=length(OneDonutFill)+length(TwoDonutFill)+length(TriMiddleFill)+length(TriDonutFill);
fprintf ("\tSe han generado %d triángulos en el One-Donut Fill (%.4f%%)\n",length(OneDonutFill),length(OneDonutFill)*100/total_triangles);
fprintf("El tiempo de ejecución para el TwoDonutFill y el TriDonutFill es de: %.4f segundos\n",Time_TDonutFill_end);
fprintf ("\tSe han generado %d triángulos en el Two-Donut Fill (%.4f%%)\n",length(TwoDonutFill)+length(TriMiddleFill),(length(TwoDonutFill)+length(TriMiddleFill))*100/total_triangles);
fprintf ("\tSe han generado %d triángulos en el Tri-Donut Fill (%.4f%%)\n",length(TriDonutFill),length(TriDonutFill)*100/total_triangles);
fprintf ("\tEn total, se han generado %d triángulos\n",total_triangles);
fprintf("\nEl tiempo de ejecución total es de: %.4f segundos\n",Time_sphere_generation_end+Time_fusion_point_cloud_end+Time_OneDonutFill_end+Time_TDonutFill_end)

%%
figure(3)
temp=[OneDonutFill;TwoDonutFill;TriDonutFill;TriMiddleFill];
trimesh(temp,Point_Cloud(:,1),Point_Cloud(:,2),Point_Cloud(:,3))
TotalTriangleMesh=temp;
% %% RECONSTRUCCIÓN MINA
% tic
% [Point_Cloud_mina,~,~]=xlsread('dataXYZ1.csv');
% Triangle_mesh_mina=temp;
% for i=length(Triangle_mesh_mina):-1:1
%     v=Triangle_mesh_mina(i,1);
%     if isequal(Point_Cloud_mina(v,:),[0 0 0])
%         Triangle_mesh_mina(i,:)=[];
%         continue
%     end
%     v=Triangle_mesh_mina(i,2);
%     if isequal(Point_Cloud_mina(v,:),[0 0 0])
%         Triangle_mesh_mina(i,:)=[];
%         continue
%     end
%     v=Triangle_mesh_mina(i,3);
%     if isequal(Point_Cloud_mina(v,:),[0 0 0])
%         Triangle_mesh_mina(i,:)=[];
%         continue
%     end
% end
% fprintf("Tiempo para crear el mejorar el llenado:\n")
% toc
% figure (4)
% temp=Point_Cloud_mina;
% trimesh(Triangle_mesh_mina,temp(:,1),temp(:,2),temp(:,3));
% title('Superficie tridimensional')
% T=Triangle_mesh_mina;
