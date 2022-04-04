clear;clc;close all
%% Establecemos valores iniciales del LIDAR
%nube de puntos
%Angle elevation range;
    beamAltitudeAngles=[15.379 13.236 11.128 9.03 6.941 4.878 2.788 0.705 -1.454 -3.448 -5.518 -7.601 -9.697 -11.789 -13.914 -16.062]*pi/180;
    %beam_altitude_angles=flip([-20,-17.3333,-14.6667,-12,-9.33333,-6.66667,-4,-1.33333,1.33333,4,6.66667,9.33333,12,14.6667,17.3333,20]*pi/180);
    nBeams=length(beamAltitudeAngles);
%Offset angle of AZBLK's points
    %beam_azimuth_angles=regexprep(num2str(linspace(-1.24,-0.857,16)),'\s+',',')
    beamAzimuthAngles=[-1.24,-1.2145,-1.1889,-1.1634,-1.1379,-1.1123,-1.0868,-1.0613,-1.0357,-1.0102,-0.98467,-0.95913,-0.9336,-0.90807,-0.88253,-0.857]*pi/180;
%Mode 1023 azimuth. Every azimuth increments in 88 ticks
    nAz=1024;
    %It should be noted that, we obtain 1024*16 = 16384 points in every LIDAR scan
    nPointsPerDonut=nAz*nBeams;
    angleBetweenAzimuths=-2*pi/nAz;
%Angle rotation,n° donuts,angular velocity,time per scan.
    rotAngle=-33.6*pi/180;
    rotMatrix=[cos(rotAngle) -sin(rotAngle)  0;
                sin(rotAngle) cos(rotAngle)   0;
                0               0               1];

%% Creamos funciones de utilidad
eqLine=@(m,x0,y0,x) m*(x-x0)+y0;
rotAxisX=@(alfa) [1 0         0;
                  0 cos(alfa) -sin(alfa);
                  0 sin(alfa) cos(alfa)];
rotAxisY=@(alfa) [cos(alfa)   0   sin(alfa);
                  0           1   0;
                  -sin(alfa)  0   cos(alfa)];
rotAxisZ=@(alfa) [cos(alfa) -sin(alfa)  0;
                  sin(alfa) cos(alfa)   0;
                  0           0     1];
    
%% Obtenemos la nube de puntos y la malla triangular
cloudDonut1=readmatrix("Don1.csv");
cloudDonut2=readmatrix("Don2.csv");
cloudDonut3=readmatrix("Don3.csv");
cloudDonut4=readmatrix("Don4.csv");
cloudDonut5=readmatrix("Don5.csv");
cloudDonut6=readmatrix("Don6.csv");
% meshDonut1=readmatrix("meshDonut1.csv");
% meshDonut2=readmatrix("meshDonut2.csv");
% meshDonut3=readmatrix("meshDonut3.csv");
% meshDonut4=readmatrix("meshDonut4.csv");
% meshDonut5=readmatrix("meshDonut5.csv");
% meshDonut6=readmatrix("meshDonut6.csv");
cloud=[cloudDonut1;cloudDonut2;cloudDonut3;cloudDonut4;cloudDonut5;cloudDonut6];
%cloud=cloud./vecnorm(cloud,2,2);
%mesh=[meshDonut1;meshDonut2+nPointsPerDonut;meshDonut3+nPointsPerDonut*2;meshDonut4+nPointsPerDonut*3;meshDonut5+nPointsPerDonut*4;meshDonut6+nPointsPerDonut*5];
%% Hacemos el OverlapRemoving
%definimos el plano Referencial
xRef=[sin(beamAltitudeAngles(16)) sin(beamAltitudeAngles(1))];
%Definimos los parametros de los demas planos que delimitan la zona con traslape
m=[0;tan(pi/2+rotAngle);tan(pi/2+rotAngle*2);tan(pi/2+rotAngle*3);tan(pi/2+rotAngle*4);tan(pi/2+rotAngle*5)];
x0B1=zeros(1,6);
y0B1=zeros(1,6);
x0B16=zeros(1,6);
y0B16=zeros(1,6);
for i=2:6
    point=rotMatrix^(i-1)*[sin(beamAltitudeAngles(1));-cos(beamAltitudeAngles(1));0];
    x0B1(i)=point(1);
    y0B1(i)=point(2);
    point=rotMatrix^(i-1)*[sin(beamAltitudeAngles(16));-cos(beamAltitudeAngles(16));0];
    x0B16(i)=point(1);
    y0B16(i)=point(2);
end
%%
for i=2:6%donuts
    for j=1:nPointsPerDonut
        index=(i-1)*nPointsPerDonut+j;
        %obtenmos los valores de X Y Z normalizados
        x=cloud(index,1)/norm(cloud(index,:));
        y=cloud(index,2)/norm(cloud(index,:));
        z=cloud(index,3)/norm(cloud(index,:));
        %verificamos si el punto se suprime o no
        if x<xRef(1)
            %esta a la izquierda
            if i~=2% La 2da Donut solo se elimina en la zoneRef
                %Hallamos el Y correspondiente segun la Eq. recta
                yTemp=eqLine(m(i-1),x0B16(i-1),y0B16(i-1),x);
                if y<yTemp
                    %Eliminamos el punto
                    cloud(index,:)=[0,0,0];
                end
            end
        elseif x>=xRef(1) && x<=xRef(2)
            %está en la zona Referencial
            cloud(index,:)=[0,0,0];
        else
            %está a la derecha
            if i~=2% La 2da Donut solo se elimina en la zoneRef
                %Hallamos el Y correspondiente segun la Eq. recta
                yTemp=eqLine(m(i-1),x0B1(i-1),y0B1(i-1),x);
                if y>yTemp
                    %Eliminamos el punto
                    cloud(index,:)=[0,0,0];
                end
            end
        end
    end
end
figure(1)
hold
for i=1:6
    index=(i-1)*nPointsPerDonut+1:nPointsPerDonut*i;
    scatter3(cloud(index,1),cloud(index,2),cloud(index,3))
end
%% De la malla triangular, eliminamos los triángulos que tienen puntos [0,0,0]
%lo hacemos para cada DonutMesh para facilitar el procesamiento
%se pudo haber usado fprintf y eval... pero que fastidio, preferí copiar y pegar :p
for i=length(meshDonut1):-1:1
    vTemp=meshDonut1(i,1);
    x=cloud(vTemp,1);
    y=cloud(vTemp,2);
    z=cloud(vTemp,3);
    if (x==0 && y==0&&z==0)
        %delete triangle
        meshDonut1(i,:)=[];
    else
        vTemp=meshDonut1(i,2);
        x=cloud(vTemp,1);
        y=cloud(vTemp,2);
        z=cloud(vTemp,3);
        if (x==0 && y==0&&z==0)
            %delete triangle
            meshDonut1(i,:)=[];
        else
            vTemp=meshDonut1(i,3);
            x=cloud(vTemp,1);
            y=cloud(vTemp,2);
            z=cloud(vTemp,3);
            if (x==0 && y==0&&z==0)
                %delete triangle
                meshDonut1(i,:)=[];    
            end
        end
    end
end
%donut2
for i=length(meshDonut2):-1:1
    vTemp= meshDonut2(i,1);
    x=cloud(vTemp,1);
    y=cloud(vTemp,2);
    z=cloud(vTemp,3);
    if (x==0 && y==0&&z==0)
        %delete triangle
         meshDonut2(i,:)=[];
    else
        vTemp= meshDonut2(i,2);
        x=cloud(vTemp,1);
        y=cloud(vTemp,2);
        z=cloud(vTemp,3);
        if (x==0 && y==0&&z==0)
            %delete triangle
             meshDonut2(i,:)=[];
        else
            vTemp= meshDonut2(i,3);
            x=cloud(vTemp,1);
            y=cloud(vTemp,2);
            z=cloud(vTemp,3);
            if (x==0 && y==0&&z==0)
                %delete triangle
                 meshDonut2(i,:)=[];    
            end
        end
    end
end
%donut3
for i=length(meshDonut3):-1:1
    vTemp= meshDonut3(i,1);
    x=cloud(vTemp,1);
    y=cloud(vTemp,2);
    z=cloud(vTemp,3);
    if (x==0 && y==0&&z==0)
        %delete triangle
         meshDonut3(i,:)=[];
    else
        vTemp= meshDonut3(i,2);
        x=cloud(vTemp,1);
        y=cloud(vTemp,2);
        z=cloud(vTemp,3);
        if (x==0 && y==0&&z==0)
            %delete triangle
             meshDonut3(i,:)=[];
        else
            vTemp= meshDonut3(i,3);
            x=cloud(vTemp,1);
            y=cloud(vTemp,2);
            z=cloud(vTemp,3);
            if (x==0 && y==0&&z==0)
                %delete triangle
                 meshDonut3(i,:)=[];    
            end
        end
    end
end
%donut4
for i=length(meshDonut4):-1:1
    vTemp= meshDonut4(i,1);
    x=cloud(vTemp,1);
    y=cloud(vTemp,2);
    z=cloud(vTemp,3);
    if (x==0 && y==0&&z==0)
        %delete triangle
         meshDonut4(i,:)=[];
    else
        vTemp= meshDonut4(i,2);
        x=cloud(vTemp,1);
        y=cloud(vTemp,2);
        z=cloud(vTemp,3);
        if (x==0 && y==0&&z==0)
            %delete triangle
             meshDonut4(i,:)=[];
        else
            vTemp= meshDonut4(i,3);
            x=cloud(vTemp,1);
            y=cloud(vTemp,2);
            z=cloud(vTemp,3);
            if (x==0 && y==0&&z==0)
                %delete triangle
                 meshDonut4(i,:)=[];    
            end
        end
    end
end
%donut5
for i=length(meshDonut5):-1:1
    vTemp= meshDonut5(i,1);
    x=cloud(vTemp,1);
    y=cloud(vTemp,2);
    z=cloud(vTemp,3);
    if (x==0 && y==0&&z==0)
        %delete triangle
         meshDonut5(i,:)=[];
    else
        vTemp= meshDonut5(i,2);
        x=cloud(vTemp,1);
        y=cloud(vTemp,2);
        z=cloud(vTemp,3);
        if (x==0 && y==0&&z==0)
            %delete triangle
             meshDonut5(i,:)=[];
        else
            vTemp= meshDonut5(i,3);
            x=cloud(vTemp,1);
            y=cloud(vTemp,2);
            z=cloud(vTemp,3);
            if (x==0 && y==0&&z==0)
                %delete triangle
                 meshDonut5(i,:)=[];    
            end
        end
    end
end
%donut6
for i=length(meshDonut6):-1:1
    vTemp= meshDonut6(i,1);
    x=cloud(vTemp,1);
    y=cloud(vTemp,2);
    z=cloud(vTemp,3);
    if (x==0 && y==0&&z==0)
        %delete triangle
         meshDonut6(i,:)=[];
    else
        vTemp= meshDonut6(i,2);
        x=cloud(vTemp,1);
        y=cloud(vTemp,2);
        z=cloud(vTemp,3);
        if (x==0 && y==0&&z==0)
            %delete triangle
             meshDonut6(i,:)=[];
        else
            vTemp= meshDonut6(i,3);
            x=cloud(vTemp,1);
            y=cloud(vTemp,2);
            z=cloud(vTemp,3);
            if (x==0 && y==0&&z==0)
                %delete triangle
                 meshDonut6(i,:)=[];    
            end
        end
    end
end
%% Aqui realizamos la union de los meshes(FineMesh)
%El primer mesh no se toca. En cambio, para el segundo mesh hacemos la
%misma lógica del FineMesh algorithm. 
%-Primero es ubicar el punto extremo distinto de cero (eso está listo). 
%-Luego, hallar el punto no nulo de la Donut anterior. 
%-Luego, hallar el otro vertice revisando la malla triangular respectiva.
%-Hacer esto mismo, para hallar todos los triángulos pivots.
%-Luego usar el algortimo de FillZone, el cual debe modificarse, ya que los
% puntos entre tripivots (freepoints) pueden variar debido al Downsampling
%-Para el caso de TriDonutFill, establecer bien los puntos de inicio y fin
 
%En teoría, hay que hacer unas breves modificaciones.. pero a primera
%impresión pareciera que funcionará.