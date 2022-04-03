function [v1,v2]=get_azimuth_vex_v2(point,i,R,angle_between_azimuths,n_AZBLK,n_rays,rot_matrix,k_ray,k_offset_azimuth,minimal_offset)
    %ver IDEA sem7, esta función permite obtener el ángulo alfa que se forma
    %con la vertical y el azimuth que contiene a los dos vertices del triangulo
    %pivot.
    %El circulo/anillo está en el plano xz
    %Recordar que acos:[0,pi] y asin:[-pi/2,pi/2]
    %Entradas: 
            %  point: vector fila que contiene las coordenadas del vertice 
            %         v3 del triángulo pivot
            %  i:Donut a la cual pertenecen los puntos
            %  R: Radio de la esfers
            %  angle_between_azimuths: Ángulo entre los azimuts
            %  n_AZBLK: numero de azimuts
            %  n_rays:  numero de rayos
            %  rot_matriz: matriz que hace la rotación de rot_angle
            %  k_ray:   indicar si es el rayo 1 o 16
            %  k_offset_azimuth: se establece la cantidad de azimuts
            %                    a desplazar debido a que el offset del
            %                    rayo es mayor que el angulo entre los
            %                    azimuts
            %  minimal_offset: el ángulo entre el eje del k_azimuth y el
            %                  rayo. En sí es la resta del AZBLK_offset y
            %                  el angle_betwen_azimuths
            %Salidas:
            %v1,ve: vertices que serán parte del triángulo pívot
    
    x_data=point(1);
    thetha=-acos(x_data/R);
    rot_point=rot_matrix^-(i-2)*point';
    rot_thetha=-acos(rot_point(1)/R);
   if x_data<0
        if thetha>=rot_thetha
            %los vértices pertenecerán a la Donut referencial
            y_data=point(2);
            z_data=point(3);
            offset=0;
        else
            %los vértices pertenecerán a la previa Donut
            y_data=rot_point(2);
            z_data=rot_point(3);
            thetha=rot_thetha;
            offset=n_AZBLK*n_rays*(i-2);
         end
    else
        if thetha<=rot_thetha
            %los vértices pertenecerán a la Donut referencial
            y_data=point(2);
            z_data=point(3);
            offset=0;
        else
            %los vértices pertenecerán a la previa Donut
            y_data=rot_point(2);
            z_data=rot_point(3);
            thetha=rot_thetha;
            offset=n_AZBLK*n_rays*(i-2);
        end
    end

    %Hallamos el ángulo alfa, que es el angulo del azimuth
    alfa=asin(-y_data/R/sin(thetha));
    if z_data>0
        alfa=-pi-alfa;
    end
    %Calculamos el azimuth que corresponde al alfa
    k_azimuth=floor(alfa/angle_between_azimuths)+1-k_offset_azimuth;
    if mod(alfa,angle_between_azimuths)>minimal_offset
        %caso1
        k_azimuth=k_azimuth-1; %transformamos el caso1 a caso2, es decir, 
                               %lo hacemos mas general.
    end   
    if k_azimuth<=0
        %En caso el resultado sea negativo
        k_azimuth=k_azimuth+n_AZBLK;
    end
    if k_azimuth==n_AZBLK 
        %Caso particular donde se concatenan los azimuts
        v1=n_rays*(n_AZBLK-1)+k_ray;
        v2=k_ray;
    else
        %Caso general
        v1=(k_azimuth-1)*n_rays+k_ray;
        v2=v1+n_rays;
    end
    %Hacemos el offset segun la Donut que pertenecen los vértices
    v1=v1+offset;
    v2=v2+offset;
    if x_data>0
        %definimos el sentido horario
        temp=v1;
        v1=v2;
        v2=temp;
    end
end

