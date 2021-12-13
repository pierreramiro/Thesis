function [v1,v2]=get_azimuth_vex(x_data,y_data,z_data,R,angle_between_azimuths,n_AZBLK,n_rays,k_ray,k_offset_azimuth,minimal_offset)
    %ver IDEA sem7, esta función permite obtener el ángulo alfa que se forma
    %con la vertical y el azimuth que contiene a los dos vertices del triangulo
    %pivot.
    %El circulo/anillo está en el plano xz
    %Recordar que acos:[0,pi] y asin:[-pi/2,pi/2]
    %Entradas: x_data,
            %  y_data,
            %  z_data,
            %  R,
            %  angle_between_azimuths,
            %  n_AZBLK:numero de azimuts
            %  n_rays: numero de rayos
            %  k_ray: indicar si es el rayo 1 o 16
            %  k_offset_azimuth: se establece la cantidad de azimuts
            %                    a desplazar debido a que el offset del
            %                    rayo es mayor que el angulo entre los
            %                    azimuts
            %  minimal_offset: el ángulo offset, entre los azimuts, debido
            %                  al offset de los rayos. Este es menor que el
            %                  ángulo entre los azimuts
    thetha=-acos(x_data/R);
    alfa=asin(-y_data/R/sin(thetha));
    if z_data>0
        alfa=-pi-alfa;
    end
    k_azimuth=floor(alfa/angle_between_azimuths)+1-k_offset_azimuth;
    if mod(alfa,angle_between_azimuths)>minimal_offset
        %caso1
        k_azimuth=k_azimuth-1; %transformamos el caso1 a caso2
    end   
    if k_azimuth<=0
        k_azimuth=k_azimuth+n_AZBLK;
    end
    if k_azimuth==n_AZBLK 
        v1=n_rays*(n_AZBLK-1)+k_ray;
        v2=k_ray;
    else
        v1=(k_azimuth-1)*n_rays+k_ray;
        v2=v1+n_rays;
    end

end

