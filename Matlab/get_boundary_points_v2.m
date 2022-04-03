function boundary_points=get_boundary_points_v2(point_cloud_Donut_i_without_overlap,n_rays,n_AZBLK)    
% Entradas
%     point_cloud_Donut_i_without_overlap: vector que contiene los puntos 
%                                           de la Donut, sin traslape
%     n_rays: Cantidad de rayos del LiDAR
%     n_AZBLK: Cantidad de azimut del LiDAR
% Salidas
%   boundary_points: celda de (1,n_rays) que contiene a el vertice de los
%                    puntos extremos y las coordenadas de estos

    %preallocate boundary_points
    boundary_points=cell(1,n_rays);
    %La primera fila contendrá la cantidad de puntos extremos. Dado que es
    %una dona que será cortada por otra dona, solo existirán como máximo 4
    %puntos que estarán desde la fila 2 a la 5ta, siendo el primer término
    %el indice del punto
    boundary_points(:)={zeros(4,4)};
    n_points=n_rays*n_AZBLK;
    for j=1:n_rays
        %%%%%%    PRIMER CASO    %%%%%%%
        %definimos los puntos actuales y previos
        prev_point=point_cloud_Donut_i_without_overlap(n_points-n_rays+j,:);
        act_point=point_cloud_Donut_i_without_overlap(j,:);
        %obtenemos los valores lógicos
        A=isequal(prev_point,[0,0,0]);%el punto anterior es igual a cero
        B=isequal(act_point,[0,0,0]);%el punto actual es igual a cero
        %realizamos los condicionales. El primer caso considera la
        %concatenación de las donas, asi que tiene que estar fuera del "for"
        if ~A && B
            %Flanco de bajada
            %hallamos el vertice
            vertice=n_points-n_rays+j;
            %hallamos el punto
            point=point_cloud_Donut_i_without_overlap(vertice,:);
            %procedemos a clasificar el vertice según el punto
            %evaluamos el x
            if point(1)<=0
                %evaluamos el z
                if point(3)<=0
                    slot=1;
                else
                    slot=2;
                end
            else
                %evaluamos el z
                if point(3)>=0
                    slot=3;
                else
                    slot=4;
                end
            end
            %guardamos el valor al vector
            boundary_points{j}(slot,:)=[vertice prev_point];
        elseif A && ~B
            %Flanco de subida
            %hallamos el vertice
            vertice=j;
            %hallamos el punto
            point=point_cloud_Donut_i_without_overlap(vertice,:);
            %procedemos a clasificar el vertice según el punto
            %evaluamos el x
            if point(1)<=0
                %evaluamos el z
                if point(3)<=0
                    slot=1;
                else
                    slot=2;
                end
            else
                %evaluamos el z
                if point(3)>=0
                    slot=3;
                else
                    slot=4;
                end
            end
            %guardamos el valor al vector
            boundary_points{j}(slot,:)=[vertice act_point];
        end
        %%%%%%    SIGUIENTES CASOS    %%%%%%%
        for k=2:n_AZBLK
            vertice=(k-1)*n_rays+j;
            prev_point=point_cloud_Donut_i_without_overlap(vertice-16,:);
            act_point=point_cloud_Donut_i_without_overlap(vertice,:);
            A=isequal(prev_point,[0,0,0]);
            B=isequal(act_point,[0,0,0]);
            if ~A && B
                %Flanco de bajada
                %hallamos el punto
                point=point_cloud_Donut_i_without_overlap(vertice-16,:);
                %procedemos a clasificar el vertice según el punto
                %evaluamos el x
                if point(1)<=0
                    %evaluamos el z
                    if point(3)<=0
                        slot=1;
                    else
                        slot=2;
                    end
                else
                    %evaluamos el z
                    if point(3)>=0
                        slot=3;
                    else
                        slot=4;
                    end
                end
                %guardamos el valor al vector
                boundary_points{j}(slot,:)=[vertice-16 prev_point];
            elseif A && ~B
                %Flanco de subida
                %hallamos el punto
                point=point_cloud_Donut_i_without_overlap(vertice,:);
                %procedemos a clasificar el vertice según el punto
                %evaluamos el x
                if point(1)<=0
                    %evaluamos el z
                    if point(3)<=0
                        slot=1;
                    else
                        slot=2;
                    end
                else
                    %evaluamos el z
                    if point(3)>=0
                        slot=3;
                    else
                        slot=4;
                    end
                end
                %guardamos el valor al vector
                boundary_points{j}(slot,:)=[vertice act_point];
            end
        end
    end
end