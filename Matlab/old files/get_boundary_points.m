function boundary_points=get_boundary_points(points_Donut,n_rays,n_AZBLK)    
% Entradas
%     points_Donut: vector que contiene los puntos de la Donut, sin
%                   traslape
%     n_rays: Cantidad de rayos del LiDAR
%     n_AZBLK: Cantidad de azimut del LiDAR
% Salidas
%   boundary_points: celda de (1,n_rays) que contiene a la cantidad de
%                    puntos extremos y los valores de estos

    %preallocate boundary_points
    boundary_points=cell(1,n_rays);
    %La primera fila contendrá la cantidad de puntos extremos. Dado que es
    %una dona que será cortada por otra dona, solo existirán como máximo 4
    %puntos que estarán desde la fila 2 a la 5ta, siendo el primer término
    %el indice del punto
    boundary_points(:)={zeros(5,4)};
    for j=1:n_rays
        %resetamos el contador
        count=0;
        %definimos los puntos actuales y previos
        act_point=points_Donut(j,:);
        prev_point=points_Donut((n_AZBLK-1)*n_rays+j,:);
        %obtenemos los valores lógicos
        A=isequal(act_point,[0,0,0]);
        B=isequal(prev_point,[0,0,0]);
        %realizamos los condicionales. El primer caso considera la
        %concatenación de las donas, asi que tiene que estar fuera del "for"
        if A && ~B
            count=count+1;
            boundary_points{j}(count+1,:)=[j act_point];
        elseif ~A && B
            count=count+1;
            boundary_points{j}(count+1,:)=[j prev_point];
        end
        for k=2:n_AZBLK
            index=(k-1)*n_rays+j;
            act_point=points_Donut(index,:);
            prev_point=points_Donut(index-16,:);
            A=isequal(act_point,[0,0,0]);
            B=isequal(prev_point,[0,0,0]);
            if A && ~B
                count=count+1;
                boundary_points{j}(count+1,:)=[index-16 prev_point];
            elseif ~A && B
                count=count+1;
                boundary_points{j}(count+1,:)=[index act_point];
            end
        end
        %guardamos en la primera fila el valor de la cantidad de puntos
        %extremos
        boundary_points{j}(1,1)=count;
    end
end