function T=sidetosideFill_v3(vL1_init,vL2_init,vL1_fin,vL2_fin,pasoL1,pasoL2,n_points) 
    mask=uint32(n_points-1);
    %Realizamos algunos ajustes para poderr obterner el valor magnitud de
    %los puntos libres tanto para la izquierda y derecha
    if pasoL1<0
        %...Me parece que hay que restar y volver a sumar una unidad....
        vmax=bitand(vL1_init,mask);
        vmin=bitand(vL1_fin,mask);
    else
        vmax=bitand(vL1_fin,mask);
        vmin=bitand(vL1_init,mask);
    end
    freepointsL1=double(bitand((vmax+(mask-vmin)+1),mask)/abs(pasoL1));
    if pasoL2<0
        vmax=bitand(vL2_init,mask);
        vmin=bitand(vL2_fin,mask);
    else
        vmax=bitand(vL2_fin,mask);
        vmin=bitand(vL2_init,mask);
    end
    freepointsL2=double(bitand((vmax+(mask-vmin)+1),mask)/abs(pasoL2));
    %definimos unas variables adicionales para poder realizar una mejor
    %amanipulación
    offsetL1=n_points*floor(double((vL1_init-1))/n_points);
    offsetL2=n_points*floor(double((vL2_init-1))/n_points);
    %considero que "v3 es de L1" y "v2 es de L2"
    v3=vL1_init;
    v2=vL2_init;
    v3_fin=vL1_fin;
    v2_fin=vL2_fin;
    pasov3=pasoL1;
    pasov2=pasoL2;
    T=[];
    %la siguiente variable es para el caso que el vertice en comun este en L2
    volteamos=false;
    if freepointsL2~=freepointsL1
       %Realizamos un triangulo con mismo vertice
       arista_mismo_vex=freepointsL2-freepointsL1;
       offset=offsetL2;
       if arista_mismo_vex<0
           %En caso L1 tenga mas puntos el v_comun estará en L2
           v3=vL2_init;
           v2=vL1_init;
           pasov2=pasoL1;
           offset=offsetL1;
           arista_mismo_vex=arista_mismo_vex*-1;
           volteamos=true;
       end
       for j=1:arista_mismo_vex
          v1=v2;
          %debemos realizar el offset adecuado segun la Donut con la que
          %trabajamos 
          v2=bitand((v2+mask)+pasov2,mask)+1+offset;
          T=[T;v1,v2,v3];
       end
       if volteamos
           %en caso habiamos volteado, volvemos al caso inicial
           temp=v3;
           v3=v2;
           v2=temp;
           pasov2=pasoL2;
       end
    end
    while (v2~=v2_fin)
        v1=v2;
        %Con esta formula podemos obtener la concatenación de Donuts
        v2=bitand((v3+mask)+pasov3,mask)+1+offsetL1;
        T=[T;v1,v2,v3];
        v3=v2;
        %Esta formula permite avanzar por medio del enmascaramiento sin
        %necesidad de usar condicionales (equivalente al operador modulo)
        v2=bitand((v1+mask)+pasov2,mask)+1+offsetL2;
        T=[T;v1,v2,v3];
    end
end