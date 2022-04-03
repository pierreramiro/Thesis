function T=sidetosideFill_v2(vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,i,n_points,mask,Donut_ref) 
    %Realizamos algunos ajustes para poderr obterner el valor magnitud de
    %los puntos libres tnaot para la izquierda y derecha
    if pasoLeft<0
        vmax=bitand(vleft_init,mask);
        vmin=bitand(vleft_fin,mask);
    else
        vmax=bitand(vleft_fin,mask);
        vmin=bitand(vleft_init,mask);
    end
    freepointsLeft=bitand((vmax+(mask-vmin)+1),mask)/abs(pasoLeft);
    if pasoRigth<0
        vmax=bitand(vrigth_init,mask);
        vmin=bitand(vrigth_fin,mask);
    else
        vmax=bitand(vrigth_fin,mask);
        vmin=bitand(vrigth_init,mask);
    end
    freepointsRigth=bitand((vmax+(mask-vmin)+1),mask)/abs(pasoRigth);
    %definimos unas variables adicionales para poder realizar una mejor
    %amanipulación
    v3=vleft_init;
    v2=vrigth_init;
    v3_fin=vleft_fin;
    v2_fin=vrigth_fin;
    pasov3=pasoLeft;
    pasov2=pasoRigth;
    %estoy considerando v3 a la izquierda y v2 a la derecha
    T=[];
    %la siguiente variable es para el caso que el vertice en comun este a
    %la derecha
    volteamos=false;
    if freepointsRigth~=freepointsLeft
       %Realizamos un triangulo con mismo vertice
       arista_mismo_vex=double(freepointsRigth)-double(freepointsLeft);
       if arista_mismo_vex<0
           %En caso Left tenga mas puntos el v_comun estará en Rigth
           v3=vrigth_init;
           v2=vleft_init;
           pasov2=pasoLeft;
           arista_mismo_vex=arista_mismo_vex*-1;
           volteamos=true;
       end
       for j=1:arista_mismo_vex
          v1=v2;
          %debemos realizar el offset adecuado segun la Donut con la que
          %trabajamos (i=Donut izquierda)
          %tambien de tiene una variable de Donut_ref, que indica que el
          %offset es cero, ya que estamos con la donut referencial
          offset=(not(Donut_ref)||volteamos)*((i-1-not(volteamos))*n_points);
          v2=bitand(v2+pasov2+mask,mask)+1+offset;
          T=[T;v1,v2,v3];
       end
       if volteamos
           %en caso habiamos volteado, volvemos al caso inicial
           temp=v3;
           v3=v2;
           v2=temp;
           pasov2=pasoRigth;
       end
    end
    while (v2~=v2_fin)
        v1=v2;
        %Con esta formula podemos obtener la concatenación de Donuts
        v2=bitand(v3+pasov3+mask,mask)+1+(i-1)*n_points;
        T=[T;v1,v2,v3];
        v3=v2;
        %realizamos el offset de acuerdo a si es una Donut previa (i-1) o
        %es la Donut refrencial, offset=0
        offset=(i-2)*n_points*not(Donut_ref);
        %Esta formula permite avanzar por medio del enmascaramiento sin
        %necesidad de usar condicionales (equivalente al operador modulo)
        v2=bitand(v1+pasov2+mask,mask)+1+offset;
        T=[T;v1,v2,v3];
    end
end