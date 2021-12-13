function T=sidetosideFill (vleft_init,vrigth_init,vleft_fin,vrigth_fin,pasoLeft,pasoRigth,i,n_points,mask)
    freepointsL1=(vleft_fin-vleft_init)/pasoLeft;
    freepointsL2=(vrigth_fin-vrigth_init)/pasoRigth;
    v3=vleft_init;
    v2=vrigth_init;
    v3_fin=vleft_fin;
    v2_fin=vrigth_fin;
    pasov3=pasoLeft;
    pasov2=pasoRigth;
    %estoy considerando v3 a la izquierda y v2 a la derecha
    T=[];
    volteamos=false;
    if freepointsL2~=freepointsL1
       %Realizamos un triangulo con mismo vertice
       arista_mismo_vex=freepointsL2-freepointsL1;
       if arista_mismo_vex<0
           %En caso L1 tenga mas puntos el v_comun estará en L2
           v3=vrigth_init;
           v2=vleft_init;
           pasov2=pasoLeft;
           arista_mismo_vex=arista_mismo_vex*-1;
           volteamos=true;
       end
       for j=1:arista_mismo_vex
          v1=v2;
          v2=bitand(v2+pasov2-1,mask)+1+(i-1-not(volteamos))*n_points;
          T=[T;v1,v2,v3];
       end
       if volteamos
           temp=v3;
           v3=v2;
           v2=temp;
           pasov2=pasoRigth;
       end
    end
    while (v2~=v2_fin)
        v1=v2;
        v2=bitand(v3+pasov3-1,mask)+1+(i-1)*n_points;
        T=[T;v1,v2,v3];
        v3=v2;
        v2=bitand(v1+pasov2-1,mask)+1+(i-2)*n_points;
        T=[T;v1,v2,v3];
    end
end