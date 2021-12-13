function [vmin,vmax]=get_tripivot(point,R,rot_matrix,slot,n_beams,n_AZBLK,i,k_beam,azimuth_angle,angle_between_azimuths,mask)
    %Se debe tener en cuenta los rango de acos y asin
    %acos: [0~pi]
    %atan: [-pi/2~pi/2]
    bits_slot=de2bi(slot-1,2);
    x_data=point(1);
    thetha=-acos(x_data/R);
    rot_point=rot_matrix^-(i-2)*point';
    rot_thetha=-acos(rot_point(1)/R);
    if thetha*(-1)^bits_slot(2)>=rot_thetha*(-1)^bits_slot(2)
        %los vértices pertenecerán a la Donut referencial
        y_data=point(2);
        z_data=point(3);
        offset=0;
    else
        %los vértices pertenecerán a la previa Donut
        y_data=rot_point(2);
        z_data=rot_point(3);
        thetha=rot_thetha;
        offset=n_AZBLK*n_beams*(i-2);
    end
    %Hallamos el ángulo alfa, que es el angulo del azimuth
    alfa=-pi/2 +atan(z_data/y_data);
    if y_data>0
        alfa=-pi+alfa;
    end
    %Añadimos -2pi al alfa para no tener probelmas con el bitand,
    %esto al final no perjudica ya que se hace el masking de bits
    %solo que tener en cuenta que la división entre -2pi/ang_bet_azit
    %da un total de 1024.
    alfa=alfa-2*pi;
    %nos hemos asegurado que el alfa sea siempre negativo
    %Calculamos el azimuth que corresponde al alfa
    k_azimuth=floor((alfa-azimuth_angle)/angle_between_azimuths);
    vmin=(k_azimuth)*n_beams+k_beam;
    vmin=bitand(vmin-1,mask)+1;%Evitamos problemas con el mask para los multiplos
    vmax=bitand(vmin+n_beams-1,mask)+1;
    %Realizamos offset y sentido
    vmin=vmin+offset;
    vmax=vmax+offset;
end


% %%%%--------%Hallamos el triángulo pivot
% %%%%Inicia código tripivot
%              %Point_Cloud
%             %point
%             %rot_matrix
%             %slot
%             %n_points
%             %i
%             %beam_azimuth_angle
%             %angle_between_azimuth
%             %mask
%             point=Point_Cloud(v3,:);
%             %Se debe tener en cuenta los rango de acos y asin
%             %acos: [0~pi]
%             %atan: [-pi/2~pi/2]
%             x_data=point(1);
%             thetha=-acos(x_data/R);
%             rot_point=rot_matrix^-(i-2)*point';
%             rot_thetha=-acos(rot_point(1)/R);
%             if thetha*(-1)^bits_slot(2)>=rot_thetha*(-1)^bits_slot(2)
%                 %los vértices pertenecerán a la Donut referencial
%                 y_data=point(2);
%                 z_data=point(3);
%                 offset=0;
%             else
%                 %los vértices pertenecerán a la previa Donut
%                 y_data=rot_point(2);
%                 z_data=rot_point(3);
%                 thetha=rot_thetha;
%                 offset=n_AZBLK*n_beams*(i-2);
%             end
%             %Hallamos el ángulo alfa, que es el angulo del azimuth
%             alfa=-pi/2 +atan(z_data/y_data);
%             if y_data>0
%                 alfa=-pi+alfa;
%             end
%             %Añadimos -2pi al alfa para no tener probelmas con el bitand,
%             %esto al final no perjudica ya que se hace el masking de bits
%             %solo que tener en cuenta que la división entre -2pi/ang_bet_azit
%             %da un total de 1024.
%             alfa=alfa-2*pi;
%             %nos hemos asegurado que el alfa sea siempre negativo
%             %Calculamos el azimuth que corresponde al alfa
%             %segun el slot, trabajamos con un cierto beam_elevate_angle
%             k_beam=15*(1-bits_slot(2))+1;
%             %realizamos la operación
%             k_azimuth=floor((alfa-beam_azimuth_angles(k_beam))/angle_between_azimuths);
%             v1=(k_azimuth)*n_beams+k_beam;
%             v1=bitand(v1-1,mask)+1;%Evitamos problemas con el mask para los multiplos
%             v2=bitand(v1+n_beams,mask);%No es necesario restar en el mask y luego sumar
%             %Realizamos offset y sentido
%             v1=v1+offset;
%             v2=v2+offset;
%             if bits_slot(2)
%                 %definimos el sentido horario Ya que para los slots 3 y 4
%                 %El sentido de los vértices es distinto a los de los
%                 %primeros slot. Entonces para seguir la jerarquía de los
%                 %sentidos, cambiamos aqui
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%% Esto podría ser o o importante                        %
%                 %%% Capaz, se puede definir un sentido para un lado y otro%
%                 %%% para los otros slot ()slot3 y slot4)                  %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 temp=v1;
%                 v1=v2;
%                 v2=temp;
%             end
%             Tripivot(n_tripivot,:)=[v1 v2 v3];
% %%%%%%Fin codígo tripivot