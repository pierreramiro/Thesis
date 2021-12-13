function y_n=xy_curve_eq_v1(elev_angle_1,elev_angle_2,i,rot_angle,xp,yp,zp,cross_point,c1,c2,c3,c4)
    R=1;
    c1=sin(elev_angle);
    c2=cos(elev_angle);
    c3=cos(rot_angle*(i-1));
    c4=sin(rot_angle*(i-1));
        x=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c3 + (c1*sen_alfa+c2*sen_phi*cos_alfa)*c4 );
        y=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c4 - (c1*sen_alfa+c2*sen_phi*cos_alfa)*c3 );    
        z=-R*c2*cos_phi;
        
        
        if yp>0
            phi=-acos(-zp/R/c2);
        else
            phi=-pi-acos(-zp/R/c2);
        end
        alfa=-rot_angle*phi/(2*pi);
        sen_phi=sin(phi);
        sen_alfa=sin(alfa);
        cos_alfa=cos(alfa);
        cos_phi=cos(phi);
        x=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c3 + (c1*sen_alfa+c2*sen_phi*cos_alfa)*c4 );
        y=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c4 - (c1*sen_alfa+c2*sen_phi*cos_alfa)*c3 );    
        z=-R*c2*cos_phi;
    else %statement for Curve-2
        
    end
    
end
%%With the next script, we can see the line plotted
% figure
% R=1;
% rot_angle=-0.5853;
% phi=-2*pi/0.1*(0:0.0001:0.1);
% elev_angle=AZBLK_angle(end);
% alfa=-rot_angle/2/pi*phi;
% x=-R*sin(elev_angle).*cos(alfa)+R*cos(elev_angle).*sin(phi).*sin(alfa);
% y=-R*sin(elev_angle).*sin(alfa)-R*cos(elev_angle).*sin(phi).*cos(alfa);
% z=-R*cos(elev_angle).*cos(phi);
% 
% plot3(x,y,z)
% elev_angle=AZBLK_angle(1);
% x=-R*sin(elev_angle).*cos(alfa)+R*cos(elev_angle).*sin(phi).*sin(alfa);
% y=-R*sin(elev_angle).*sin(alfa)-R*cos(elev_angle).*sin(phi).*cos(alfa);
% z=-R*cos(elev_angle).*cos(phi);
% hold
% plot3(x,y,z)