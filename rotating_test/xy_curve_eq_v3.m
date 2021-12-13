function [x,y]=xy_curve_eq_v3(elev_angle,i,rot_angle,zp)
    R=1;
    c1=sin(elev_angle);
    c2=cos(elev_angle);
    c3=cos(rot_angle*(i-1));
    c4=sin(rot_angle*(i-1));
    phi=-2*pi+acos(-zp/R/c2);
    alfa=-rot_angle*phi/(2*pi);
    sen_phi=sin(phi);
    sen_alfa=sin(alfa);
    cos_alfa=cos(alfa);
    cos_phi=cos(phi);
    x=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c3 + (c1*sen_alfa+c2*sen_phi*cos_alfa)*c4 );
    y=R*( (-c1*cos_alfa+c2*sen_phi*sen_alfa)*c4 - (c1*sen_alfa+c2*sen_phi*cos_alfa)*c3 );    
    z=-R*c2*cos_phi;
        
end