function all_points_modified=particular_surfaces_allowed_v1(points_data,up_point,down_point)
    %all_points_modified=zeros(length(points_data{1,1})*length(points_data),1);
    all_points_modified=[];
    for i=1:length(points_data)
        all_points_modified=[all_points_modified;
                             points_data{1,i}];
    end
    R=1;
    p1=[R*cos((-31*4-16.062)*pi/180) R*sin((-31*4-16.062)*pi/180)];
    p0=[(R*sin(-16.062*pi/180)-p1(2))/tan(-31*4*pi/180)+p1(1) R*sin(-16.062*pi/180)];
    p0(3)=sqrt(R*R-p0(1)*p0(1)-p0(2)*p0(2));
    p2=[-R*cos(-16.062*pi/180) R*sin(-16.062*pi/180)];
    for i=1:length(all_points_modified)
        x=all_points_modified(i,1);
        y=all_points_modified(i,2);
        z=all_points_modified(i,3);
        if (abs(z)>p0(3)) || x>p0(1)|| x<p2(1)
            all_points_modified(i,:)=[0,0,0];
        elseif x>p1(1)
            temp_y=tan(-31*4*180/pi)*(x-p1(1))+p1(2);
            if y<temp_y || y>R*sin(-16.062*pi/180)
                all_points_modified(i,:)=[0,0,0];
            end
        else
            temp_y=-sqrt(R*R-x*x);
            if y<temp_y || y>R*sin(-16.062*pi/180)
                all_points_modified(i,:)=[0,0,0];
            end
        end
    end
    all_points_modified(all_points_modified(:,1)==0,:) = [] ;
    rot_180_matrix=[-1 0 0;
                    0 -1 0;
                    0  0 1];
    all_points_modified=[all_points_modified;(rot_180_matrix*all_points_modified')'];
end