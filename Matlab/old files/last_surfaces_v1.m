function all_points_modified= last_surfaces_v1(points_data)
    %all_points_modified=zeros(length(points_data{1,1})*length(points_data),1);
    all_points_modified=[];
    for i=1:length(points_data)
        all_points_modified=[all_points_modified;
                             points_data{1,i}];
    end
    R=1;
    p1=[R*cos((-31*5-16.062)*pi/180) R*sin((-31*5-16.062)*pi/180)];
    p2=[R*cos((-180-31+15.379)*pi/180) R*sin((-180-31+15.379)*pi/180)];
    syms x
    eqn= tan(-31*5*pi/180)*(x-p1(1))+p1(2)==tan(-31*pi/180)*(x-p2(1))+p2(2);
    p0(1)=eval(solve(eqn,x));
    p0(2)=tan(-31*5*pi/180)*(p0(1)-p1(1))+p1(2);
    p0(3)=sqrt(R*R-p0(1)*p0(1)-p0(2)*p0(2));
    
    for i=1:length(all_points_modified)
        x=all_points_modified(i,1);
        y=all_points_modified(i,2);
        z=all_points_modified(i,3);
        if abs(z)>p0(3)|| y>p2(2) || y<p1(2)
            all_points_modified(i,:)=[0,0,0];
        elseif y>p0(2)
            temp_x=(y-p2(2))/tan(-31*pi/180)+p2(1);
            if x>temp_x || x<-sqrt(R*R-y*y)
                all_points_modified(i,:)=[0,0,0];
            end
        else
            temp_x=(y-p1(2))/tan(-31*5*pi/180)+p1(1);
            if x>temp_x ||x<-sqrt(R*R-y*y)
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