function all_points_modified=surfaces_allowed_v1(points_data,up_points,down_points)
    %all_points_modified=zeros(length(points_data{1,1})*length(points_data),1);
    all_points_modified=[];
    for i=1:length(points_data)
        all_points_modified=[all_points_modified;
                             points_data{1,i}];
    end
    
    for i=1:length(all_points_modified)
        x=all_points_modified(i,1);
        y=all_points_modified(i,2);
        z=all_points_modified(i,3);
        for j=1:1%length up_points
            %1st surface 
            if z>up_points(j,3) || z<down_points(j,3)
                all_points_modified(i,:)=[0,0,0];
            elseif x<up_points(j,1) || x> 1*cos((16.062)*pi/180)
                all_points_modified(i,:)=[0,0,0];
            elseif x <1*cos((15.379-62)*pi/180)
                temp_y=tan(-62*pi/180)*(x-up_points(j,1))+up_points(j,2);
                if y>1*sin(-16.062*pi/180)||y<temp_y
                    all_points_modified(i,:)=[0,0,0];
                end
            else
                temp_y=-sqrt(1-x*x);
                if y>1*sin(-16.062*pi/180)||y<temp_y
                    all_points_modified(i,:)=[0,0,0];
                end
            end
        end
        
    end
end