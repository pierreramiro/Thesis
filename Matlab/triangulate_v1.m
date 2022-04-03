function T=triangulate_v1(temp,n_elevation_points,n_AZBLK)
    n_triangles_per_column=(n_elevation_points-1)*2;
    T_default=zeros((n_elevation_points-1)*(n_AZBLK)*2,3);
    for i=1:n_AZBLK-1
        for j=1:n_elevation_points-1
            T_default(n_triangles_per_column*(i-1)+j*2-1,:)=[ n_elevation_points*(i-1)+j
                                                            n_elevation_points*i+j
                                                            n_elevation_points*i+j+1];
            T_default(n_triangles_per_column*(i-1)+j*2,:)=[ n_elevation_points*(i-1)+j
                                                          n_elevation_points*i+j+1
                                                          n_elevation_points*(i-1)+j+1];
        end
    end
    i=n_AZBLK;
    for j=1:n_elevation_points-1
        T_default(n_triangles_per_column*(i-1)+j*2-1,:)=[ n_elevation_points*(i-1)+j
                                                            j
                                                            j+1];
        T_default(n_triangles_per_column*(i-1)+j*2,:)=[ n_elevation_points*(i-1)+j
                                                          j+1
                                                          n_elevation_points*(i-1)+j+1];
    end
    
    T=T_default;
    
    for i=length(T):-1:1
        if isequal(temp(T(i,1),:),[0,0,0])
            T(i,:)=[];
        elseif isequal(temp(T(i,2),:),[0,0,0])
            T(i,:)=[];
        elseif isequal(temp(T(i,3),:),[0,0,0])
            T(i,:)=[];
        end
    end    
end