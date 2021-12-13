function resultado=check_if_PC_zero(T,point_cloud)
    temp=zeros(length(T)*3,3);
    [row,col]=size(T);
    for i=1:row
        point1=point_cloud(T(i,1),:);
        point2=point_cloud(T(i,2),:);
        point3=point_cloud(T(i,3),:);
        temp((i-1)*3+1:(i-1)*3+3,:)=[point1;point2;point3];
    end
    resultado=any(any(temp==0));
end