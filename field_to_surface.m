function S = field_to_surface(thres,field,grid)
S = [];
for i=1:size(field,1)
    %if (field(i,1)>=0&&field(i,1)<=thres)
    if (field(i,1)<=thres)
        S = [S;grid(i,:)];
    end
end
end