function out = boundary_select(pts_2D)
    if (size(unique(pts_2D(:,1)),1)~=1&&size(unique(pts_2D(:,2)),1)~=1)
        K = convhull(pts_2D(:,1),pts_2D(:,2));
        out = pts_2D(K,:);
    else
        out = pts_2D;
    end
end