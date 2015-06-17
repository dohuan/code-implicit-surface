function out = surface_refiner(points)
%% Use k-mean clustering and convex hull to extract ONLY the surface points
% points: 3-D coordinates
z_line = unique(points(:,end),'rows','stable');
n = size(points,1);
nz = size(z_line,1);
for i=1:nz
    pts(i).z = z_line(i);
    pts(i).xy = [];
end
for i=1:n
    ix = find(z_line==points(i,end));
    pts(ix).xy = [pts(ix).xy;points(i,1:2)];
end

for i=1:nz
    plot(pts(i).xy(:,1),pts(i).xy(:,2),'bo');
    pause
end
%%                  For each cluster, apply convex hull
end