function exportMesh( pointCloud, filePath )
%%                    Author: Huan Do  - dohuan@msu.edu
% pointCloud: Nx3 matrix of N coordinates of 3-D points

fid = fopen(filePath,'w');
fprintf(fid,'ply\n');
fprintf(fid,'format ascii 1.0\n');
fprintf(fid,'element vertex %i\n',size(pointCloud,1));
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');
fprintf(fid,'end_header\n');
nume = size(pointCloud,1);
for i=1:nume
    fprintf('progressing...%d %%\n',round(i/nume*100));
    fprintf(fid,'%f %f %f\n',pointCloud(i,1),pointCloud(i,2),pointCloud(i,3));
end
fclose(fid);
fprintf('File %s Was Successfully Written\n', filePath);
end