% Write result to text file
fid = fopen('./cloud_point.txt','w');
for i=1:size(predict.S_est,1)
	 fprintf(fid,[num2str(predict.S_est(i,1)) ',' ...
				  num2str(predict.S_est(i,2)) ',' ...
				  num2str(predict.S_est(i,3))]);
     fprintf(fid,'\n');
end
fclose(fid);

% Plot the mesh
load run_all_061815
h = pointCloud2mesh(predict.S_est,[0 0 1],0.2);
coord = h.vertices;
conn = h.triangles;


hold on;
for i=1:size(conn,1)
	plot3(coord(conn(i,[1 2]),1),coord(conn(i,[1 2]),2),coord(conn(i,[1 2]),3),'b*-');
	plot3(coord(conn(i,[2 3]),1),coord(conn(i,[2 3]),2),coord(conn(i,[2 3]),3),'b*-');
	plot3(coord(conn(i,[1 3]),1),coord(conn(i,[1 3]),2),coord(conn(i,[1 3]),3),'b*-');
end
scatter3(predict.S_est(:,1),predict.S_est(:,2),predict.S_est(:,3),'ro');
hold off

% -----------------------
t = 0:pi/10:2*pi;
[x,y,z] = cylinder(sin(pi/10*t));
coord = [reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
h = pointCloud2mesh(coord,[0 0 1],0.6);
conn = h.triangles;

hold on;
for i=1:size(conn,1)
	plot3(coord(conn(i,[1 2]),1),coord(conn(i,[1 2]),2),coord(conn(i,[1 2]),3),'b*-');
	plot3(coord(conn(i,[2 3]),1),coord(conn(i,[2 3]),2),coord(conn(i,[2 3]),3),'b*-');
	plot3(coord(conn(i,[1 3]),1),coord(conn(i,[1 3]),2),coord(conn(i,[1 3]),3),'b*-');
end
hold off


% --------------
load ./results/HH_HPCC_est051115

[S1,S2,S3] = meshgrid(option.x_mesh,option.y_mesh,option.z_mesh); % [middle shortest longest]
S_temp = [S1(:),S2(:),S3(:)];
S_test = [];
thres_min = 1.6;
for j=1:size(est,1)
	if (est(j,1)>=0&&est(j,1)<=thres_min)
		S_test = [S_test;S_temp(j,:)];
	end
	fprintf('Process:%0.2f%%\n',j/size(est,1)*100);
end
S_est = surface_refiner(S_est);


% --------------
load ./Patient_Data/AA1

temp = AA(:,1);
AA(:,1) = AA(:,3);
AA(:,3) = temp;

h = pointCloud2mesh(AA,[0 0 1],0.2);
coord = h.vertices;
conn = h.triangles;

hold on;
for i=1:size(conn,1)
	plot3(coord(conn(i,[1 2]),1),coord(conn(i,[1 2]),2),coord(conn(i,[1 2]),3),'b*-');
	plot3(coord(conn(i,[2 3]),1),coord(conn(i,[2 3]),2),coord(conn(i,[2 3]),3),'b*-');
	plot3(coord(conn(i,[1 3]),1),coord(conn(i,[1 3]),2),coord(conn(i,[1 3]),3),'b*-');
end
hold off

% -------------
load ./Patient_Data/AA1
[the,rho,z] = cart2pol(AA(:,1),AA(:,2),AA(:,3));
scatter3(the,rho,z);
