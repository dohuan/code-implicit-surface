% --- plot 3D of predicted whole model
cm = colormap(pink);
max_est = max(est);
colorSize = size(cm,1);
sur_thres = 0.007;
for i=1:size(est,1)
	if (est(i)==0)
		cid(i,:) = cm(1,:);
		vis{i} = 'on';
	elseif (est(i)>0)
		if (est(i)<=sur_thres)
			cid(i,:) = cm(ceil(est(i)/sur_thres*64),:);
			vis{i} = 'on';
		else
			cid(i,:) = cm(end,:);
			vis{i} = 'off';
		end
	else
		cid(i,:) = cm(end,:);
		vis{i} = 'off';
	end
end
markerSize = 25*ones(size(est,1),1);
scatter3(S(:,1),S(:,2),S(:,3),markerSize,cid,'.');

% --- Plot 3D of predicted and true model (show surface points ONLY)

cm = colormap(pink);
max_est = max(est);
sur_thres = 0.0001;
S_ = [];
for i=1:size(est,1)
	if (est(i)>=0&&est(i)<=sur_thres)
		S_ = [S_;S(i,:)];
	end
end
cid = repmat(cm(1,:),size(S_,1),1);
markerSize = 25*ones(size(S_,1),1);
scatter3(S_(:,1),S_(:,2),S_(:,3),markerSize,cid);

cid = repmat(cm(1,:),size(X,1),1);
markerSize = 25*ones(size(X,1),1);
figure(2)
scatter3(X(:,1),X(:,2),X(:,3),markerSize,cid);

% --- Test Yufei code about creating the field
Psi = [10 1 1 10]; % psi_1, psi_2x, psi_2y, sig_t
sig_w = 1e-5; % small number to avoid numerical problem
covfunc = {'covSum', {'covSEard','covNoise'}};
loghyper = [log(Psi(1));log(Psi(2));log(Psi(3));log(Psi(4));log(sig_w)];

num_grid = 21;
num_time = 30; % sec
ts = 1; % sampling time
XMIN = -2; XMAX = 2;
YMIN = -2; YMAX = 2;
ZMIN = -100; ZMAX = 100;
x_grid = linspace(XMIN,XMAX,num_grid)';
y_grid = linspace(YMIN,YMAX,num_grid)';
t_grid = 1:ts:ts*num_time;
s_grid = zeros(num_grid^2*num_time,3);

for i = 1:num_time
    s_grid((i-1)*num_grid^2+1:i*num_grid^2,1) = t_grid(i);
    for j = 1:num_grid
        s_grid((i-1)*num_grid^2+(j-1)*num_grid+1:(i-1)*num_grid^2+j*num_grid,2:3)...
            = [x_grid(j)*ones(num_grid,1) y_grid];
    end
end

Z_grid = reshape(chol(feval(covfunc{:},loghyper,s_grid))'...
    *randn(num_grid^2*num_time,1),num_grid,num_grid,num_time);
	
% --- Print file to .ply extension

% --- Read data from "Centerlines_J.xlsx" and merger on-surface and inner into .mat file
% --- J1
clear
fileDir = './Centerlines_J.xlsx';
sheetName = 'J1';
range = 'A2:D190';
[raw,dataName] = xlsread(fileDir,sheetName,range);
raw(:,end) = - raw(:,end);
load JJ1
data.on_surface = [JJ,zeros(size(JJ,1),1)];
data.inner_line = raw;
save('JJ1_inner','data');

% --- J2
clear
fileDir = './Centerlines_J.xlsx';
sheetName = 'J2';
range = 'A2:D246';
[raw,dataName] = xlsread(fileDir,sheetName,range);
raw(:,end) = - raw(:,end);
load JJ2
data.on_surface = [JJ,zeros(size(JJ,1),1)];
data.inner_line = raw;
save('JJ2_inner','data');

% --- J3
clear
fileDir = './Centerlines_J.xlsx';
sheetName = 'J3';
range = 'A2:D251';
[raw,dataName] = xlsread(fileDir,sheetName,range);
raw(:,end) = - raw(:,end);
load JJ3
data.on_surface = [JJ,zeros(size(JJ,1),1)];
data.inner_line = raw;
save('JJ3_inner','data');

% --- J4
clear
fileDir = './Centerlines_J.xlsx';
sheetName = 'J2';
range = 'A2:D252';
[raw,dataName] = xlsread(fileDir,sheetName,range);
raw(:,end) = - raw(:,end);
load JJ4
data.on_surface = [JJ,zeros(size(JJ,1),1)];
data.inner_line = raw;
save('JJ4_inner','data');


% --- Compute Hausdorff distance
[hd,~] = HausdorffDist(S_,X(1:cutoff,:));


% --- Test Hausdorff function
X =[0,0,0;1,0,0;1,1,0;0,1,0];
Y =[0,0,1;1,0,1;1,0,0;0,2,0];
scatter3(X(:,1),X(:,2),X(:,3),'r');
hold on
scatter3(Y(:,1),Y(:,2),Y(:,3),'b*');
hold off
[hd,~] = HausdorffDist(X,Y)



scatter3(S_test_est(:,1),S_test_est(:,2),S_test_est(:,3))
scatter3(S_est(:,1),S_est(:,2),S_est(:,3))

figure(3)
hold on
plot(est_train)
plot([0 7e4],[thres_train_min thres_train_min],'k-')

figure(4)
hold on
plot(est_test,'r')
plot(est_train,'b')
plot([0 7e4],[thres_train_min thres_train_min],'k-')
plot([0 7e4],[thres_train_min+thres_offset thres_train_min+thres_offset],'k--')


a = 1;
hold on
plot(predict(a).est_test,'r')
plot(predict(a).est_train,'b')
plot([0 7e4],[predict(a).thres_train predict(a).thres_train],'k-')

% --- Plot a slice of 3-D organ
z_slice = unique(predict.S_est(:,end),'rows','stable');
for i=1:size(z_slice,1)
	figure(i);
	index_slice = find(predict.S_est(:,3)==z_slice(i));
	%S_plot = predict.S_est(index_slice,:);
	%scatter3(S_plot(:,1),S_plot(:,2),S_plot(:,3));
	S_plot = predict.S_est(index_slice,1:2);
	plot(S_plot(:,1),S_plot(:,2),'bo');
end



z_slice = 2000;
a = 1;
index_slice = find(S_est_min_(:,3)==S_est_min_(z_slice,3));
S_plot = S_est_min_(index_slice,:);
scatter3(S_plot(:,1),S_plot(:,2),S_plot(:,3));


% --- show convex hull and k-mean cluster (patient BB)
z_slice = 1000; % 650: two segments
a = 1;
index_slice = find(S_test_est(:,3)==S_test_est(z_slice,3));
S_plot = S_test_est(index_slice,:);
plot(S_plot(:,1),S_plot(:,2),'bo');
S_xy = S_plot(:,1:2);
T = kmeans(S_xy,2);
group1 = [];
group2 = [];
pts_thres = 3;
ch_pts1 = [];
ch_pts2 = [];
hold on
for j=1:size(T,1)
	if (T(j)==1)
		group1 = [group1;S_xy(j,:)];
	else
		group2 = [group2;S_xy(j,:)];
	end
end
plot(group1(:,1),group1(:,2),'b*');
plot(group2(:,1),group2(:,2),'r*');
if (size(group1,1)>pts_thres&&size(group2,1)>pts_thres)
	temp1 = boundary_select(group1);
	temp2 = boundary_select(group2);
	%xy_temp = [temp1;temp2];
	ch_pts1 = [ch_pts1;temp1];
	ch_pts2 = [ch_pts2;temp2];
elseif(size(group1,1)<=pts_thres&&size(group2,1)>pts_thres)
	xy_temp = boundary_select(group2);
	ch_pts2 = [ch_pts2;xy_temp];
elseif (size(group2,1)<=pts_thres&&size(group1,1)>pts_thres)
	xy_temp = boundary_select(group1);
	ch_pts1 = [ch_pts1;xy_temp];
end
plot(ch_pts1(:,1),ch_pts1(:,2),'k--','LineWidth',1.5);
plot(ch_pts2(:,1),ch_pts2(:,2),'k--','LineWidth',1.5);
xlabel('(mm)');
ylabel('(mm)');
set(gca,'FontSize',16);



	
	
plot(est_test,'r')
hold on
plot([0 7e4],[thres_min thres_min],'k-')
plot(est_train)
plot([0 7e4],[opt_thres opt_thres],'r-')

plot(predict.est_test,'r')
hold on
plot([0 7e4],[predict.thres_test predict.thres_test],'r-','LineWidth',2)
plot(predict.est_train)
plot([0 7e4],[predict.thres_train.up predict.thres_train.up],'b-','LineWidth',2)

% test surface_refiner.m
load PP13_test
out = surface_refiner(predict.S_est);


scatter3(predict.S_est(:,1),predict.S_est(:,2),predict.S_est(:,3))
hold on
scatter3(predict.S_true(:,1),predict.S_true(:,2),predict.S_true(:,3),'ro')



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


% ---------------
plot(predict(1).est_train,'b.');
axis tight
hold on
plot([1 7e4],[predict(1).thres_train.up predict(1).thres_train.up],'k-','LineWidth',1.5);
set(gca,'FontSize',16);
xlabel('points on the grid');
ylabel('posterior mean');

% ---------------
plot(predict(1).est_train,'b.');
axis tight
hold on
plot(predict(1).est_test,'r.');
plot([1 7e4],[predict(1).thres_train.up predict(1).thres_train.up],'k-','LineWidth',1.5);
plot([1 7e4],[predict(1).thres_test predict(1).thres_test],'k--','LineWidth',1.5);
set(gca,'FontSize',16);
xlabel('points on the grid');
ylabel('posterior mean');


% ----------------
field = reshape(predict(1).est_train,option.gridsize,option.gridsize,option.gridsize);


% ---------- Create contour plot include +std and -std
z_line = unique(predict.S_est(:,end),'rows','stable');
index = [5 10 15 20 25];
slice_list = z_line(index);
hold on
rec = [70 40;170 40;170 110;70 110;70 40];
for i=1:size(slice_list,1)
	ix = find(predict.S_est(:,3)==slice_list(i));
	xy = predict.S_est(ix,1:2);
	xy = boundary_select(xy);
	z = slice_list(i)*ones(size(xy,1),1);
	plot3(xy(:,1),xy(:,2),z,'r-','LineWidth',2);
	plot3(rec(:,1),rec(:,2),slice_list(i)*ones(size(rec,1),1),'k','LineWidth',1);
	
	if (isempty(predict.S_est_up)==0)
		ix = find(predict.S_est_up(:,3)==slice_list(i));
		xy = predict.S_est_up(ix,1:2);
		xy = boundary_select(xy);
		z = slice_list(i)*ones(size(xy,1),1);
		plot3(xy(:,1),xy(:,2),z,'b-.','LineWidth',2);
	end
	
	if (isempty(predict.S_est_down)==0)
		ix = find(predict.S_est_down(:,3)==slice_list(i));
		xy = predict.S_est_down(ix,1:2);
		xy = boundary_select(xy);
		z = slice_list(i)*ones(size(xy,1),1);
		plot3(xy(:,1),xy(:,2),z,'g--','LineWidth',2);
	end
end
hold off
set(gca,'Fontsize',16);




% ---------- Create contour plot include +std and -std (BEST and WORST)
predict_best = predict(1);
set(0,'defaultfigurecolor',[1 1 1])
%subplot(2,1,1)
figure(1)
z_line = unique(predict_best.S_est(:,end),'rows','stable');
index = [5 10 15 20 25];
slice_list = z_line(index);
hold on

rec = [70 30;140 30;140 90;70 90;70 30];
for i=1:size(slice_list,1)
	ix = find(predict_best.S_est(:,3)==slice_list(i));
	xy = predict_best.S_est(ix,1:2);
	xy = boundary_select(xy);
	z = slice_list(i)*ones(size(xy,1),1);
	plot3(xy(:,1),xy(:,2),z,'r-','LineWidth',2);
	plot3(rec(:,1),rec(:,2),slice_list(i)*ones(size(rec,1),1),'k','LineWidth',1);
	
	if (isempty(predict_best.S_est_up)==0)
		ix = find(predict_best.S_est_up(:,3)==slice_list(i));
		xy = predict_best.S_est_up(ix,1:2);
		xy = boundary_select(xy);
		z = slice_list(i)*ones(size(xy,1),1);
		plot3(xy(:,1),xy(:,2),z,'b-.','LineWidth',2);
	end
	
	if (isempty(predict_best.S_est_down)==0)
		ix = find(predict_best.S_est_down(:,3)==slice_list(i));
		xy = predict_best.S_est_down(ix,1:2);
		xy = boundary_select(xy);
		z = slice_list(i)*ones(size(xy,1),1);
		plot3(xy(:,1),xy(:,2),z,'g--','LineWidth',2);
	end
end
hold off
xlabel('(mm)','FontSize',14);
ylabel('(mm)','FontSize',14);
zlabel('(mm)','FontSize',14);
set(gca,'Fontsize',16);


predict_worst = predict(2);
set(0,'defaultfigurecolor',[1 1 1])
%subplot(2,1,2)
figure(2)
z_line = unique(predict_worst.S_est(:,end),'rows','stable');
index = [5 10 15 20 25];
slice_list = z_line(index);
hold on
rec = [70 40;170 40;170 110;70 110;70 40];
for i=1:size(slice_list,1)
	ix = find(predict_worst.S_est(:,3)==slice_list(i));
	xy = predict_worst.S_est(ix,1:2);
	xy = boundary_select(xy);
	z = slice_list(i)*ones(size(xy,1),1);
	plot3(xy(:,1),xy(:,2),z,'r-','LineWidth',2);
	plot3(rec(:,1),rec(:,2),slice_list(i)*ones(size(rec,1),1),'k','LineWidth',1);
	
	if (isempty(predict_worst.S_est_up)==0)
		ix = find(predict_worst.S_est_up(:,3)==slice_list(i));
		xy = predict_worst.S_est_up(ix,1:2);
		xy = boundary_select(xy);
		z = slice_list(i)*ones(size(xy,1),1);
		plot3(xy(:,1),xy(:,2),z,'b-.','LineWidth',2);
	end
	
	if (isempty(predict_worst.S_est_down)==0)
		ix = find(predict_worst.S_est_down(:,3)==slice_list(i));
		xy = predict_worst.S_est_down(ix,1:2);
		xy = boundary_select(xy);
		z = slice_list(i)*ones(size(xy,1),1);
		plot3(xy(:,1),xy(:,2),z,'g--','LineWidth',2);
	end
end
hold off
xlabel('(mm)','FontSize',14);
ylabel('(mm)','FontSize',14);
zlabel('(mm)','FontSize',14);
set(gca,'Fontsize',16);






% --------- Print result to screen
hold on
fprintf('name band_f band_t band_x band_y band_z Haust_dist\n');
for i=1:size(predict,2)
	fprintf('%s ',predict(i).name);
	fprintf('%.1f ',predict(i).band_f);
	fprintf('%.1f ',predict(i).band_t);
	fprintf('%.1f ',predict(i).band_x);
	fprintf('%.1f ',predict(i).band_y);
	fprintf('%.1f ',predict(i).band_z);
	fprintf('%.1f ',predict(i).Haus_dist);
	fprintf('\n');
	plot(Pat_list(i).numScan,predict(i).Haus_dist,'ks','LineWidth',2);
	text(Pat_list(i).numScan+0.2,predict(i).Haus_dist+0.2,Pat_list(i).name(2:end));
	error(i) = predict(i).Haus_dist;
end
fprintf('Mean of Haus dist: %f\n',mean(error));
fprintf('Std of Haus dist: %f\n',sqrt(var(error)));
xlim([2 8]);
xlabel('Number of scans','FontSize',14);
ylabel('Hausdoff distance','FontSize',14);
box on
set(gca,'FontSize',16);



% ----------- Plot linear fit and print out R^2
Pat_list = patient_list([],0);
saveFile = './results/2015824_1142(final)/';
for i=1:size(Pat_list,2)
	M = load([saveFile Pat_list(i).name]); 
	predict(i).Haus_dist = M.predict.Haus_dist;
	clear M;
end
hold on
for i=1:size(predict,2)
	x(i) = Pat_list(i).numScan;
	y(i) = predict(i).Haus_dist;
	plot(Pat_list(i).numScan,predict(i).Haus_dist,'ks','LineWidth',2);
	text(Pat_list(i).numScan+0.2,predict(i).Haus_dist+0.2,Pat_list(i).ID);
end
p = polyfit(x,y,1);
yfit =  p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
t = 0:1:8;
yplot = p(1)*t + p(2);
std_error = sqrt(SSresid/length(y));
fprintf('R square: %.3f\n',rsq);
fprintf('Standard Error: %.3f\n',std_error);
plot(t,yplot,'k','LineWidth',2);
xlim([2 8]);
xlabel('Number of scans','FontSize',14);
ylabel('Hausdoff distance','FontSize',14);
box on
set(gca,'FontSize',16);



% --------- Plot Haus dist vs. std of scanning periods
time_data = xlsread('./Patient_Data/Dt_patients_truncated.xlsx','Sheet1');
count = 1;
Pat_list = patient_list([],0);
saveFile = './results/2015824_1142(final)/';
hold on
for i=1:size(Pat_list,2)
	M = load([saveFile Pat_list(i).name]); 
	predict(i).Haus_dist = M.predict.Haus_dist;
	clear M;
	for j=1:Pat_list(i).numScan
		var_temp(j) = time_data(count,1);
		count = count + 1;
	end
	predict(i).var_scan = sqrt(var(var_temp));
	% --- Plot
	plot(predict(i).var_scan,predict(i).Haus_dist,'ks','LineWidth',2);
	text(predict(i).var_scan+0.2,predict(i).Haus_dist+0.2,Pat_list(i).name);
	x(i) = predict(i).var_scan;
	y(i) = predict(i).Haus_dist;
end
p = polyfit(x,y,1);
yfit =  p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
t = min(x):1:max(x);
yplot = p(1)*t + p(2);
std_error = sqrt(SSresid/length(y));
fprintf('R square: %.3f\n',rsq);
fprintf('Standard Error: %.3f\n',std_error);
plot(t,yplot,'k','LineWidth',2);
xlim([2 8]);
xlabel('std of scan periods','FontSize',14);
ylabel('Hausdoff distance','FontSize',14);
box on
set(gca,'FontSize',16);
hold off


% --------- Plot Haus dist vs. std of scanning periods (with Haus dist MANUALLY inserted)

Pat_list = patient_list([],0);

a(1).Haus_dist = 17.86;
a(2).Haus_dist = 11.28;
a(3).Haus_dist = 8.32;
a(4).Haus_dist = 26.06;
a(5).Haus_dist = 9.85;
a(6).Haus_dist = 6.68;
a(7).Haus_dist = 8.32;
for i=1:size(Pat_list,2)
	%M = load([saveFile Pat_list(i).name]); 
	predict(i).Haus_dist = a(i).Haus_dist;
	clear M;
end
hold on
for i=1:size(predict,2)
	x(i) = Pat_list(i).numScan;
	y(i) = predict(i).Haus_dist;
	plot(Pat_list(i).numScan,predict(i).Haus_dist,'ks','LineWidth',2);
	text(Pat_list(i).numScan+0.2,predict(i).Haus_dist+0.2,Pat_list(i).ID);
end
p = polyfit(x,y,1);
yfit =  p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
t = 0:1:8;
yplot = p(1)*t + p(2);
std_error = sqrt(SSresid/length(y));
fprintf('R square: %.3f\n',rsq);
fprintf('Standard Error: %.3f\n',std_error);
plot(t,yplot,'k','LineWidth',2);
xlim([2 8]);
xlabel('Number of scans','FontSize',14);
ylabel('Hausdoff distance','FontSize',14);
box on
set(gca,'FontSize',16);


% --------- check symmetric problem due to numerical floating error ---------
temp = predict.var_test-predict.var_test';
count = 1;
for i=1:size(temp,1)
	for j=1:size(temp,2)
		if (temp(i,j)~=0)
			a(count) = temp(i,j);
			count = count + 1;
		end
	end
end

for i=1:size(predict.var_test,1)
	for j=1:size(predict.var_test,2)
		Sig(i,j) = round(predict.var_test(i,j)*1e5)/1e5;
	end
end


% ----------
hold on
scatter3(predict.S_true(:,1),predict.S_true(:,2),predict.S_true(:,3),'k+');
scatter3(predict.CB{1}(:,1),predict.CB{1}(:,2),predict.CB{1}(:,3),'b*');
scatter3(predict.CB{2}(:,1),predict.CB{2}(:,2),predict.CB{2}(:,3),'ro');
scatter3(predict.CB{3}(:,1),predict.CB{3}(:,2),predict.CB{3}(:,3),'gx');

% ---------- Plot true cloud along with credible band clouds
hold on
scatter3(predict.S_true(:,1),predict.S_true(:,2),predict.S_true(:,3),'k+');
for i=1:size(predict.CB,2)
	scatter3(predict.CB{i}(:,1),predict.CB{i}(:,2),predict.CB{i}(:,3),'b.');
	%alpha(h,.5);
end
hold off

% ------- Plot diagonal of the test covariance matrix
plot(diag(predict.var_test));

% --------------------------
change_name('P01');
change_name('P02');
change_name('P03');
change_name('P04');
change_name('P11');
change_name('P12');
change_name('P14');
change_name('P15');
change_name('P21');
change_name('P22');
change_name('P23');
change_name('P24');
change_name('P25');
change_name('P26');
change_name('P31');
change_name('P32');
change_name('P33');
change_name('P34');
change_name('P41');
change_name('P42');
change_name('P43');


% ---------- Check scan of patient H
for i=1:7
	figure(i);
	dix = (i-1)*1000 + 1;
	uix = i*1000;
	scatter3(X(dix:uix,2),X(dix:uix,3),X(dix:uix,4));
end

% ------------ Plot times between scans
time_data = xlsread('./Patient_Data/Dt_patients_truncated.xlsx','Sheet1');
time_data = time_data(:,5);
numScan = [7 6 5 5 4 6 4];
mark = {'bd-' 'rs-' 'gv-' 'c^-' 'mp-' 'k+-' 'yo-'};
hold on
count = 1;
temp = [];
for i=1:7
	clear temp;
	for j=1:numScan(i)
		temp(j,1) = time_data(count);
		count = count + 1;
	end
	plot(temp,mark{i},'MarkerSize',10,'LineWidth',1.5);
end
legend('H','I','J','K','P10','P12','P13');
box on
xlabel('Scan','FontSize',16);
ylabel('Scan Time (days)','FontSize',16);
set(gca,'FontSize',16);




% --------------- Read max dia and compute confidence interval
file_name = './results/max_dia_PCL_101515_onesheet.xlsx';
confidence_reg(file_name);



% -------------------
x = randn(1000,1);                      % Create Data
%SEM = std(x)/sqrt(length(x));               % Standard Error
SEM = std(x);               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;  
[h,n] = hist(x,100);
h = h./sum(h);
bar(n,h);
hold on
plot([mean(x) mean(x)],[0 max(h)],'LineWidth',2);
plot([CI(1) CI(1)],[0 max(h)],':','LineWidth',2);
plot([CI(2) CI(2)],[0 max(h)],':','LineWidth',2);
hold off
psum = 0;
for j=1:size(n,2)
	if (n(j)>=CI(1)&&n(j)<=CI(2))
		psum = psum + h(j);
	end
end
fprintf('Psum is: %.2f\n',psum);




% ----------- Use quantile (90% CI)
x = randn(1000,1);  
lowbound = quantile(x,.05);
highbound = quantile(x,.95);
CI = [lowbound highbound];
[h,n] = hist(x,100);
h = h./sum(h);
bar(n,h);
hold on
plot([mean(x) mean(x)],[0 max(h)],'LineWidth',2);
plot([CI(1) CI(1)],[0 max(h)],':','LineWidth',2);
plot([CI(2) CI(2)],[0 max(h)],':','LineWidth',2);
hold off
psum = 0;
for j=1:size(n,2)
	if (n(j)>=CI(1)&&n(j)<=CI(2))
		psum = psum + h(j);
	end
end
fprintf('Psum is: %.2f\n',psum);



% -------------- Plot to check results
if (size(Pat_list,2)==1)
	figure('name','Cross-section view');
	z_line = unique(predict.S_true(:,end),'rows','stable');
	
	for i=1:size(z_line,1)
		pts_true(i).z = z_line(i);
		pts_true(i).xy = [];
	end
	for i=1:size(predict.S_true,1)
		ix = find(z_line==predict.S_true(i,end));
		pts_true(ix).xy = [pts_true(ix).xy; predict.S_true(i,1:2)];
	end
	hold on
	for j=1:option.CB_run
		
		point_tmp = predict.CB{j};
		pts_CB = [];
		for i=1:size(z_line,1)
			pts_CB(i).z = z_line(i);
			pts_CB(i).xy = [];
		end
		for i=1:size(point_tmp,1)
			ix = find(z_line==point_tmp(i,end));
			pts_CB(ix).xy = [pts_CB(ix).xy; point_tmp(i,1:2)];
		end
		
		for k=1:length(z_line)
			subplot(4,5,k)
			plot(pts_CB(k).xy(:,1),pts_CB(k).xy(:,2),'r.','MarkerSize',2);
		end
	end
	for k=1:length(z_line)
		subplot(4,5,k)
		title(sprintf('z=%d',z_line(k)))
		plot(pts_true(k).xy(:,1),pts_true(k).xy(:,2),'r.','MarkerSize',2);
	end
	hold off
end


% === Poly fit of threshold values ===
for i=1:7
	subplot(2,4,i)
	t=  [];
	for j=1:Pat_list(i).numScan-1
		load(['./Patient_Data/truncated_data/' Pat_list(i).name num2str(j) '_inner']);
		t = [t;data.time_stamp];
	end
	if (Pat_list(i).numScan-1<4)
		p = polyfit(t,thres(i).a,2);
	else
		p = polyfit(t,thres(i).a,3);
	end
	x = min(t):1:max(t);
	y = polyval(p,x);
	load(['./Patient_Data/truncated_data/' Pat_list(i).name num2str(Pat_list(i).numScan) '_inner']);
	x_test = data.time_stamp;
	y_test = polyval(p,x_test);
	
	hold on;
	plot(t,thres(i).a,'bo-');
	plot(x,y,'r-');
	plot(x_test,y_test,'sk');
	hold off
	title(Pat_list(i).name);
	box on;
end

% === GP of threshold values ===
for i=1:7
	subplot(2,4,i)
	x=  [];
	for j=1:Pat_list(i).numScan-1
		load(['./Patient_Data/truncated_data/' Pat_list(i).name num2str(j) '_inner']);
		x = [x;data.time_stamp];
	end
	y = thres(i).a;
	
	load(['./Patient_Data/truncated_data/' Pat_list(i).name num2str(Pat_list(i).numScan) '_inner']);
	x_test = (min(x):1:data.time_stamp)';
	
	covfunc = @covSEard; 
	likfunc = @likGauss;
	meanfunc = @meanConst;
	
	hyp.cov(1) = log(5);
	hyp.cov(2) = log(100);
	hyp.mean = mean(y);
	hyp.lik = log(2);
	%hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);
	[y_est, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, x_test);
	
	hold on;
	plot(x,y,'bo-');
	plot(x_test,y_est,'r-');
	plot(x_test(end),y_est(end),'sk');
	hold off
	title(Pat_list(i).name);
	box on;
end

% ----- Test non-uniform grid
theta = linspace(0,2*pi,1000);
x = cos(theta);
y = sin(theta);
N = 100;
bar(1/N*hist(x,100));
figure(2)
bar(1/N*hist(y,100));



% -------- Non-uniform grid with independent x and y-axis
x_mesh = linspace(x_min,x_max,option.gridsize);
y_mesh = linspace(y_min,y_max,option.gridsize);

z_mesh = linspace(z_min,z_max,option.gridsize);
z_window = 0.4*(z_max-z_min)/option.gridsize;
%step = .2*(x_max-y_min)/option.gridsize;
histlength = 100;
iter = 300;
for i=1:length(z_mesh)
	ix = find(X_all(:,3)<(z_mesh(i)+z_window) & (z_mesh(i)-z_window)<X_all(:,3));
	xy = X_all(ix,[1 2]);
	
	%figure(2)
	subplot(4,5,i);
	hold on
	plot(xy(:,1),xy(:,2),'bo');
	
	[xh,xhp] = hist(xy(:,1),histlength);
	[yh,yhp] = hist(xy(:,2),histlength);
	xh = 1/histlength*xh;
	yh = 1/histlength*yh;
	
	%x_mesh_ = x_mesh;
	%y_mesh_ = y_mesh;
	
	x_mesh_ = linspace(min(xy(:,1)),max(xy(:,1)),option.gridsize);
	y_mesh_ = linspace(min(xy(:,2)),max(xy(:,2)),option.gridsize);
	
	
	%figure(1)
	%subplot(4,5,i);
	%hold on
	%plot(xy(:,1),xy(:,2),'bo');
	%[S1,S2] = meshgrid(x_mesh_,y_mesh_);
	%plot(S1,S2,'r.');
	
	step = .3*(max(xy(:,1))-min(xy(:,1)))/option.gridsize;
	
	for k=1:iter
		for j=2:length(x_mesh_)-1
			[~,tmp] = sort(abs(xhp-x_mesh_(j)));
			tmp = sort(tmp(1:2));
			movestep = step*sin(atan((xh(tmp(2))-xh(tmp(1)))/(xhp(tmp(2))-xhp(tmp(1)))));
			if (x_mesh_(j)+movestep>x_mesh_(1) & x_mesh_(j)+movestep<x_mesh_(end))
				x_mesh_(j) = x_mesh_(j) + movestep;
			end
			%x_mesh_(j) = x_mesh_(j) + movestep;
			
			
			[~,tmp] = sort(abs(yhp-y_mesh_(j)));
			tmp = sort(tmp(1:2));
			movestep = step*sin(atan((xh(tmp(2))-xh(tmp(1)))/(xhp(tmp(2))-xhp(tmp(1)))));
			if (y_mesh_(j)+movestep>y_mesh_(1) & y_mesh_(j)+movestep<y_mesh_(end))
				y_mesh_(j) = y_mesh_(j) + movestep;
			end
			%y_mesh_(j) = y_mesh_(j) + movestep;
		end
	end
	x_mesh_ = sort(x_mesh_);
	y_mesh_ = sort(y_mesh_);
	
	x_mesh_z(:,i) = x_mesh_;
	y_mesh_z(:,i) = y_mesh_;
	
	%figure(2)
	subplot(4,5,i);
	[S1,S2] = meshgrid(x_mesh_,y_mesh_);
	plot(S1,S2,'r.');
	
	hold off
end




% -------- Non-uniform grid with connected components
density = 20;
x_mesh = linspace(x_min,x_max,density);
y_mesh = linspace(y_min,y_max,density);
[S1,S2] = meshgrid(x_mesh,y_mesh);
xy_mesh = [S1(:) S2(:)];

z_mesh = linspace(z_min,z_max,option.gridsize);
z_window = 0.4*(z_max-z_min)/option.gridsize;
meshsum = 0;
for i=1:length(z_mesh)
	ix = find(X_all(:,3)<(z_mesh(i)+z_window) & (z_mesh(i)-z_window)<X_all(:,3));
	xy = X_all(ix,[1 2]);
	
	xy_label = zeros(size(xy_mesh,1),1);
	
	for j=1:size(xy,1)
		[~,tmp] = sort(sum((repmat(xy(j,:),size(xy_mesh,1),1)-xy_mesh).^2,2));
		tmp = sort(tmp(1:4));
		xy_label(tmp(1)) = xy_label(tmp(1)) | 1;
		xy_label(tmp(2)) = xy_label(tmp(2)) | 1;
		xy_label(tmp(3)) = xy_label(tmp(3)) | 1;
		xy_label(tmp(4)) = xy_label(tmp(4)) | 1;
	end
	
	labeltmp = reshape(xy_label,density,density);
	se = strel('disk',4);
	labeltmp = imclose(labeltmp,se);
	xy_label = labeltmp(:);
	
	subplot(4,5,i);
	hold on
	plot(xy(:,1),xy(:,2),'bo');
	for j=1:size(xy_mesh,1)
		if (xy_label(j)==1)
			plot(xy_mesh(j,1),xy_mesh(j,2),'r.');
		end
	end
	hold off
	meshsum = meshsum + sum(xy_label);
end



% ------- Run EM with different initial values of A and Sig_w

figure(1)
hold on

count = 1;
A0   = (GPIS(2).IS.est-GPIS(1).IS.est)/(GPIS(2).scantime-GPIS(1).scantime);
SW0  = GPIS(1).IS.var;

mu0  = GPIS(1).IS.est;
Sig0 = GPIS(1).IS.var;
logL = zeros(option.EM_run,1);

while (count<option.EM_run+1)
    if (count ==1)
        [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
        KF_update(A0,SW0,mu0,Sig0,GPIS,pat_info,spatial_grid);
    else
        [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
        KF_update(A,SW,mu0,Sig0,GPIS,pat_info,spatial_grid);
    end
    [A,SW, logL(count)] = EM_update(mean_T,cov_T,cov_T_,mean_T0,cov_T0,...
    cov_T_10,Sig0,GPIS,pat_info);
    count = count + 1;
end

plot(logL,'b','LineWidth',2);

% -----------
count = 1;
A0   = (GPIS(3).IS.est-GPIS(2).IS.est)/(GPIS(3).scantime-GPIS(2).scantime);
SW0  = GPIS(2).IS.var;

mu0  = GPIS(1).IS.est;
Sig0 = GPIS(1).IS.var;
logL = zeros(option.EM_run,1);

while (count<option.EM_run+1)
    if (count ==1)
        [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
        KF_update(A0,SW0,mu0,Sig0,GPIS,pat_info,spatial_grid);
    else
        [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
        KF_update(A,SW,mu0,Sig0,GPIS,pat_info,spatial_grid);
    end
    [A,SW, logL(count)] = EM_update(mean_T,cov_T,cov_T_,mean_T0,cov_T0,...
    cov_T_10,Sig0,GPIS,pat_info);
    count = count + 1;
end

plot(logL,'r--','LineWidth',2);

% -----------
count = 1;
A0   = 3*(GPIS(2).IS.est-GPIS(1).IS.est)/(GPIS(2).scantime-GPIS(1).scantime);
SW0  = 3*GPIS(1).IS.var;

mu0  = GPIS(1).IS.est;
Sig0 = GPIS(1).IS.var;
logL = zeros(option.EM_run,1);

while (count<option.EM_run+1)
    if (count ==1)
        [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
        KF_update(A0,SW0,mu0,Sig0,GPIS,pat_info,spatial_grid);
    else
        [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
        KF_update(A,SW,mu0,Sig0,GPIS,pat_info,spatial_grid);
    end
    [A,SW, logL(count)] = EM_update(mean_T,cov_T,cov_T_,mean_T0,cov_T0,...
    cov_T_10,Sig0,GPIS,pat_info);
    count = count + 1;
end

plot(logL,'g-.','LineWidth',2);

box on
xlabel('iterations','FontSize',16);
ylabel('Log likelihood','FontSize',16);
set(gca,'FontSize',16);






% ------- Run EM with different initial values of A and Sig_w
% ONLY USE FOR PATIENT H
plotstyle{1} = 'bs-';
plotstyle{2} = 'ro--';
plotstyle{3} = 'bd-';
plotstyle{4} = 'g+-';
plotstyle{5} = 'r*:';

for i =1:5
	count = 1;
	A0   = (GPIS(i+1).IS.est-GPIS(i).IS.est)/(GPIS(i+1).scantime-GPIS(i).scantime);
	SW0  = GPIS(i).IS.var;

	mu0  = GPIS(1).IS.est;
	Sig0 = GPIS(1).IS.var;
	logL(:,i) = zeros(option.EM_run,1);
	
	A_track(:,count,i) = A0;
	Sigw_track(:,:,count,i) = SW0;
	normA_track = [];
	normSigw_track = [];
	while (count<option.EM_run+1)
		if (count ==1)
			[mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
			KF_update(A0,SW0,mu0,Sig0,GPIS,pat_info,spatial_grid);
		else
			[mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,~] =...
			KF_update(A,SW,mu0,Sig0,GPIS,pat_info,spatial_grid);
		end
		[A,SW, logL(count,i)] = EM_update(mean_T,cov_T,cov_T_,mean_T0,cov_T0,...
		cov_T_10,Sig0,GPIS,pat_info);
		count = count + 1;
		
		A_track(:,count,i) = A;
		normA_track = [normA_track;norm(A_track(:,count,i)-A_track(:,count-1,i),'fro')];
		%normA_track = [normA_track;mean(A_track(:,count))];
		
		Sigw_track(:,:,count,i) = SW;
		normSigw_track = [normSigw_track;norm(Sigw_track(:,:,count,i)-Sigw_track(:,:,count-1,i),'fro')];
		%normSigw_track = [normSigw_track;mean(mean(Sigw_track(:,:,count)))];
	end
	%figure(i)
	%plot(logL(:,i),plotstyle{i},'LineWidth',2);
	figure(1)
	hold on
	plot(normA_track,plotstyle{i},'LineWidth',2);
	
	figure(2)
	hold on
	plot(normSigw_track,plotstyle{i},'LineWidth',2);
end


figure(1)
box on
xlabel('iterations','FontSize',16);
ylabel('Norms of the increments of A','FontSize',16);
set(gca,'FontSize',16);

figure(2)
box on
xlabel('iterations','FontSize',16);
ylabel('Norms of the increments of \Sigma_w','FontSize',16);
set(gca,'FontSize',16);




% --- (cont of above code)
close all
A_mean = mean(A_track,3);
A_mean = A_mean(:,end);
Sigw_mean = mean(Sigw_track,4);
Sigw_mean = Sigw_mean(:,:,end);
for i=1:count
	for j=1:5
		Atemp = A_track;
		Atemp(:,:,j) = [];
		
		Sigwtemp = Sigw_track;
		Sigwtemp(:,:,:,j) = [];
		
		for k=1:4
			tmpA(k,1) = norm(A_track(:,i,j)-Atemp(:,i,k));
			tmpSigw(k,1) = norm(Sigw_track(:,:,i,j)-Sigwtemp(:,:,i,k),'fro');
		end
		distA(j,i) = max(tmpA);
		distSigw(j,i) = max(tmpSigw);
		
		distAmean(j,i) = norm(A_track(:,i,j)-A_mean);
		distSigwmean(j,i) = norm(Sigw_track(:,:,i,j)-Sigw_mean,'fro');
	end
end


for j=1:5
	figure(1)
	hold on
	plot(distA(j,:),plotstyle{j},'LineWidth',2);
	
	figure(2)
	hold on
	plot(distSigw(j,:),plotstyle{j},'LineWidth',2);
	
	figure(3)
	hold on
	plot(distAmean(j,:),plotstyle{j},'LineWidth',2);
	
	figure(4)
	hold on
	plot(distSigwmean(j,:),plotstyle{j},'LineWidth',2);
	
end

figure(1)
box on
xlabel('iterations','FontSize',16);
ylabel('Evolution of differences of A','FontSize',16);
set(gca,'FontSize',16);

figure(2)
box on
xlabel('iterations','FontSize',16);
ylabel('Evolution of differences of \Sigma_w','FontSize',16);
set(gca,'FontSize',16);

figure(3)
box on
xlabel('iterations','FontSize',16);
ylabel('Evolution of differences of A to A_\infty','FontSize',16);
set(gca,'FontSize',16);
axis tight

figure(4)
box on
xlabel('iterations','FontSize',16);
ylabel('Evolution of differences of \Sigma_w to \Sigma_{w\infty}','FontSize',16);
set(gca,'FontSize',16);
axis tight


% --- (cont of above code)---interpolate A(x)

%A_ = (100-5)*(A-min(A))./(max(A)-min(A))+5;

V = reshape(A,5,5,5);
%X = reshape(spatial_grid(:,1),5,5,5);
%Y = reshape(spatial_grid(:,1),5,5,5);
%Z = reshape(spatial_grid(:,1),5,5,5);

iscale = 5;

gridtmp_x = linspace(min(spatial_grid(:,1)),max(spatial_grid(:,1)),iscale);
gridtmp_y = linspace(min(spatial_grid(:,2)),max(spatial_grid(:,2)),iscale);
gridtmp_z = linspace(min(spatial_grid(:,3)),max(spatial_grid(:,3)),iscale);

[X,Y,Z] = meshgrid(gridtmp_x,gridtmp_y,gridtmp_z);

iscale = 50;

gridtmp_x = linspace(min(spatial_grid(:,1)),max(spatial_grid(:,1)),iscale);
gridtmp_y = linspace(min(spatial_grid(:,2)),max(spatial_grid(:,2)),iscale);
gridtmp_z = linspace(min(spatial_grid(:,3)),max(spatial_grid(:,3)),iscale);

[X_,Y_,Z_] = meshgrid(gridtmp_x,gridtmp_y,gridtmp_z);

V_ = interp3(X,Y,Z,V,X_,Y_,Z_);
V_ = V_(:);
C = colormap('lines');
for i=1:numel(V_)
	ix = round((64-1)*(V_(i)-min(V_))./(max(V_)-min(V_))+1); % color index
	V_ci(i,:) = C(ix,:);
end

z_attitude = 0; zoffset = .07;
Xtmp = X_all(5001:6000,:);
ix_slice = find(Xtmp(:,3)>z_attitude-zoffset&Xtmp(:,3)<z_attitude+zoffset);
X_pre = Xtmp(ix_slice,1:2);
K = convhull(X_pre);
X_pre = X_pre(K,:);

Xtmp = X_all(6001:7000,:);
ix_slice = find(Xtmp(:,3)>z_attitude-zoffset&Xtmp(:,3)<z_attitude+zoffset);
X_cur = Xtmp(ix_slice,1:2);
K = convhull(X_cur);
X_cur = X_cur(K,:);





h=scatter3(X_(:),Y_(:),Z_(:),10*ones(size(V_ci,1),1),V_ci,'filled','MarkerFaceAlpha',.2);
%h.FaceAlpha = .2;

xslice = [];%[-2,0,2]; 
yslice = [];%[-2,0,2]; 
zslice = [-2,-1,0];%[-2,0,2]; 

h=slice(X_,Y_,Z_,V_,xslice,yslice,zslice);
shading interp
hold on
scatter3(X_all(5001:6000,1),X_all(5001:6000,2),X_all(5001:6000,3),'ks');
scatter3(X_all(6001:7000,1),X_all(6001:7000,2),X_all(6001:7000,3),'ko');
plot3(X_pre(:,1),X_pre(:,2),0.01*ones(size(X_pre,1),1),'w-','LineWidth',2);
plot3(X_cur(:,1),X_cur(:,2),0.01*ones(size(X_cur,1),1),'w--','LineWidth',2);




% ---------------- Scalar-version likelihood function
f = @(x,y) log(x)+1/x*(1+2*y+y^2)
y = linspace(-10,10,100);
x = linspace(0.01,0.1,100);
[S1,S2]=meshgrid(x,y);
value = [S1(:) S2(:)];
for i =1:size(value,1)
	z(i,1) = f(value(i,1),value(i,2));
end
z_ = reshape(z,100,100);
 surf(z_,'EdgeColor','none');
 xlabel('\Sigma_w')
 ylabel('A')
 
 
 
 
 % ---------- Convert all patients' last scans to PLY
patlist = patient_list_speed([],0);
for i=1:length(patlist) 
	loadfilename = ['./Patient_Data/truncated_data/' patlist(i).name num2str(patlist(i).numScan) '_inner' ];
	savefilename = ['./Patient_Data/PLY_true_truncated/' patlist(i).name '_true.ply'];
	load(loadfilename);
	exportMesh(data.on_surface,savefilename);
	clear data;
end



% ---------- Fix data error with LAST-3 setup
% --- Subtract all scan dates with the last
data_path = './Patient_Data/truncated_data_last3/';
patList = patient_list_speed([],0,'last3');
for i=1:length(patList)
	fprintf('Fixing patient: %s\n',patList(i).name);
	for j=1:patList(i).numScan
		file_name = [data_path patList(i).name num2str(j) '_inner'];
		load(file_name)
		if (j==1)
			timelast = data.time_stamp;
		end
		data.time_stamp = data.time_stamp-timelast;
		save(file_name,'data');
	end
end


% --- Plot slice of 3D structure
z_slice = unique(S_est(:,end),'rows','stable');
for i=1:size(z_slice,1)
	figure(i);
	index_slice = find(S_est(:,3)==z_slice(i));
	%S_plot = predict.S_est(index_slice,:);
	%scatter3(S_plot(:,1),S_plot(:,2),S_plot(:,3));
	S_plot = S_est(index_slice,1:2);
	plot(S_plot(:,1),S_plot(:,2),'bo');
end




% -------------------------
figure(1)
hold on
C = colormap('lines');
for i=1:5
	loadFile = ['./Patient_Data/truncated_data/KK' num2str(i) '_inner'];
	load(loadFile);
	C_temp = repmat(C(i,:), [size(data.on_surface,1) 1]);
	subplot(1,5,i)
	h(i) = scatter3(data.on_surface(:,1), data.on_surface(:,2), data.on_surface(:,3),[],C_temp, 'filled');
	%h(i) = scatter3(data.on_surface(:,1), data.on_surface(:,2), data.on_surface(:,3), 'filled');
	view([0 0])
	h_legend{i} = ['scan ' num2str(i)];
	xlabel(['scan ' num2str(i)]);
	zlim([140 340]);
	box on
end
%legend(h,'scan 1','scan 2','scan 4','scan 5','scan 6','scan 7');
%legend(h_legend);
%[~, I] = unique(C);
%p = findobj(gca,'Type','Patch');
%legend(p(I),h_legend);
hold off
=======

% --- save GPIS and std_info
saveFile = './results/GPIS-saved/P3';
for i=1:length(GPIS)-1
	GPIS_truncated(i).IS.hyp = GPIS(i).IS.hyp;
end
save(saveFile,'GPIS_truncated','std_info','time_info','pat_info');






% --- Compute hyper parameter
file = './results/hyperParameter.txt';
fileID = fopen(file,'w');
fprintf(fileID,'ID,sig_t,sig_x,sig_y,sig_z \n');
data(1).ID = 'H';
data(1).num = 6;
data(2).ID = 'I';
data(2).num = 5;
data(3).ID = 'J';
data(3).num = 4;
data(4).ID = 'K';
data(4).num = 4;
data(5).ID = 'P0';
data(5).num = 3;
data(6).ID = 'P2';
data(6).num = 5;
data(7).ID = 'P3';
data(7).num = 3;

for i=1:length(data)
	loadFile = ['./results/GPIS-saved/' data(i).ID];
	load(loadFile);
	hyptrack = [];
	for j=1:length(GPIS_truncated)
		hyptmp = exp(GPIS_truncated(j).IS.hyp.cov(1:4));
		hyptmp(1) = hyptmp(1)*(time_info.max-time_info.min)+time_info.min;
		hyptmp(2) = hyptmp(2)*std_info(1).std;
		hyptmp(3) = hyptmp(3)*std_info(2).std;
		hyptmp(4) = hyptmp(4)*std_info(3).std;
		
		hyptrack = [hyptrack; hyptmp];
	end
	hyp = mean(hyptrack);
	fprintf(fileID,'%s, %.2f, %.2f, %.2f, %.2f \n',data(i).ID, hyp(1), hyp(2), hyp(3), hyp(4));
end
fclose(fileID);





% --- Plot P2 posterior mean

hold on
plot(GPIS(5).IS.est,'+')
plot(est,'.')
plot([0 length(est)],[thres_value thres_value],'k--','LineWidth',2)
plot([0 length(est)],[thres(end) thres(end)],'k-','LineWidth',2)
hold off
axis tight
box on
xlabel('Points');
ylabel('Potential field');
set(gca,'FontSize',16);



% --- Plot figure of GP field evolves with time
load P2_test_postMean





























