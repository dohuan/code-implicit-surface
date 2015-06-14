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
z_slice = 2000;
a = 1;
index_slice = find(predict(a).S_est(:,3)==predict(a).S_est(z_slice,3));
S_plot = predict(a).S_est(index_slice,:);
scatter3(S_plot(:,1),S_plot(:,2),S_plot(:,3));

z_slice = 2000;
a = 1;
index_slice = find(S_est_min_(:,3)==S_est_min_(z_slice,3));
S_plot = S_est_min_(index_slice,:);
scatter3(S_plot(:,1),S_plot(:,2),S_plot(:,3));


plot(est_test,'r')
hold on
plot([0 7e4],[thres_min thres_min],'k-')
plot(est_train)
plot([0 7e4],[opt_thres opt_thres],'r-')

plot(predict.est_test,'r')
hold on
plot([0 7e4],[predict.thres_test predict.thres_test],'r-')
plot(predict.est_train)
plot([0 7e4],[predict.thres_train predict.thres_train],'b-')