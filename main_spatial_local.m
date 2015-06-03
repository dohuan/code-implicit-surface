close all
clear
clc
%%
addpath(genpath('./gpml'))
load JJ1_inner
tic

pts_mode = 1; % 0: only on surface, 1: inner line included, 2: outter surface included

max_ = max(data.on_surface);
option.x_max = max_(1);
option.y_max = max_(2);
option.z_max = max_(3);
min_ = min(data.on_surface);
option.x_min = min_(1);
option.y_min = min_(2);
option.z_min = min_(3);

gridsize = 50;
option.x_mesh = linspace(option.x_min,option.x_max,gridsize);
option.y_mesh = linspace(option.y_min,option.y_max,gridsize);
option.z_mesh = linspace(option.z_min,option.z_max,gridsize);
[S1,S2,S3] = meshgrid(option.x_mesh,option.y_mesh,option.z_mesh);
S = [S1(:),S2(:),S3(:)];

covfunc  = @covSEard; 
likfunc  = @likGauss;
meanfunc = @meanOne;

hyp.cov(1) = log(5);   % bandwidth of x
hyp.cov(2) = log(5);   % bandwidth of y
hyp.cov(3) = log(10);  % bandwidth of z
hyp.cov(4) = log(1);
hyp.lik = log(0.03);

if (pts_mode==0)
    % --- RANDOMLY downsampling data
    cutoff = 3000;
    index = randperm(size(data.on_surface,1));

    X = data.on_surface(index(1:cutoff),1:3);
    y = data.on_surface(index(1:cutoff),end);
    nt = size(X,1);
elseif(pts_mode==1)
    % --- scale radii of inner data to be \in [0,1]
    min_inner = min(data.inner_line(:,end));
    max_inner = max(data.inner_line(:,end));
    inner_temp = (data.inner_line(:,end)-min_inner)/(max_inner-min_inner);
    
    cutoff = 3000;
    index = randperm(size(data.on_surface,1));

    X = data.on_surface(index(1:cutoff),1:3);
    y = data.on_surface(index(1:cutoff),end);
    X = [X;data.inner_line(:,1:3)];
    %y = [y;data.inner_line(:,end)];
    y = [y;inner_temp];
    nt = size(X,1);
elseif(pts_mode==2)
    
end

%[est, var] = gp(hyp, @infExact, [], covfunc, likfunc, X, y, S);
[est, var] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, y, S);
ellapsedTime = toc;
fprintf('Ellapsed time: %d (seconds)\n',ellapsedTime);

%%                      Comment if LOCAL use
%save('run_spatial_0309_hpcc1')
%%

%%                      Comment if HPCC use
% --- Plot 3D of predicted and true model (show surface points ONLY)
cm = colormap(pink);
max_est = max(est);
sur_thres = 0.002;
S_ = [];
S_index = [];
for i=1:size(est,1)
	if (est(i)>=0&&est(i)<=sur_thres)
		S_ = [S_;S(i,:)];
        S_index = [S_index,i];
	end
end
cid = repmat(cm(1,:),size(S_,1),1);
markerSize = 25*ones(size(S_,1),1);
scatter3(S_(:,1),S_(:,2),S_(:,3),markerSize,cid);

cid = repmat(cm(1,:),size(X,1),1);
markerSize = 25*ones(size(X,1),1);
figure(2)
scatter3(X(:,1),X(:,2),X(:,3),markerSize,cid);
%%
