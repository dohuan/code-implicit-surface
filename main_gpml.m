close all
clear
clc
%%
addpath(genpath('./gpml'))
load JJ1
tic

max_ = max(JJ);
option.x_max = max_(1);
option.y_max = max_(2);
option.z_max = max_(3);
min_ = min(JJ);
option.x_min = min_(1);
option.y_min = min_(2);
option.z_min = min_(3);

gridsize = 40;
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


% --- RANDOMLY downsampling data
cutoff = 500;
index = randperm(size(JJ,1));

X = JJ(index(1:cutoff),:);
nt = size(X,1);
y = zeros(nt,1);

%[est, var] = gp(hyp, @infExact, [], covfunc, likfunc, X, y, S);
[est, var] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, y, S);
ellapsedTime = toc;
fprintf('Ellapsed time: %d (seconds)\n',ellapsedTime);

% --- Plot 3D of predicted and true model (show surface points ONLY)

cm = colormap(pink);
max_est = max(est);
sur_thres = 0.001;
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

