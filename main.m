close all
clear
clc
%%
tic
addpath(genpath('./gpml'))
gridsize = 40;
timesize = 4;
vidObj = VideoWriter('simulation_1.avi');
vidObj.FrameRate =  1;
open(vidObj);

% --- Data pre-processing
max_ = [0 0 0];
min_ = [1000 1000 1000];
cutoff = 500;
X = [];
for i=1:timesize
    file_name = sprintf('JJ%d',i);
    load(file_name)
    %data{i} = JJ;
    max_temp = max(JJ);
    min_temp = min(JJ);
    if (max_(1)<max_temp(1))
        max_(1) = max_temp(1);
    end
    if (max_(2)<max_temp(2))
        max_(2) = max_temp(2);
    end
    if (max_(3)<max_temp(3))
        max_(3) = max_temp(3);
    end
    if (min_(1)>min_temp(1))
        min_(1) = min_temp(1);
    end
    if (min_(2)>min_temp(2))
        min_(2) = min_temp(2);
    end
    if (min_(3)>min_temp(3))
        min_(3) = min_temp(3);
    end
    index = randperm(size(JJ,1));
    X_temp = JJ(index(1:cutoff),:);
    time_temp = i*ones(size(X_temp,1),1);
    %X = [X;[X_temp time_temp]];
    X = [X;[time_temp X_temp]];
    X_true{i} = [time_temp X_temp];
end
clear JJ

nt = size(X,1);
y  = zeros(nt,1);

%% Config GPML
option.x_max = max_(1);
option.y_max = max_(2);
option.z_max = max_(3);

option.x_min = min_(1);
option.y_min = min_(2);
option.z_min = min_(3);

ts = 1; % sampling time
option.x_mesh = linspace(option.x_min,option.x_max,gridsize);
option.y_mesh = linspace(option.y_min,option.y_max,gridsize);
option.z_mesh = linspace(option.z_min,option.z_max,gridsize);
t_grid = 1:ts:ts*timesize;

[S1,S2,S3] = meshgrid(option.x_mesh,option.y_mesh,option.z_mesh); % [middle shortest longest]
S_temp = [S1(:),S2(:),S3(:)];
S = zeros(gridsize^3*timesize,3);
% for i=1:timesize
%     S((i-1)*gridsize^2+1:i*gridsize^2,1) = t_grid(i);
%     for j = 1:gridsize
%         S((i-1)*gridsize^3+(j-1)*gridsize+1:...
%             (i-1)*gridsize^3+j*gridsize,2:4)...
%             = [x_grid(j)*ones(gridsize,1) y_grid];
%     end
% end

for i=1:timesize
    S((i-1)*gridsize^3+1:i*gridsize^3,1) = t_grid(i);
    S((i-1)*gridsize^3+1:i*gridsize^3,2:4) = S_temp;
end

covfunc  = @covSEard; 
likfunc  = @likGauss;
meanfunc = @meanOne;

hyp.cov(1) = log(1);   % bandwidth of time
hyp.cov(2) = log(5);   % bandwidth of x
hyp.cov(3) = log(5);   % bandwidth of y
hyp.cov(4) = log(10);  % bandwidth of z

hyp.cov(5) = log(1);   % \sig_f
hyp.lik = log(0.03);

%% Run GPML
[est, var] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, y, S);
field_est = reshape(est,[gridsize gridsize gridsize timesize]);
%field_est = reshape(est,[],timesize);

% --- Visualize the result
% --- Initialize video writer
axis tight
set(gca,'nextplot','replacechildren');
h = figure(1);
set(h,'position',[50 50 1600 900]);

cm = colormap(pink);
sur_thres = 0.002;
for i=1:timesize
    S_ = [];
    field_temp = reshape(field_est(:,:,:,i),[],1);
    %field_temp = reshape(field_est(:,i),[],1);
    for j=1:size(field_temp,1)
        if (field_temp(j)>=0&&field_temp(j)<=sur_thres)
            S_ = [S_;S((i-1)*gridsize^3+j,:)];
        end
    end
    
    cid = repmat(cm(1,:),size(S_,1),1);
    markerSize = 25*ones(size(S_,1),1);
    cid_true = repmat(cm(1,:),size(X_true{i},1),1);
    markerSize_true = 25*ones(size(X_true{i},1),1);
    
    subplot(2,3,1)
    scatter3(S_(:,2),S_(:,3),S_(:,4),markerSize,cid);
    view([1,0,0])
    title('X view: est')
    
    subplot(2,3,4)
    scatter3(X_true{i}(:,2),X_true{i}(:,3),X_true{i}(:,4),markerSize_true,cid_true);
    view([1,0,0])
    title('X view: true')
    
    subplot(2,3,2)
    scatter3(S_(:,2),S_(:,3),S_(:,4),markerSize,cid);
    view([0,1,0])
    title('Y view: est')
    
    subplot(2,3,5)
    scatter3(X_true{i}(:,2),X_true{i}(:,3),X_true{i}(:,4),markerSize_true,cid_true);
    view([0,1,0])
    title('Y view: true')
    
    subplot(2,3,3)
    scatter3(S_(:,2),S_(:,3),S_(:,4),markerSize,cid);
    view([0,0,1])
    title('Z view: est')
    
    subplot(2,3,6)
    scatter3(X_true{i}(:,2),X_true{i}(:,3),X_true{i}(:,4),markerSize_true,cid_true);
    view([0,0,1])
    title('Z view:true')
    
    suptitle(sprintf('t=%d',i));
    
    pause(1)
    %file_name = sprintf('./fig_%d.jpg',i);
    %saveas(h,file_name);
    currFrame = getframe(figure(1));
    writeVideo(vidObj,currFrame);
end

close(vidObj);
ellapsedTime = toc;
fprintf('Ellapsed time: %d (mins)\n',ellapsedTime/60);