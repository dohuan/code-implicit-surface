close all
clear
clc
%%
tic
addpath(genpath('./gpml'))

Pat_list(1).name = 'H';
Pat_list(1).numScan = 7;
Pat_list(2).name = 'J';
Pat_list(2).numScan = 5;
Pat_list(3).name = 'K';
Pat_list(3).numScan = 5;

gridsize = 40;
pts_mode = 1; % 0: only on surface, 1: inner line included, 2: outter surface included
cutoff = 1000;
edgeLimit = 0.004;

for p=1:size(Pat_list,2)
    fprintf('Processing patient %s...\n',Pat_list(p).name);
    max_ = [0 0 0];
    min_ = [1000 1000 1000];
    cutoff = 1000;
    X = [];
    y = [];
    timesize = Pat_list(p).numScan-1; % last scan for prediction
    for i=1:timesize
        %file_name = sprintf('JJ%d_inner',i);
        file_name = ['./Patient_Data/HPCC_data/'...
            Pat_list(p).name Pat_list(p).name num2str(i) '_inner'];
        load(file_name)
        %data{i} = JJ;
        
        max_temp = max(data.on_surface);
        min_temp = min(data.on_surface);
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
        if (pts_mode==0)
            index = randperm(size(data.on_surface,1));
            X_temp = data.on_surface(index(1:cutoff),1:3);
            time_temp = i*ones(size(X_temp,1),1);
            X = [X;[time_temp X_temp]];
        else
            % --- scale radii of inner data to be \in [-1,0]
            min_inner = min(data.inner_line(:,end));
            max_inner = max(data.inner_line(:,end));
            inner_temp = (data.inner_line(:,end)-max_inner)/(max_inner-min_inner);
            inner_temp = inner_temp*edgeLimit;
            %inner_temp = data.inner_line(:,end);
            
            index = randperm(size(data.on_surface,1));
            X_temp = data.on_surface(index(1:cutoff),1:3);
            y_temp = data.on_surface(index(1:cutoff),end);
            X_temp = [X_temp;data.inner_line(:,1:3)];
            y_temp = [y_temp;inner_temp];
            y = [y;y_temp];
            time_temp = i*ones(size(X_temp,1),1);
            X = [X;[time_temp X_temp]];
        end
    end
    
    clear data
    %nt = size(X,1);
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
    
    for i=1:timesize
        S((i-1)*gridsize^3+1:i*gridsize^3,1) = t_grid(i);
        S((i-1)*gridsize^3+1:i*gridsize^3,2:4) = S_temp;
    end
    
    covfunc  = @covSEard;
    likfunc  = @likGauss;
    
    %meanfunc = @meanOne;
    %meanfunc = {'meanProd',{'meanZero','meanOne'}};
    
    meanfunc = @meanConst;
    hyp.mean = edgeLimit;
    
    hyp.cov(1) = log(1);   % bandwidth of time
    hyp.cov(2) = log(5);   % bandwidth of x
    hyp.cov(3) = log(5);   % bandwidth of y
    hyp.cov(4) = log(10);  % bandwidth of z
    
    hyp.cov(5) = log(1);   % \sig_f
    hyp.lik = log(0.03);
    
    % --- Create X vector for prediction at t = last scan
    S_test = [Pat_list(p).numScan*ones(size(S_temp,1),1) S_temp];
    [est, var] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, y, S_test);
    ellapsedTime = toc;

    fprintf('Ellapsed time: %d (mins)\n',ellapsedTime/60);
    save_name = [Pat_list(p).name Pat_list(p).name '_est'];
    save(save_name,'est','var');
    
end