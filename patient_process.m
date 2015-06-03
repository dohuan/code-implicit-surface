function out = patient_process(pat_info,option)

fprintf('Progressing patient %s...\n',pat_info.name);

max_ = [0 0 0];
min_ = [1000 1000 1000];
%thres_offset = 10;
%time_ratio = exp((pat_info.numScan-1)/(pat_info.numScan-2))/...
%                              exp((pat_info.numScan)/(pat_info.numScan-1));
time_ratio = 1;
X = [];
y = [];
%timesize = pat_info.numScan-1; % last scan for prediction
timesize = pat_info.numScan-2; % last scan for prediction
for i=1:timesize
    
    file_name = ['./Patient_Data/HPCC_data/'...
        pat_info.name num2str(i) '_inner'];
    load(file_name)
    
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
    if (option.pts_mode==0)
        index = randperm(size(data.on_surface,1));
        X_temp = data.on_surface(index(1:option.cutoff),1:3);
        y_temp = data.on_surface(index(1:option.cutoff),end);
        %time_temp = i*ones(size(X_temp,1),1);
        time_temp = data.time_stamp*ones(size(X_temp,1),1);
        X = [X;[time_temp X_temp]];
        y = [y;y_temp];
    else
        % --- scale radii of inner data to be \in [-1,0]
        min_inner = min(data.inner_line(:,end));
        max_inner = max(data.inner_line(:,end));
        inner_temp = (data.inner_line(:,end)-max_inner)/(max_inner-min_inner);
        inner_temp = inner_temp*option.edgeLimit;
        %inner_temp = data.inner_line(:,end);
        
        index = randperm(size(data.on_surface,1));
        X_temp = data.on_surface(index(1:option.cutoff),1:3);
        y_temp = data.on_surface(index(1:option.cutoff),end);
        X_temp = [X_temp;data.inner_line(:,1:3)];
        y_temp = [y_temp;inner_temp];
        y = [y;y_temp];
        %time_temp = i*ones(size(X_temp,1),1);
        time_temp = data.time_stamp*ones(size(X_temp,1),1);
        X = [X;[time_temp X_temp]];
    end
    if (i==timesize)
        S_offset = data.on_surface(index(1:option.cutoff),1:3);
    end
end

%--- Load the validation scan to the train scan for prediction
file_name = ['./Patient_Data/HPCC_data/'...
    pat_info.name num2str(timesize+1) '_inner'];
load(file_name)
if (option.pts_mode==0)
    index = randperm(size(data.on_surface,1));
    X_temp = data.on_surface(index(1:option.cutoff),1:3);
    y_temp = data.on_surface(index(1:option.cutoff),end);
    %time_temp = (timesize+1)*ones(size(X_temp,1),1);
    time_temp = data.time_stamp*ones(size(X_temp,1),1);
    
    X_test = [X;[time_temp X_temp]];
    y_test = [y;y_temp];
else
    inner_temp = (data.inner_line(:,end)-max_inner)/(max_inner-min_inner);
    inner_temp = inner_temp*option.edgeLimit;
    index = randperm(size(data.on_surface,1));
    X_temp = data.on_surface(index(1:option.cutoff),1:3);
    y_temp = data.on_surface(index(1:option.cutoff),end);
    X_temp = [X_temp;data.inner_line(:,1:3)];
    y_temp = [y_temp;inner_temp];
    
    y_test = [y;y_temp];
    %time_temp = (timesize+1)*ones(size(X_temp,1),1);
    time_temp = data.time_stamp*ones(size(X_temp,1),1);
    X_test = [X;[time_temp X_temp]];
end

clear data

%%                          Run implicit surface GPML

% --- Config temporal
%ts = 1; % sampling time
%t_grid = 1:ts:ts*timesize;
% --- Config spatial
x_max = max_(1);
y_max = max_(2);
z_max = max_(3);

x_min = min_(1);
y_min = min_(2);
z_min = min_(3);

x_mesh = linspace(x_min,x_max,option.gridsize);
y_mesh = linspace(y_min,y_max,option.gridsize);
z_mesh = linspace(z_min,z_max,option.gridsize);
[S1,S2,S3] = meshgrid(x_mesh,y_mesh,z_mesh); % [middle shortest longest]
S_temp = [S1(:),S2(:),S3(:)];

% S = zeros(option.gridsize^3*timesize,3);
% for j=1:timesize
%     S((j-1)*option.gridsize^3+1:j*option.gridsize^3,1) = t_grid(j);
%     S((j-1)*option.gridsize^3+1:j*option.gridsize^3,2:4) = S_temp;
% end

covfunc  = @covSEard;
likfunc  = @likGauss;

meanfunc = @meanConst;
hyp.mean = option.edgeLimit;

%hyp.cov(1) = log(80);   % bandwidth of time
hyp.cov(1) = log(pat_info.band_t);   % bandwidth of time
hyp.cov(2) = log(5);   % bandwidth of x
hyp.cov(3) = log(5);   % bandwidth of y
hyp.cov(4) = log(10);  % bandwidth of z

hyp.cov(5) = log(1);   % \sig_f
hyp.lik = log(0.03);

%% --- Find optimal hyper-parameters from initial guess
hyp = minimize(hyp, @gp, -8, @infExact, meanfunc, covfunc, likfunc, X_test, y_test);
% exp(hyp.cov)
% exp(hyp.mean)
% exp(hyp.lik)

%%                  Do greedy search for threshold value
% --- Create spatio-temporal grid for latest scan in the training
file_name = ['./Patient_Data/HPCC_data/'...
    pat_info.name num2str(timesize+1)...
    '_inner'];
load(file_name);

%Grid_train = [(pat_info.numScan-1)*ones(size(S_temp,1),1) S_temp];
Grid_train = [data.time_stamp*ones(size(S_temp,1),1) S_temp];
[est_train, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, y, Grid_train);

index = randperm(size(data.on_surface,1));
S_validate = data.on_surface(index(1:option.cutoff),1:3);

[thres_train_min,Haus_min,~] = thresCal(pat_info.name,est_train,...
                                                 S_validate,S_temp,option,1,0);

out.Hause_min_train = Haus_min;
out.thres_train = thres_train_min;


%%                  Calculate the threshold offset
% Grid_offset = [timesize*ones(size(S_temp,1),1) S_temp];
% [est_offset, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, y, Grid_offset);
% 
% [thres_offset_min,~,~] = thresCal(pat_info.name,est_offset,...
%                                                    S_offset,S_temp,option);
% 
% thres_offset = abs(thres_offset_min-thres_train_min);

%% Load the true scan

file_name = ['./Patient_Data/HPCC_data/'...
    pat_info.name num2str(pat_info.numScan)...
    '_inner'];
load(file_name);
%%                      Predict the test scan
% --- Create spatio-temporal grid for prediction at t = last scan
fprintf('\nPredicting ...\n');
%Grid_test = [pat_info.numScan*ones(size(S_temp,1),1) S_temp];
Grid_test = [data.time_stamp*ones(size(S_temp,1),1) S_temp];
[est_test,~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X_test, y_test, Grid_test);

S_test_est = [];
for j=1:size(est_test,1)
%     if (est_test(j,1)>=(option.thres_min-option.thres_step)...
%             &&est_test(j,1)<=thres_train_min)
    if (est_test(j,1)>=(option.thres_min-option.thres_step)...
        &&est_test(j,1)<=thres_train_min*time_ratio) 
        S_test_est = [S_test_est;S_temp(j,:)];
    end
end

if(isempty(S_test_est)==1)
    fprintf('Prediction is empty!\n');
end

out.S_est = S_test_est;
out.est_train = est_train;
out.est_test = est_test;

index = randperm(size(data.on_surface,1));
out.S_true = data.on_surface(index(1:option.cutoff),1:3);
out.Haus_dist = HausdorffDist(S_test_est,out.S_true);

