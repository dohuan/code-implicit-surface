function out = patient_process(pat_info,option)

fprintf('Progressing patient %s...\n',pat_info.name);

max_ = [0 0 0];
min_ = [1000 1000 1000];

time_convert = 1; % convert month to day
X = [];
y = [];

%timesize = pat_info.numScan-1; % last scan for prediction
data_path = './Patient_Data/HPCC_data/';
for i=1:pat_info.numScan
    
    file_name = [data_path pat_info.name num2str(i) '_inner'];
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
        time_temp = time_convert*data.time_stamp*ones(size(X_temp,1),1);
        % --- Update X and y
        X = [X;[time_temp X_temp]];
        y = [y;y_temp];
    else
        % --- scale radii of inner data to be \in [-1,0]
        min_inner = min(data.inner_line(:,end));
        max_inner = max(data.inner_line(:,end));
        inner_temp = (data.inner_line(:,end)-max_inner)/(max_inner-min_inner);
        %inner_temp = inner_temp*option.edgeLimit;
        %inner_temp = data.inner_line(:,end);
        
        index = randperm(size(data.on_surface,1));
        X_temp = data.on_surface(index(1:option.cutoff),1:3);
        y_temp = data.on_surface(index(1:option.cutoff),end);
        X_temp = [X_temp;data.inner_line(:,1:3)];
        y_temp = [y_temp;inner_temp];
        time_temp = time_convert*data.time_stamp*ones(size(X_temp,1),1);
        
        % --- Update X and y
        X = [X;[time_temp X_temp]];
        y = [y;y_temp];
    end
end

clear data

%%                          Standardize data
std_info(1).mean = mean(X(:,1));     % t
std_info(1).std = sqrt(var(X(:,1)));
std_info(2).mean = mean(X(:,2));     % x
std_info(2).std = sqrt(var(X(:,2)));
std_info(3).mean = mean(X(:,3));     % y
std_info(3).std = sqrt(var(X(:,3)));
std_info(4).mean = mean(X(:,4));     % z
std_info(4).std = sqrt(var(X(:,4)));
% --- Standardize data
for i=1:size(X,2)
    X(:,i) = (X(:,i)-std_info(i).mean)./std_info(i).std;
end

time_line = unique(X(:,1)); % after standardized

X_test = X(1:(pat_info.numScan-1)*option.cutoff,:);
X_train = X(1:(pat_info.numScan-2)*option.cutoff,:);
                                
y_test  = y(1:(pat_info.numScan-1)*option.cutoff,:);
y_train = y(1:(pat_info.numScan-2)*option.cutoff,:);

S_validate = X((pat_info.numScan-2)*option.cutoff+1:...
                                    (pat_info.numScan-1)*option.cutoff,2:4);
S_true = X((pat_info.numScan-1)*option.cutoff+1:...
                                    (pat_info.numScan)*option.cutoff,2:4);

%%                          Run implicit surface GPML

% --- Config spatial
x_max = (max_(1)-std_info(2).mean)/std_info(2).std;
y_max = (max_(2)-std_info(3).mean)/std_info(3).std;
z_max = (max_(3)-std_info(4).mean)/std_info(4).std;

x_min = (min_(1)-std_info(2).mean)/std_info(2).std;
y_min = (min_(2)-std_info(3).mean)/std_info(3).std;
z_min = (min_(3)-std_info(4).mean)/std_info(4).std;

x_mesh = linspace(x_min,x_max,option.gridsize);
y_mesh = linspace(y_min,y_max,option.gridsize);
z_mesh = linspace(z_min,z_max,option.gridsize);
[S1,S2,S3] = meshgrid(x_mesh,y_mesh,z_mesh); % [middle shortest longest]
S_temp = [S1(:),S2(:),S3(:)];

covfunc  = @covSEard;
likfunc  = @likGauss;

meanfunc = @meanConst;
hyp.mean = option.edgeLimit;

hyp.cov(1) = log((pat_info.band_t-std_info(1).mean)/std_info(1).std);   % bandwidth of time
hyp.cov(2) = log(0.25);   % bandwidth of x
hyp.cov(3) = log(0.25);   % bandwidth of y
hyp.cov(4) = log(0.05);  % bandwidth of z

hyp.cov(5) = log(1);   % \sig_f
hyp.lik = log(0.03);

%% --- Find optimal hyper-parameters from initial guess

hyp = minimize(hyp, @gp, -10, @infExact, meanfunc, covfunc, likfunc, X_train, y_train);

% exp(hyp.cov)
% exp(hyp.mean)
% exp(hyp.lik)

%%                  Do greedy search for threshold value
% --- Create spatio-temporal grid for latest scan in the training

Grid_train = [time_line(pat_info.numScan-1)*ones(size(S_temp,1),1) S_temp];
[est_train, ~] = gp(hyp, @infExact,meanfunc,covfunc,likfunc,X_train,y_train,Grid_train);

[thres_train_min,Haus_min,Haus_track] = thresCal(pat_info.name,est_train,...
                                                 S_validate,S_temp,option,1,0);

out.Hause_min_train = Haus_min;
out.thres_train = thres_train_min;

%%                      Predict the test scan
% --- Create spatio-temporal grid for prediction at t = last scan
fprintf('\nPredicting ...\n');

Grid_test = [time_line(pat_info.numScan)*ones(size(S_temp,1),1) S_temp];
[est_test,~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X_test, y_test, Grid_test);

S_test_est = [];
for j=1:size(est_test,1)
    %if (est_test(j,1)>=(option.thres_min-option.thres_step)...
    if (est_test(j,1)>=0&&est_test(j,1)<=thres_train_min) 
        S_test_est = [S_test_est;S_temp(j,:)];
    end
end

if(isempty(S_test_est)==1)
    fprintf('Prediction is empty!\n');
end

out.est_train = est_train;
out.est_test = est_test;
%%                  Convert 3-D model to unstandardized coordinates
for i=1:3
    S_test_est(:,i) = S_test_est(:,i).*std_info(i+1).std + std_info(i+1).mean;
    S_true(:,i) = S_true(:,i).*std_info(i+1).std + std_info(i+1).mean;
end
out.S_true = S_true;
out.S_est = S_test_est;
out.Haus_dist = HausdorffDist(S_test_est,S_true);

