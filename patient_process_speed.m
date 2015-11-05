function out = patient_process_speed(pat_info,option)

fprintf('Progressing patient %s...\n',pat_info.name);

max_ = [0 0 0];
min_ = [1000 1000 1000];

time_convert = 30; % convert a month to days

data_path = './Patient_Data/truncated_data/';
out.name = pat_info.name;
GPIS(pat_info.numScan) = struct('X',[],'y',[],'IS',[],'scantime',[]);
X_all = []; % --- Used for standardizing data
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
    X_temp = data.on_surface(:,1:3);
    GPIS(i).y = data.on_surface(:,end);
    GPIS(i).scantime = time_convert*data.time_stamp;
    %time_temp = time_convert*data.time_stamp*ones(size(X_temp,1),1);
    %GPIS(i).X = [time_temp X_temp];
    GPIS(i).X = X_temp;
    X_all = [X_all; GPIS(i).X];
end

clear data
if (option.ifStand==1)
    std_info(1).mean = mean(X_all(:,1));     % t
    std_info(1).std = sqrt(var(X_all(:,1)));
    std_info(2).mean = mean(X_all(:,2));     % x
    std_info(2).std = sqrt(var(X_all(:,2)));
    std_info(3).mean = mean(X_all(:,3));     % y
    std_info(3).std = sqrt(var(X_all(:,3)));
    std_info(4).mean = mean(X_all(:,4));     % z
    std_info(4).std = sqrt(var(X_all(:,4)));
    
    % --- Standardize data
    for k=1:pat_info.numScan
        GPIS(k).X = (GPIS(k).X(:,i)-std_info(i).mean)./std_info(i).std;
    end
end

% --- Config spatial
if (option.ifStand==1)
    x_max = (max_(1)-std_info(2).mean)/std_info(2).std;
    y_max = (max_(2)-std_info(3).mean)/std_info(3).std;
    z_max = (max_(3)-std_info(4).mean)/std_info(4).std;
    
    x_min = (min_(1)-std_info(2).mean)/std_info(2).std;
    y_min = (min_(2)-std_info(3).mean)/std_info(3).std;
    z_min = (min_(3)-std_info(4).mean)/std_info(4).std;
else
    x_max = max_(1);
    y_max = max_(2);
    z_max = max_(3);
    
    x_min = min_(1);
    y_min = min_(2);
    z_min = min_(3);
end

x_mesh = linspace(x_min,x_max,option.gridsize);
y_mesh = linspace(y_min,y_max,option.gridsize);
z_mesh = linspace(z_min,z_max,option.gridsize);
[S1,S2,S3] = meshgrid(x_mesh,y_mesh,z_mesh); % [middle shortest longest]
spatial_grid = [S1(:),S2(:),S3(:)];

% --- Compute GPIS
for i=1:pat_info.numScan
    GPIS(i) = compute_GPIS_spatial(GPIS(i), spatial_grid, pat_info, option);
end

% --- Compute GR (growth rate) field: f'(t)=(GPIS(t+1)-GPIS(t))/\delta_t
GR(pat_info.numScan-1) = struct('est',[],'var',[]);
for i=1:pat_info.numScan-2
    GR(i).est = (GPIS(i+1).IS.est-GPIS(i).IS.est)./(GPIS(i+1).time-GPIS(i).time);
    GR(i).var = (GPIS(i+1).IS.var-GPIS(i).IS.var)./(GPIS(i+1).time-GPIS(i).time);
end

end






function data = compute_GPIS_spatial(data,spatial_grid,pat_info,option)
covfunc  = {'covSum', {'covSEard','covNoise'}};
likfunc  = @likGauss;
meanfunc = @meanOne;
if (option.ifStand == 1)
    %hyp.cov(1) = log(abs((pat_info.band_t)/std_info(1).std));   % bandwidth of time
    hyp.cov(1) = log(abs((pat_info.band_x)/std_info(2).std));   % bandwidth of x
    hyp.cov(2) = log(abs((pat_info.band_y)/std_info(3).std));   % bandwidth of y
    hyp.cov(3) = log(abs((pat_info.band_z)/std_info(4).std));  % bandwidth of z
else
    %hyp.cov(1) = log(pat_info.band_t);
    hyp.cov(1) = log(pat_info.band_x);
    hyp.cov(2) = log(pat_info.band_y);
    hyp.cov(3) = log(pat_info.band_z);
end
hyp.cov(4) = log(pat_info.band_f);   % \sig_f
hyp.cov(5) = log(0.05);              % \sig_w: measurement noise .05
hyp.lik = log(0.05);                 % initial prior probability for hyper p(\theta) 0.03

hyp = minimize(hyp, @gp, -10, @infExact, meanfunc, covfunc, ...
                                         likfunc, data.X, data.y);
[est, var] = ...
       gaussian_process(hyp, covfunc,meanfunc, data.X,data.y,spatial_grid);
data.IS.est = est;
data.IS.var = var;
data.IS.hyp = hyp;

end

function data = compute_GPIS_temporal()

end