function out = patient_process_speed_1(pat_info,option)

fprintf('Progressing patient %s...\n',pat_info.name);

max_ = [0 0 0];
min_ = [1000 1000 1000];

time_convert = 30; % convert a month to days

%data_path = './Patient_Data/truncated_data_LAST-3/';
data_path = './Patient_Data/truncated_data/';
out.name = pat_info.name;
out.method = 2;
GPIS(pat_info.numScan) = struct('X',[],'y',[],'IS',[],'scantime',[]);
X_all = [];         % --- Used for standardizing data
time_all = [];      % --- Used for standardizing data

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
    time_all = [time_all;GPIS(i).scantime];
    
end

spatial_offset = 0;
max_=max_+spatial_offset;
min_=min_-spatial_offset;

time_info.max = max(time_all);
time_info.min = min(time_all);
clear data
if (option.ifStand==1)
    std_info(1).mean = mean(X_all(:,1));
    std_info(1).std = sqrt(var(X_all(:,1)));
    std_info(2).mean = mean(X_all(:,2));
    std_info(2).std = sqrt(var(X_all(:,2)));
    std_info(3).mean = mean(X_all(:,3));
    std_info(3).std = sqrt(var(X_all(:,3)));
    %std_info(4).mean = mean(time_all);
    %std_info(4).std = sqrt(var(time_all));
    
    % --- Standardize spatial data and transform time
    for i=1:pat_info.numScan
        for k=1:3
            GPIS(i).X(:,k) = (GPIS(i).X(:,k)-std_info(k).mean)./std_info(k).std;
        end
        GPIS(i).scantime = (GPIS(i).scantime-min(time_all))/(max(time_all)-min(time_all));
    end
end

% --- Config spatial grid
if (option.ifStand==1)
    x_max = (max_(1)-std_info(1).mean)/std_info(1).std;
    y_max = (max_(2)-std_info(2).mean)/std_info(2).std;
    z_max = (max_(3)-std_info(3).mean)/std_info(3).std;
    
    x_min = (min_(1)-std_info(1).mean)/std_info(1).std;
    y_min = (min_(2)-std_info(2).mean)/std_info(2).std;
    z_min = (min_(3)-std_info(3).mean)/std_info(3).std;
    
    X_all(:,1) = (X_all(:,1)-std_info(1).mean)./std_info(1).std;
    X_all(:,2) = (X_all(:,2)-std_info(2).mean)./std_info(2).std;
    X_all(:,3) = (X_all(:,3)-std_info(3).mean)./std_info(3).std;
else
    x_max = max_(1);
    y_max = max_(2);
    z_max = max_(3);
    
    x_min = min_(1);
    y_min = min_(2);
    z_min = min_(3);
end

x_mesh = linspace(x_min,x_max,option.xygridsize);
y_mesh = linspace(y_min,y_max,option.xygridsize);
z_mesh = linspace(z_min,z_max,option.zgridsize);

[S1,S2,S3] = meshgrid(x_mesh,y_mesh,z_mesh); % [middle shortest longest]
spatial_grid_pool = [S1(:),S2(:),S3(:)];

for s = 1:option.grid_size
    for i=1:length(GPIS)
        GPIS(i).IS = [];
    end
    fprintf(['Run grid: ' num2str(s) '\n']);
    spatial_grid = spatial_grid_pool(s:option.grid_size:end,:);
    % --- Non-uniform grid ---
    % [S1,S2] = meshgrid(x_mesh,y_mesh);
    % xy_mesh = [S1(:) S2(:)];
    % z_window = 0.4*(z_max-z_min)/option.zgridsize;
    % spatial_grid = [];
    % for i=1:length(z_mesh)
    % 	ix = X_all(:,3)<(z_mesh(i)+z_window) & (z_mesh(i)-z_window)<X_all(:,3);
    % 	xy = X_all(ix,[1 2]);
    %
    % 	xy_label = zeros(size(xy_mesh,1),1);
    %
    % 	for j=1:size(xy,1)
    % 		[~,tmp] = sort(sum((repmat(xy(j,:),size(xy_mesh,1),1)-xy_mesh).^2,2));
    % 		tmp = sort(tmp(1:4));
    % 		xy_label(tmp(1)) = xy_label(tmp(1)) | 1;
    % 		xy_label(tmp(2)) = xy_label(tmp(2)) | 1;
    % 		xy_label(tmp(3)) = xy_label(tmp(3)) | 1;
    % 		xy_label(tmp(4)) = xy_label(tmp(4)) | 1;
    % 	end
    %
    % 	labeltmp = reshape(xy_label,option.xygridsize,option.xygridsize);
    % 	se = strel('disk',4);
    % 	labeltmp = imclose(labeltmp,se);
    % 	xy_label = labeltmp(:);
    %
    % 	for j=1:size(xy_mesh,1)
    % 		if (xy_label(j)==1)
    % 			%plot(xy_mesh(j,1),xy_mesh(j,2),'r.');
    %             spatial_grid = [spatial_grid;[xy_mesh(j,:) z_mesh(i)]];
    % 		end
    % 	end
    % end
    
    
    % ========== Compute spatio-temporal up to T-1 as observations ========
    
    for t=1:pat_info.numScan-1
        GPIS(1:end-1) = compute_GPIS_temporal(GPIS(t).scantime,GPIS(1:pat_info.numScan-1), ...
            spatial_grid, pat_info, std_info,time_info,option);
        thres(t,1) = thresCal(pat_info.name,GPIS(t).IS.est,GPIS(t).X,...
            spatial_grid,option,0);
        thres_time(t,1) = GPIS(t).scantime;
    end
    
    thres_value = GP_thres(thres_time,thres,GPIS(end).scantime);

    %out.thres = thres;
    
    % ============================= Run EM ================================
    fprintf('Starting EM...\n');
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
    
    % --- Prediction as predict step of KF:
    delta_t_test = GPIS(end).scantime-GPIS(end-1).scantime;
    est = mean_T(:,end) + delta_t_test*A;
    sig = cov_T(:,:,end) + delta_t_test^2*SW;
    
    
    %%
    % % ===== Compute GPIS ======
    %
    % for t=2:pat_info.numScan
    %     GPIS(t) = compute_GPIS_temporal(GPIS(t).scantime,GPIS(1:t), ...
    %         spatial_grid, pat_info, std_info,time_info,option);
    % end
    % % =========================
    %
    % % --- Compute GR (growth rate) field: f'(t)=(GPIS(t+1)-GPIS(t))/\delta_t
    % GR(pat_info.numScan-1) = struct('est',[],'var',[]);
    % X = [];
    % y = [];
    % for i=1:pat_info.numScan-1
    %     GR(i).est = (GPIS(i+1).IS.est-GPIS(i).IS.est)./(GPIS(i+1).scantime-GPIS(i).scantime); % E(\hat{A})
    %     GR(i).var = (GPIS(i+1).IS.var+GPIS(i).IS.var)./(GPIS(i+1).scantime-GPIS(i).scantime)^2; % var(\hat{A}):= \Sig_W
    %     X_temp = [GPIS(i).scantime*ones(size(spatial_grid,1),1) spatial_grid];
    %     X = [X;X_temp];
    %     y = [y;GR(i).est];
    % end
    %
    % % --- Compute the threshold value
    % [thres_temp,~,~] = thresCal(pat_info.name,GPIS(end-1).IS.est,GPIS(end-1).X,...
    %     spatial_grid,option,0);
    %
    % %out = GPIS(end);
    % out.thres = thres_temp.up;
    %
    % %  !!! HERE: have to estimate GR(end) with spatial-temporal GP
    %
    % % GR  : 1,...,L-1
    % % GPIS: 1,...,L
    % out.est = GPIS(end-1).IS.est + ...
    %     (GPIS(end).scantime-GPIS(end-1).scantime)*GR(end-1).est;
    % out.var = GPIS(end-1).IS.var + ...
    %     (GPIS(end).scantime-GPIS(end-1).scantime)^2*GR(end-1).var;
    %%
    
    
    S_est = field_to_surface(thres_value,est,spatial_grid);
    %S_est = surface_refiner(S_est);
    S_true = GPIS(end).X;
    
    % --- Create CB ONLY for the LAST prediction
    temp = chol(sig);
    var_revised = temp'*temp;
    %CB = mvnrnd(out.est,out.var,option.CB_run);
    CB_field = mvnrnd(est,var_revised,option.CB_run);
    for i=1:option.CB_run
        CB_temp = field_to_surface(thres_value,CB_field(i,:)',spatial_grid);
        %CB{i} = surface_refiner(CB_temp);
        CB{i} = CB_temp;
        
    end
    
    if (option.ifStand == 1)
        %%                  Convert 3-D model to unstandardized coordinates
        for i=1:3
            S_est(:,i) = S_est(:,i).*std_info(i).std + std_info(i).mean;
            S_true(:,i) = S_true(:,i).*std_info(i).std + std_info(i).mean;
            for j=1:option.CB_run
                CB{j}(:,i) = CB{j}(:,i).*std_info(i).std + std_info(i).mean;
            end
        end
        %     out.band_t = exp(hyp.cov(1))*std_info(1).std;
        %     out.band_x = exp(hyp.cov(2))*std_info(2).std;
        %     out.band_y = exp(hyp.cov(3))*std_info(3).std;
        %     out.band_z = exp(hyp.cov(4))*std_info(4).std;
        %     out.band_f = exp(hyp.cov(5));
    else
        %     out.band_t = exp(hyp.cov(1));
        %     out.band_x = exp(hyp.cov(2));
        %     out.band_y = exp(hyp.cov(3));
        %     out.band_z = exp(hyp.cov(4));
        %     out.band_f = exp(hyp.cov(5));
    end
    
    if (s==1)
        out.CB = CB;
        out.S_est = S_est;
    else
        for i=1:option.CB_run
            out.CB{i} = [out.CB{i};CB{i}];
        end
        out.S_est = [out.S_est;S_est];
    end
    clear CB;
    
end
out.S_est = surface_refiner(out.S_est);
for i=1:option.CB_run
    out.CB{i} = surface_refiner(out.CB{i});
end
out.thres = thres_value;
out.S_true = S_true;
[out.Haus,~] = HausdorffDist(out.S_est,S_true);

end


function S = field_to_surface(thres,field,grid)
S = [];
for i=1:size(field,1)
    %if (field(i,1)>=0&&field(i,1)<=thres)
    if (field(i,1)<=thres)
        S = [S;grid(i,:)];
    end
end
end

function out = GP_thres(x,y,x_test)
covfunc = @covSEard;
likfunc = @likGauss;
meanfunc = @meanConst;

hyp.cov(1) = log(5);
hyp.cov(2) = log(100);
hyp.mean = mean(y);
hyp.lik = log(2);
%hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);
[out, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, x_test);
end





