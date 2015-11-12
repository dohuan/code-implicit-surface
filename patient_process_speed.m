function out = patient_process_speed(pat_info,option)

fprintf('Progressing patient %s...\n',pat_info.name);

max_ = [0 0 0];
min_ = [1000 1000 1000];

time_convert = 30; % convert a month to days

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
    GPIS(i) = compute_GPIS_spatial(GPIS(i), spatial_grid, pat_info, std_info,option);
end

% --- Compute GR (growth rate) field: f'(t)=(GPIS(t+1)-GPIS(t))/\delta_t
GR(pat_info.numScan-1) = struct('est',[],'var',[]);
X = [];
y = [];
for i=1:pat_info.numScan-1
    GR(i).est = (GPIS(i+1).IS.est-GPIS(i).IS.est)./(GPIS(i+1).scantime-GPIS(i).scantime); % E(\hat{A})
    GR(i).var = (GPIS(i+1).IS.var+GPIS(i).IS.var)./(GPIS(i+1).scantime-GPIS(i).scantime)^2; % var(\hat{A}):= \Sig_W
    X_temp = [GPIS(i).scantime*ones(size(spatial_grid,1),1) spatial_grid];
    X = [X;X_temp];
    y = [y;GR(i).est];
end

% --- Compute the threshold value
[thres_temp,~,~] = thresCal(pat_info.name,GPIS(end-1).IS.est,GPIS(end-1).X,...
                                                    spatial_grid,option,0);

%out = GPIS(end);
out.thres = thres_temp.up;

%  !!! HERE: have to estimate GR(end) with spatial-temporal GP

% GR  : 1,...,L-1
% GPIS: 1,...,L
out.est = GPIS(end-1).IS.est + ...
                   (GPIS(end).scantime-GPIS(end-1).scantime)*GR(end-1).est;
out.var = GPIS(end-1).IS.var + ...
                 (GPIS(end).scantime-GPIS(end-1).scantime)^2*GR(end-1).var;

out.S_est = field_to_surface(out.thres,out.est,spatial_grid);
out.S_est = surface_refiner(out.S_est);
out.S_true = GPIS(end).X;

% --- Create CB
CB = mvnrnd(out.est,out.var,option.CB_run);
for i=1:option.CB_run
    CB_temp = field_to_surface(out.thres,CB(i,:)',spatial_grid);
    out.CB{i} = surface_refiner(CB_temp);
    
end

if (option.ifStand == 1)
%%                  Convert 3-D model to unstandardized coordinates
    for i=1:3
        out.S_est(:,i) = out.S_est(:,i).*std_info(i).std + std_info(i).mean;
        out.S_true(:,i) = out.S_true(:,i).*std_info(i).std + std_info(i).mean;
        for j=1:option.CB_run
            out.CB{j}(:,i) = out.CB{j}(:,i).*std_info(i).std + std_info(i).mean;
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

[out.Haus,~] = HausdorffDist(out.S_est,out.S_true);

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



