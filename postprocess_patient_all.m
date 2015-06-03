close all
clear
clc
%%
addpath(genpath('./HausdorffDist'))

Pat_list(1).name = 'H';
Pat_list(1).numScan = 7;
% Pat_list(2).name = 'J';
% Pat_list(2).numScan = 5;
% Pat_list(3).name = 'K';
% Pat_list(3).numScan = 5;

cutoff = 2000;
thres_step  = 1e-5;
thres_min = 1.5e-3;
thres_max = 1.9e-3;
thres_range = (thres_min:thres_step:thres_max)';

for p=1:size(Pat_list,2)
    % --- Load true cloud points
    fprintf('Progressing patient %s...\n',Pat_list(p).name);
    file_name = ['./Patient_Data/HPCC_data/'...
        Pat_list(p).name Pat_list(p).name num2str(Pat_list(p).numScan)...
        '_inner'];
    load(file_name); % load 'data' after this line
    % --- Load estimated cloud points
    file_name = ['./HPCC_est/' Pat_list(p).name Pat_list(p).name...
        '_HPCC_est051115'];
    load(file_name); 
    
    [S1,S2,S3] = meshgrid(option.x_mesh,option.y_mesh,option.z_mesh); % [middle shortest longest]
    S_temp = [S1(:),S2(:),S3(:)];
    
    
    index = randperm(size(data.on_surface,1));
    S_true = data.on_surface(index(1:cutoff),1:3);
    Haus_track = zeros(size(thres_range,1),1);
    S_est_min = [];
    Haus_min = 1000;
    thres_min_temp = 0;
    for i=1:size(thres_range,1)
        S_est = [];
        for j=1:size(est,1)
            if (est(j,1)>=(thres_min-thres_step)&&est(j,1)<=thres_range(i))
                S_est = [S_est;S_temp(j,:)];
            end
        end
        [Haus_temp,~] = HausdorffDist(S_est,S_true);
        Haus_track(i) = Haus_temp;
        if (Haus_min>Haus_temp)
            Haus_min = Haus_temp;
            S_est_min = S_est;
            thres_min_temp = thres_range(i);
        end
        fprintf('Patient %s %.2f%%...\n',...
            Pat_list(p).name,i/size(thres_range,1)*100);
    end
    predict(p).S_est = S_est_min;
    predict(p).thres = thres_min_temp;
end

for i=1:size(predict,2)
    cm = colormap(pink);
    cid = repmat(cm(1,:),size(predict(i).S_est,1),1);
    markerSize = 25*ones(size(predict(i).S_est,1),1);
    scatter3(predict(i).S_est(:,1),predict(i).S_est(:,2),predict(i).S_est(:,3)...
        ,markerSize,cid);
end


