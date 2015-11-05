function out = patient_process_speed(pat_info,option)

fprintf('Progressing patient %s...\n',pat_info.name);

max_ = [0 0 0];
min_ = [1000 1000 1000];

time_convert = 30; % convert month to day

data_path = './Patient_Data/truncated_data/';
out.name = pat_info.name;
GPIS(pat_info.numScan) = struct('X',[],'y',[],'IS',[]);
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
    time_temp = time_convert*data.time_stamp*ones(size(X_temp,1),1);
    GPIS(i).X = [time_temp X_temp];
    X_all = [X_all; GPIS(i).X];
end

clear data
if (option.ifStand==1)
    %%                          Standardize data
    std_info(1).mean = mean(X_all(:,1));     % t
    std_info(1).std = sqrt(var(X_all(:,1)));
    std_info(2).mean = mean(X_all(:,2));     % x
    std_info(2).std = sqrt(var(X_all(:,2)));
    std_info(3).mean = mean(X_all(:,3));     % y
    std_info(3).std = sqrt(var(X_all(:,3)));
    std_info(4).mean = mean(X_all(:,4));     % z
    std_info(4).std = sqrt(var(X_all(:,4)));

    % --- Standardize data
    for k=1:pat_info
    for i=1:size(X,2)
        X(:,i) = (X(:,i)-std_info(i).mean)./std_info(i).std;
    end
    
end


end