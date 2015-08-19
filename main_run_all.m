close all
clear
clc
%%
% Estimate the point cloud and do thres search in 1 .m file
tic
addpath(genpath('./gpml'))
addpath(genpath('./HausdorffDist'))

ifPlot = 0;
ifSave = 0;
%Pat_list = patient_list([],0);
Pat_list = patient_list('II',1);
option = Configuration();

if (ifSave==1)
    c = clock;
    folderName = ['./results/' num2str(c(1)) ...
        num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5)) '/'];
    mkdir(folderName);

    plySaveFolder = ['./results/PLY' num2str(c(1)) ...
            num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5)) '/'];
    mkdir(plySaveFolder);
end
for i=1:size(Pat_list,2)
    predict = patient_process(Pat_list(i),option);
    %predict(i) = patient_process_1(Pat_list(i),option);
    
    if (ifSave==1)
        saveFile = [folderName Pat_list(i).name];
        save(saveFile);

        % --- Export to PLY files
        for j=1:option.CB_run
            plySave = [plySaveFolder Pat_list(i).name '_'  num2str(j) '.ply'];
            exportMesh(predict.CB{j},plySave);
        end

        clear predict
    end
end

if (ifPlot==1)
    %% Visualize results
    for i=1:size(Pat_list,2)
        fig_name = sprintf('Patient %s',Pat_list(i).name);
        figure('name',fig_name);
        
        subplot(1,2,1)
        scatter3(predict(i).S_est(:,1),predict(i).S_est(:,2),predict(i).S_est(:,3));
        title('Predicted');
        subplot(1,2,2)
        scatter3(predict(i).S_true(:,1),predict(i).S_true(:,2),predict(i).S_true(:,3));
        title('True');
    end
    for i=1:size(Pat_list,2)
        fprintf('Estimated time bandwith of patient %s: %.2f\n',Pat_list(i).name,predict(i).band_t);
        fprintf('Haus distance of patient %s: %.2f\n',Pat_list(i).name,predict(i).Haus_dist);
    end
end

time_run = toc;
fprintf('\nRun time: %.2f minutes',time_run/60);

%save('./results/run_all_071715_local')
