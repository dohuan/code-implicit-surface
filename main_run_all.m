close all
clear
clc
%%
% Estimate the point cloud and do thres search in 1 .m file
tic
addpath(genpath('./gpml'))
addpath(genpath('./HausdorffDist'))

% Pat_list(1).name = 'BB';
% Pat_list(1).numScan = 3;
% Pat_list(1).band_t = 550;
% Pat_list(2).name = 'DD';
% Pat_list(2).numScan = 3;
% Pat_list(2).band_t = 700;
% Pat_list(3).name = 'HH';
% Pat_list(3).numScan = 7;
% Pat_list(3).band_t = 85;
% Pat_list(4).name = 'II';
% Pat_list(4).numScan = 6;
% Pat_list(4).band_t = 80;
% Pat_list(5).name = 'JJ';
% Pat_list(5).numScan = 5;
% Pat_list(5).band_t = 85;
% Pat_list(6).name = 'KK';
% Pat_list(6).numScan = 5;
% Pat_list(6).band_t = 85;
% Pat_list(7).name = 'PP10';
% Pat_list(7).numScan = 4;
% Pat_list(7).band_t = 350;
% Pat_list(8).name = 'PP12';
% Pat_list(8).numScan = 6;
% Pat_list(8).band_t = 85;
% Pat_list(9).name = 'PP13';
% Pat_list(9).numScan = 4;
% Pat_list(9).band_t = 80;
% Pat_list(10).name = 'PP14';
% Pat_list(10).numScan = 3;
% Pat_list(10).band_t = 80;

% Pat_list(1).name = 'PP10';
% Pat_list(1).numScan = 4;
% Pat_list(1).band_t = 450; % 350 12.87

Pat_list(1).name = 'II';
Pat_list(1).numScan = 5;
Pat_list(1).band_t = 500;

opt = Configuration();

parfor i=1:size(Pat_list,2)
	predict(i) = patient_process(Pat_list(i),opt);
end

%% Visualize results

for i=1:size(Pat_list,2)
	fig_name = sprintf('Patient %s',Pat_list(i).name);
    figure('name',fig_name);
    
%     cm = colormap(pink);
%     cid = repmat(cm(1,:),size(predict(i).S_est,1),1);
%     markerSize = 25*ones(size(predict(i).S_est,1),1);
%     scatter3(predict(i).S_est(:,1),predict(i).S_est(:,2),predict(i).S_est(:,3)...
%         ,markerSize,cid);
    subplot(1,2,1)
    scatter3(predict(i).S_est(:,1),predict(i).S_est(:,2),predict(i).S_est(:,3));
    title('Predicted');
    subplot(1,2,2)
    scatter3(predict(i).S_true(:,1),predict(i).S_true(:,2),predict(i).S_true(:,3));
    title('True');
    
    fprintf('Haus distance of patient %s: %.2f\n',Pat_list(i).name,predict(i).Haus_dist);
end
time_run = toc;
fprintf('\nRun time: %.2f minutes',time_run/60);

%save('run_all_060115_3')