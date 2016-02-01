close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
% Estimate the point cloud and do thres search in 1 .m file
tic
addpath(genpath('./gpml'))
addpath(genpath('./HausdorffDist'))

ifSave = 1;
ifPlot = 0;
%Pat_list = patient_list_speed([],0);
Pat_list = patient_list_speed('II',1);
option = Configuration();

savetag = 'save1';

if (ifSave==1)
    c = clock;
    folderName = ['./results/' savetag  num2str(c(1)) ...
        num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5)) '/'];
    mkdir(folderName);
    
    plySaveFolder = ['./results/PLY' savetag  num2str(c(1)) ...
        num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5)) '/'];
    mkdir(plySaveFolder);
    %fileID = fopen([folderName 'report.txt'],'w');
    
end

% --- Config parallel computing
%N = 10;
%poolobj = gcp('nocreate');
%if isempty(poolobj)
%    poolsize = 0;
%else
%    poolsize = poolobj.NumWorkers;
%end

%poolsize = matlabpool('size');
%if poolsize == 0
%    matlabpool('local',N);
%else
%    if poolsize~=N
%        matlabpool(close);
%        matlabpool('local',N);
%    end
%end

%predict(size(Pat_list,2)) = struct;
for i=1:size(Pat_list,2)
    %[predict, thres(i).a] = patient_process_speed(Pat_list(i),option);
    predict = patient_process_speed_1(Pat_list(i),option);
    if (ifSave==1)
        saveFile = [folderName Pat_list(i).name];
        %save(filename,'dumvar','-v7.3');
        %save(saveFile,'-v7.3');
        
        % --- Export to PLY files
        for j=1:option.CB_run
            plySave = [plySaveFolder Pat_list(i).name '_'  num2str(j) '.ply'];
            exportMesh(predict.CB{j},plySave);
        end
        exportMesh(predict.S_est,[plySaveFolder Pat_list(i).name '_est.ply']);
        % --- Save images
        h = figure(i);
        hold on
        scatter3(predict.S_true(:,1),predict.S_true(:,2),predict.S_true(:,3),'k+');
        for j=1:size(predict.CB,2)
            scatter3(predict.CB{j}(:,1),predict.CB{j}(:,2),predict.CB{j}(:,3),'b.');
        end
        hold off
        view(30,24);
        saveas(h,[folderName Pat_list(i).name '.jpg']);
        close(h);
        % --- Save predict report to text file
        fileID = fopen([folderName 'report' Pat_list(i).name '.txt'],'w');
        if (predict.method ==1)
            fprintf(fileID,'Patient ID,sig_t,sig_f,sig_x,sig_y,sig_z,Haust_dist, thres of test\n');
            fprintf(fileID,'%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f\n',...
                predict.name,predict.band_t,predict.band_f,predict.band_x,...
                predict.band_y,predict.band_z,predict.Haus_dist,predict.thres_test);
        else
            fprintf(fileID,'Patient ID, Haust_dist, threshold value\n');
            fprintf(fileID,'%s,%.2f,%.4f\n',...
                predict.name, predict.Haus,predict.thres);
        end
        fclose(fileID);
        %clear predict;
    end
end

%matlabpool close;

time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);
fprintf('Here is the music!\n')
EndSound = load('handel');
sound(EndSound.y,EndSound.Fs);
