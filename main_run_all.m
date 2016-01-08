close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
% Estimate the point cloud and do thres search in 1 .m file
tic
addpath(genpath('./gpml'))
addpath(genpath('./HausdorffDist'))

ifSave = 0;
ifPlot = 0;
%Pat_list = patient_list_speed([],0);
Pat_list = patient_list_speed('HH',1);
option = Configuration();

if (ifSave==1)
    c = clock;
    folderName = ['./results/' num2str(c(1)) ...
        num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5)) '/'];
    mkdir(folderName);
    
    plySaveFolder = ['./results/PLY' num2str(c(1)) ...
        num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5)) '/'];
    mkdir(plySaveFolder);
    %fileID = fopen([folderName 'report.txt'],'w');
    
end

% --- Config parallel computing
% N = 10;
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     poolsize = 0;
% else
%     poolsize = poolobj.NumWorkers;
% end
% 
% if poolsize == 0
%     parpool('local',N);
% else
%     if poolsize~=N
%         delete(poolobj);
%         parpool('local',N);
%     end
% end

%predict(size(Pat_list,2)) = struct;
for i=1:size(Pat_list,2)
    %[predict, thres(i).a] = patient_process_speed(Pat_list(i),option);
    [predict, thres(i).a] = patient_process_speed(Pat_list(i),option);
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

%if (ifSave==1)
%    fclose(fileID);
%end

% if (ifPlot==1)
%     %% --- Visualize results
%     for i=1:size(Pat_list,2)
%         fig_name = sprintf('Patient %s',Pat_list(i).name);
%         figure('name',fig_name);
%         
%         subplot(1,2,1)
%         scatter3(predict(i).S_est(:,1),predict(i).S_est(:,2),predict(i).S_est(:,3));
%         title('Predicted');
%         subplot(1,2,2)
%         scatter3(predict(i).S_true(:,1),predict(i).S_true(:,2),predict(i).S_true(:,3));
%         title('True');
%     end
%     
%     % --- Plot side views of the TRUE and CB
%     if (size(Pat_list,2)==1)
%         figure('name','TRUE and CB')
%         subplot(1,2,1)
%         xlabel('X-Z view')
%         hold on
%         for i=1:option.CB_run
%             %h = scatter(predict.CB{i}(:,1),predict.CB{i}(:,3),'filled','SizeData',30);
%             plot(predict.CB{i}(:,1),predict.CB{i}(:,3),'b.','MarkerSize',5);
%             %pH = arrayfun(@(x) allchild(x),h);
%             %set(pH,'FaceAlpha',.01);
%         end
%         plot(predict.S_true(:,1),predict.S_true(:,3),'r.','MarkerSize',5);
%         hold off
%         
%         subplot(1,2,2)
%         hold on
%         for i=1:option.CB_run
%             plot(predict.CB{i}(:,2),predict.CB{i}(:,3),'b.','MarkerSize',5);
%             %alpha(h,0.1);
%         end
%         xlabel('Y-Z view')
%         plot(predict.S_true(:,2),predict.S_true(:,3),'r.','MarkerSize',5);
%         hold off
%         
%         % --- Plot side views of TRUE and MEAN
%         figure('name','TRUE and MEAN')
%         subplot(1,2,1)
%         xlabel('X-Z view')
%         hold on
%         plot(predict.S_true(:,1),predict.S_true(:,3),'r.','MarkerSize',5);
%         plot(predict.S_est(:,1),predict.S_est(:,3),'b.','MarkerSize',5);
%         hold off
%         
%         subplot(1,2,2)
%         xlabel('Y-Z view')
%         hold on
%         plot(predict.S_true(:,2),predict.S_true(:,3),'r.','MarkerSize',5);
%         plot(predict.S_est(:,2),predict.S_est(:,3),'b.','MarkerSize',5);
%         hold off
%     end
%     
%     
%     %for i=1:size(Pat_list,2)
%     %fprintf('Estimated time bandwith of patient %s: %.2f\n',Pat_list(i).name,predict(i).band_t);
%     %fprintf('Haus distance of patient %s: %.2f\n',Pat_list(i).name,predict(i).Haus_dist);
%     %end
%     
% end

time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);
fprintf('Here is the music!\n')
EndSound = load('handel');
sound(EndSound.y,EndSound.Fs);
%save('./results/run_all_071715_local')
