close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
% --- Run patient H with different window horizon
tic
addpath(genpath('./gpml'))
addpath(genpath('./HausdorffDist'))

ifSave = 1;
ifPlot = 0;
Pat_list = patient_list_speed('HH',1);
option = Configuration();

savetag = 'EM_patient_H';
horizon_window = [6 5 4 3];

for t=1:length(horizon_window)
    Pat_list.start_ix = Pat_list.numScan - horizon_window(t);
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
        predict = patient_process_speed_H(Pat_list(i),option);
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
    Haus_track(t) = predict.Haus;
end
%matlabpool close;

x = horizon_window;
y = Haus_track;

close all
figure(1);
title('1st order fitting: patient H');
hold on
p = polyfit(x,y,1);
yfit =  p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
t = 0:1:8;
yplot = p(1)*t + p(2);
std_error = sqrt(SSresid/length(y));
fprintf('1st order fitting:\n');
fprintf('R square: %.3f\n',rsq);
fprintf('Standard Error: %.3f\n',std_error);
plot(t,yplot,'k','LineWidth',2);
%xlim([2 8]);
xlabel('Number of scans','FontSize',14);
ylabel('Hausdoff distance','FontSize',14);
box on
set(gca,'FontSize',16);
saveas(figure(1),[folderName 'Haus_vs_scans_1storder.eps']);

figure(2)
title('2nd order fitting: patient H');
hold on
p2 = polyfit(x,y,2);
yfit = p2(1)*x.^2+p2(2)*x+p2(3);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
t = 0:1:8;
yplot = p2(1)*t.^2 + p2(2)*t + p2(3);
std_error = sqrt(SSresid/length(y));
fprintf('2nd order fitting:\n');
fprintf('R square: %.3f\n',rsq);
fprintf('Standard Error: %.3f\n',std_error);
plot(t,yplot,'k','LineWidth',2);
xlim([2 8]);
xlabel('Number of scans','FontSize',14);
ylabel('Hausdoff distance','FontSize',14);
box on
set(gca,'FontSize',16);
saveas(figure(2),[folderName 'Haus_vs_scans_2ndorder.eps']);


time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);
fprintf('Here is the music!\n')
EndSound = load('handel');
sound(EndSound.y,EndSound.Fs);
