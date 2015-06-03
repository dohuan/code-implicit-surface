function [thres_min,Haus_min,Haus_track] = thresCal(pat_name,est,S_true,spatial_grid,option,flag,thres_up)
%%

Haus_track = zeros(size(option.thres_range,1),1);
%S_est_min = [];
Haus_min = 1000;
if (flag==1)
    for i=1:size(option.thres_range,1)
        S_est = [];
        
        low_bound = option.thres_min-option.thres_step;
        high_bound = option.thres_range(i);
        
        for j=1:size(est,1)
            if (est(j,1)>=low_bound && est(j,1)<=high_bound)
                S_est = [S_est;spatial_grid(j,:)];
            end
        end
        if (isempty(S_est)==0)
            [Haus_temp,~] = HausdorffDist(S_est,S_true);
            Haus_track(i) = Haus_temp;
            if (Haus_min>Haus_temp)
                Haus_min = Haus_temp;
                %S_est_min = S_est;
                thres_min = option.thres_range(i);
            end
        end
        fprintf('Greedy search 1: Patient %s %.2f%%...\n',...
            pat_name,i/size(option.thres_range,1)*100);
        
    end
else
    index_temp = option.range<=thres_up;
    thres_range = option.thres_range(index_temp);
    
    for i=1:size(thres_range,1)
        S_est = [];
        
        low_bound = thres_range(i);
        high_bound = thres_up;
        
        for j=1:size(est,1)
            if (est(j,1)>=low_bound && est(j,1)<=high_bound)
                S_est = [S_est;spatial_grid(j,:)];
            end
        end
        if (isempty(S_est)==0)
            [Haus_temp,~] = HausdorffDist(S_est,S_true);
            Haus_track(i) = Haus_temp;
            if (Haus_min>Haus_temp)
                Haus_min = Haus_temp;
                %S_est_min = S_est;
                thres_min = thres_range(i);
            end
        end
        fprintf('Greedy search 2: Patient %s %.2f%%...\n',...
            pat_name,i/size(thres_range,1)*100);
    end
end
