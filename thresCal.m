function [thres,Haus_min,Haus_track] = thresCal(pat_name,est,S_true,spatial_grid,option,flag)
%%
% flag = 0 : search for upper threshold ONLY
% flag = 1 : search for upper AND lower threshold
isPlot = 1;


thres_lowBound = min(est);
thres_highBound = thres_lowBound + option.thres_span;
thres_range = linspace(thres_lowBound,thres_highBound,option.thres_size)';

Haus_track = zeros(size(thres_range,1),1);
S_est_min = [];
Haus_min = 1000;

for i=1:size(thres_range,1)
    S_est = [];
    
    low_bound = thres_lowBound;
    high_bound = thres_range(i);
    
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
            S_est_min = S_est;
            thres_min_up = thres_range(i);
        end
    end
    fprintf('Greedy search (up thres): Patient %s %.2f%%...\n',...
        pat_name,i/size(thres_range,1)*100);
    
end
thres.up = thres_min_up;
if (flag==1)
    index_temp = thres_range<=thres_min_up;
    thres_range = thres_range(index_temp);
    Haus_track_ = zeros(size(thres_range,1),1);
    S_est_min_ = [];
    Haus_min_ = 1000;
    for i=1:size(thres_range,1)
        S_est = [];
        
        low_bound = thres_range(i);
        high_bound = thres_min_up;
        
        for j=1:size(est,1)
            if (est(j,1)>=low_bound && est(j,1)<=high_bound)
                S_est = [S_est;spatial_grid(j,:)];
            end
        end
        if (isempty(S_est)==0)
            [Haus_temp,~] = HausdorffDist(S_est,S_true);
            Haus_track_(i) = Haus_temp;
            if (Haus_min_>Haus_temp)
                Haus_min_ = Haus_temp;
                S_est_min_ = S_est;
                thres_min_low = thres_range(i);
            end
        end
        fprintf('Greedy search (lower thres): Patient %s %.2f%%...\n',...
            pat_name,i/size(thres_range,1)*100);
    end
    thres.low = thres_min_low;
end
if (isPlot==1)
    plot(est);
    hold on;
    plot([0 7e4],[thres_min_up thres_min_up],'k-');
    if (flag==1)
        plot([0 7e4],[thres_min_low thres_min_low],'k--');
    end
    figure(2)
    scatter3(S_est_min(:,1),S_est_min(:,2),S_est_min(:,3));
    figure(3)
    scatter3(S_est_min_(:,1),S_est_min_(:,2),S_est_min_(:,3));
end
