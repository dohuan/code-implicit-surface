function [S_test,thres_min] = bin_search(est_train,est_test,thres,spatial_grid,option,flag)
%%
% flag = 0 : search for upper threshold ONLY
% flag = 1 : search for upper AND lower threshold

binvec_train = zeros(size(est_train,1),1);
index = est_train<=thres.up;
binvec_train(index) = 1;
Jac_min = 1000;
Jac_track = zeros(option.bin_size,1);
for i=1:option.bin_size
    binvec_test = zeros(size(est_test,1),1);
    index_temp = est_test<=option.bin_range(i);
    binvec_test(index_temp) = 1;
    Jac_track(i) = pdist([binvec_train';binvec_test'],'Jaccard');
    
    if (Jac_min>Jac_track(i))
        Jac_min = Jac_track(i);
        thres_min = option.bin_range(i);
    end
end

if (flag==1)
    
else
    S_test = [];
    for j=1:size(est_test,1)
        if (est_test(j,1)>=0&&est_test(j,1)<=thres_min)
            S_test = [S_test;spatial_grid(j,:)];
        end
    end
end
end

function [thres,Jac_track] = Jac_search(est,option)
binvec = zeros(size(est,1),1);
index = est<=thres.up;
binvec_train(index) = 1;
Jac_min = 1000;
Jac_track = zeros(option.bin_size,1);

for i=1:option.bin_size
    binvec_test = zeros(size(est_test,1),1);
    index_temp = est_test<=option.bin_range(i);
    binvec_test(index_temp) = 1;
    Jac_track(i) = pdist([binvec_train';binvec_test'],'Jaccard');
    
    if (Jac_min>Jac_track(i))
        Jac_min = Jac_track(i);
        thres_min = option.bin_range(i);
    end
end

end