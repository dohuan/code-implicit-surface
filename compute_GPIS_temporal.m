function data = compute_GPIS_temporal(test_time,data,spatial_grid,pat_info,std_info,time_info,option)

covfunc  = {'covSum', {'covSEard','covNoise'}};
likfunc  = @likGauss;
meanfunc = @meanOne;
hyp.cov(1) = (pat_info.band_t-time_info.min)/(time_info.max-time_info.min);
if (option.ifStand == 1)
    %hyp.cov(1) = log(abs((pat_info.band_t)/std_info(1).std));   % bandwidth of time
    hyp.cov(2) = log(abs((pat_info.band_x)/std_info(1).std));  % bandwidth of x
    hyp.cov(3) = log(abs((pat_info.band_y)/std_info(2).std));  % bandwidth of y
    hyp.cov(4) = log(abs((pat_info.band_z)/std_info(3).std));  % bandwidth of z
else
    %hyp.cov(1) = log(pat_info.band_t);
    hyp.cov(2) = log(pat_info.band_x);
    hyp.cov(3) = log(pat_info.band_y);
    hyp.cov(4) = log(pat_info.band_z);
end
hyp.cov(5) = log(pat_info.band_f);   % \sig_f
hyp.cov(6) = log(0.05);              % \sig_w: measurement noise .05
hyp.lik = log(0.05);                 % initial prior probability for hyper p(\theta) 0.03


end