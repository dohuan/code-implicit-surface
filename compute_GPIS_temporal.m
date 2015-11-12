function data = compute_GPIS_temporal(test_time,data,spatial_grid,pat_info,std_info,time_info,option)
%% Input data ONLY up to time t
nt = length(data); % available data up to time t
X = [];
y = [];
for i=1:nt
	time_tmp = data(i).scantime*ones(size(data(i).X,1),1);
	X = [X;[time_tmp data(i).X]];
	y = [y;data(i).y];
end
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

hyp = minimize(hyp, @gp, -10, @infExact, meanfunc, covfunc, ...
                                         likfunc, X, y);
test_grid = [];
for i=1:length(test_time)
	test_grid = [test_grid;[test_time(i)*ones(size(spatial_grid,1),1) spatial_grid]];
end
[est, var] = ...
       gaussian_process(hyp, covfunc,meanfunc, data.X,data.y,test_grid);
data.IS.est = est;
data.IS.var = var;
data.IS.hyp = hyp;
end