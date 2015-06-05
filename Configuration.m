function option = Configuration()
    option.gridsize = 40;
    option.pts_mode = 0; % 0: only on surface, 1: inner line included, 2: outter surface included
    option.cutoff = 1000;
    option.edgeLimit = 4e-5;
    
    option.thres_step  = 1e-5;
    option.thres_min = -0.5e-3;
    option.thres_max = 0.5e-3;
    option.thres_range = (option.thres_min:option.thres_step:option.thres_max)';
    
    option.num_worker = 10;
    option.ifStand = 1; % 0: not standardize data, 1: standardize
end