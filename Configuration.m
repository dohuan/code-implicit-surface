function option = Configuration()
    option.gridsize = 40;
    option.pts_mode = 0; % 0: only on surface, 1: inner line included, 2: outter surface included
    option.cutoff = 1000;
    option.edgeLimit = 4e-3; % 4e-3
    
    option.thres_step  = 1e-4;
    option.thres_min = 0.85;
    option.thres_max = 0.9;
    option.thres_range = (option.thres_min:option.thres_step:option.thres_max)';
    
    option.num_worker = 10;
    
    option.ifStand = 0; % 0: not standardize data, 1: standardize
    
    option.band_x = 5;
    option.band_y = 5;
    option.band_z = 10;
    option.band_f = 1;
end