function option = Configuration()
    option.gridsize = 40;
    option.pts_mode = 0; % 0: only on surface, 1: inner line included, 2: outter surface included
    option.cutoff = 1000;
    option.edgeLimit = 1; % 4e-3
    
    option.range_size = 500;
    option.thres_min = 0.4;
    option.thres_max = 0.6;
    option.thres_range = ...
            linspace(option.thres_min,option.thres_max,option.range_size)';
    
    option.bin_size = 2000;
    option.bin_min = 0;
    option.bin_max = 1;
    option.bin_range = ...
        linspace(option.bin_min,option.bin_max,option.bin_size)';
    
    option.num_worker = 10;
    
    option.ifStand = 1; % 0: not standardize data, 1: standardize data
    
    option.band_x = 2.5;
    option.band_y = 2.5;
    option.band_z = 8;
    option.band_f = 1;
end