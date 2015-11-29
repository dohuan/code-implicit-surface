function option = Configuration()
    option.gridsize = 10; % 20
    option.pts_mode = 0; % 0: only on surface, 1: inner line included, 2: outter surface included
    option.cutoff = 1000;
    option.edgeLimit = 1; % 4e-3
    
    option.thres_size = 100; % 500
    option.thres_span = 0.3; % 0.2
    
    option.bin_size = 2000;
    option.bin_min = 0;
    option.bin_max = 1;
    option.bin_range = ...
        linspace(option.bin_min,option.bin_max,option.bin_size)';
    
    option.num_worker = 10;
    
    option.ifStand = 1; % 0: not standardize data, 1: standardize data
    
    option.CB_run = 100; % number of MC run for credible band
    
    option.EM_run = 10;
    
end