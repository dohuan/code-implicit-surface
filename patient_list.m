function Pat_list = patient_list()
mode = 1; % 0: run all 1: test one patient

if (mode==0)
    Pat_list(1).name = 'BB';
    Pat_list(1).numScan = 3;
    Pat_list(1).band_t = 360;
    Pat_list(1).band_x = 1.5;
    Pat_list(1).band_y = 1.5;
    Pat_list(1).band_z = 1.5;
    Pat_list(1).band_f = 1;
    
    Pat_list(2).name = 'DD';
    Pat_list(2).numScan = 3;
    Pat_list(2).band_t = 360;
    Pat_list(2).band_x = 1.9;
    Pat_list(2).band_y = 1.9;
    Pat_list(2).band_z = 1.9;
    Pat_list(2).band_f = 1;
    
    Pat_list(3).name = 'HH';
    Pat_list(3).numScan = 7;
    Pat_list(3).band_t = 360;
    Pat_list(3).band_x = 1.1;
    Pat_list(3).band_y = 1.1;
    Pat_list(3).band_z = 1.1;
    Pat_list(3).band_f = 1;
    
    Pat_list(4).name = 'II';
    Pat_list(4).numScan = 6;
    Pat_list(4).band_t = 360;
    Pat_list(4).band_x = 1.1;
    Pat_list(4).band_y = 1.1;
    Pat_list(4).band_z = 1.1;
    Pat_list(4).band_f = 1;
    
    Pat_list(5).name = 'JJ';
    Pat_list(5).numScan = 5;
    Pat_list(5).band_t = 360;
    Pat_list(5).band_x = 1.1;
    Pat_list(5).band_y = 1.1;
    Pat_list(5).band_z = 1.1;
    Pat_list(5).band_f = 1;
    
    Pat_list(6).name = 'KK';
    Pat_list(6).numScan = 5;
    Pat_list(6).band_t = 360;
    Pat_list(6).band_x = 1.1;
    Pat_list(6).band_y = 1.1;
    Pat_list(6).band_z = 1.1;
    Pat_list(6).band_f = 1;
    
    Pat_list(7).name = 'PP10';
    Pat_list(7).numScan = 4;
    Pat_list(7).band_t = 360;
    Pat_list(7).band_x = 0.9;
    Pat_list(7).band_y = 0.9;
    Pat_list(7).band_z = 0.9;
    Pat_list(7).band_f = 1;
    
    Pat_list(8).name = 'PP12';
    Pat_list(8).numScan = 6;
    Pat_list(8).band_t = 360;
    Pat_list(8).band_x = 1.1;
    Pat_list(8).band_y = 1.1;
    Pat_list(8).band_z = 1.1;
    Pat_list(8).band_f = 1;
    
    Pat_list(9).name = 'PP13';
    Pat_list(9).numScan = 4;
    Pat_list(9).band_t = 360;
    Pat_list(9).band_x = 1.1;
    Pat_list(9).band_y = 1.1;
    Pat_list(9).band_z = 1.1;
    Pat_list(9).band_f = 1;
    
    Pat_list(10).name = 'PP14';
    Pat_list(10).numScan = 3;
    Pat_list(10).band_t = 360;
    Pat_list(10).band_x = 1.1;
    Pat_list(10).band_y = 1.1;
    Pat_list(10).band_z = 1.1;
    Pat_list(10).band_f = 1;
    
else
    Pat_list(1).name = 'JJ';
    Pat_list(1).numScan = 5;
    Pat_list(1).band_t = 360;
    Pat_list(1).band_x = 1.1;
    Pat_list(1).band_y = 1.1;
    Pat_list(1).band_z = 1.1;
    Pat_list(1).band_f = 1;
    
    Pat_list(2).name = 'PP10';
    Pat_list(2).numScan = 4;
    Pat_list(2).band_t = 360;
    Pat_list(2).band_x = 0.9;
    Pat_list(2).band_y = 0.9;
    Pat_list(2).band_z = 0.9;
    Pat_list(2).band_f = 1;
    
%     Pat_list(1).name = 'JJ';
%     Pat_list(1).numScan = 5;
%     Pat_list(1).band_t = 360;
%     Pat_list(1).band_x = 1.1;
%     Pat_list(1).band_y = 1.1;
%     Pat_list(1).band_z = 1.1;
%     Pat_list(1).band_f = 1;
    
%     Pat_list(2).name = 'DD';
%     Pat_list(2).numScan = 3;
%     Pat_list(2).band_t = 360;
%     Pat_list(2).band_x = 1.9;
%     Pat_list(2).band_y = 1.9;
%     Pat_list(2).band_z = 1.9;
%     Pat_list(2).band_f = 1;
    
%     Pat_list(1).name = 'BB';
%     Pat_list(1).numScan = 3;
%     %Pat_list(1).band_t = 360;
%     Pat_list(1).band_t = 180;
%     Pat_list(1).band_x = 1.5;
%     Pat_list(1).band_y = 1.5;
%     Pat_list(1).band_z = 1.5;
%     Pat_list(1).band_f = 1;
    
end

end