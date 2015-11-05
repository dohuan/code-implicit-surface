function Pat_list = patient_list(pat_name, mode)
% mode = 0; % 0: run all 1: test one patient

% Pat_list(1).name = 'BB';
% Pat_list(1).numScan = 3;
% Pat_list(1).band_t = 360;
% Pat_list(1).band_x = 1.5;
% Pat_list(1).band_y = 1.5;
% Pat_list(1).band_z = 1.5;
% Pat_list(1).band_f = 1;
% 
% Pat_list(2).name = 'DD';
% Pat_list(2).numScan = 3;
% Pat_list(2).band_t = 360;
% Pat_list(2).band_x = 1.9;
% Pat_list(2).band_y = 1.9;
% Pat_list(2).band_z = 1.9;
% Pat_list(2).band_f = 1;

Pat_list(1).name = 'HH';
Pat_list(1).ID = 'H';
Pat_list(1).numScan = 7; %7
Pat_list(1).band_t = 210; %
Pat_list(1).band_x = 1.05; % 1.05
Pat_list(1).band_y = 1.05; % 1.05
Pat_list(1).band_z = 1.05; % 1.05
Pat_list(1).band_f = 0.95; % 0.95

Pat_list(2).name = 'II';
Pat_list(2).ID = 'I';
Pat_list(2).numScan = 6;
Pat_list(2).band_t = 270;
Pat_list(2).band_x = 1.0;
Pat_list(2).band_y = 1.0;
Pat_list(2).band_z = 1.0;
Pat_list(2).band_f = 0.9;

Pat_list(3).name = 'JJ';
Pat_list(3).ID = 'J';
Pat_list(3).numScan = 5;
Pat_list(3).band_t = 400;
Pat_list(3).band_x = 1.0;
Pat_list(3).band_y = 1.0;
Pat_list(3).band_z = 1.0;
Pat_list(3).band_f = 0.95;

Pat_list(4).name = 'KK';
Pat_list(4).ID = 'K';
Pat_list(4).numScan = 5;
Pat_list(4).band_t = 270;
Pat_list(4).band_x = 1.0; % 1.0 1.1 .95 .8 1.2
Pat_list(4).band_y = 1.0;
Pat_list(4).band_z = 1.0;
Pat_list(4).band_f = 0.9;

Pat_list(5).name = 'P0';
Pat_list(5).ID = 'P10';
Pat_list(5).numScan = 4;
Pat_list(5).band_t = 360;
Pat_list(5).band_x = 0.9;
Pat_list(5).band_y = 0.9;
Pat_list(5).band_z = 0.9;
Pat_list(5).band_f = 1;

Pat_list(6).name = 'P2';
Pat_list(6).ID = 'P12';
Pat_list(6).numScan = 6;
Pat_list(6).band_t = 360;
Pat_list(6).band_x = 1.0;
Pat_list(6).band_y = 1.0;
Pat_list(6).band_z = 1.0;
Pat_list(6).band_f = 0.95;

Pat_list(7).name = 'P3';
Pat_list(7).ID = 'P13';
Pat_list(7).numScan = 4;
Pat_list(7).band_t = 270;
Pat_list(7).band_x = 1.0; % 1.1
Pat_list(7).band_y = 1.0; % 1.1
Pat_list(7).band_z = 1.0; % 1.1 0.65
Pat_list(7).band_f = 1.0; % 1.0

% Pat_list(10).name = 'PP14';
% Pat_list(10).numScan = 3;
% Pat_list(10).band_t = 360;
% Pat_list(10).band_x = 1.1;
% Pat_list(10).band_y = 1.1;
% Pat_list(10).band_z = 1.1;
% Pat_list(10).band_f = 1;

if (mode==1)
    for i=1:size(Pat_list,2)
        if (strcmp(Pat_list(i).name,pat_name)==1)
            Pat_list = Pat_list(i);
            break;
        end
    end
end

end