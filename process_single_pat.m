function process_single_pat(pat_name)
%% Process and save each patient separately

tic
addpath(genpath('./gpml'))
addpath(genpath('./HausdorffDist'))

Pat_list = patient_list(pat_name,1);
option = Configuration();

predict = patient_process(Pat_list,option);
c = clock;
saveFile = ['./results/' pat_name '_' num2str(c(1)) num2str(c(2)) num2str(c(3)) '_' num2str(c(4)) num2str(c(5))];

save(saveFile);
end