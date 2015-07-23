clear
clc
close all
%% take output as cloud point, convert to PLY for 3-D rendering
cloudFile = './Results/run_all_062115_local(final)';
load(cloudFile);

for i=1:size(predict,2)
    fileSaveName = ['./PLY_data/' Pat_list(i).name '_test.ply'];
    exportMesh(predict(i).S_true,fileSaveName);
    
    fileSaveName = ['./PLY_data/' Pat_list(i).name '_est.ply'];
    exportMesh(predict(i).S_est,fileSaveName);
end