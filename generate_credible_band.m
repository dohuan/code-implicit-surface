function generate_credible_band()
    est_list = make_list();
    for i=1:size(est_list,2)
        load(est_list{i});
        
    end
end

function out = make_list()
    out{1} = 'BB_201587_917';
%     out{2} = 'DD_201587_919';
%     out{3} = 'HH_201587_101';
%     out{4} = 'II_201587_1022';
%     out{5} = 'JJ_201587_1030';
%     out{6} = 'KK_201587_1038';
%     out{7} = 'PP10_201587_1043';
%     out{8} = 'PP12_201587_115';
%     out{9} = 'PP13_201587_118';
%     out{10} = 'PP14_201587_1111';
end