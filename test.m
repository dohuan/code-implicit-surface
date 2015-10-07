figure('name','Cross-section view');
z_line = unique(predict.S_true(:,end),'rows','stable');

for i=1:size(z_line,1)
    pts_true(i).z = z_line(i);
    pts_true(i).xy = [];
end
for i=1:size(predict.S_true,1)
    ix = find(z_line==predict.S_true(i,end));
    pts_true(ix).xy = [pts_true(ix).xy; predict.S_true(i,1:2)];
end
hold on
for j=1:option.CB_run
    
    point_tmp = predict.CB{j};
    pts_CB = [];
    for i=1:size(z_line,1)
        pts_CB(i).z = z_line(i);
        pts_CB(i).xy = [];
    end
    for i=1:size(point_tmp,1)
        ix = find(z_line==point_tmp(i,end));
        pts_CB(ix).xy = [pts_CB(ix).xy; point_tmp(i,1:2)];
    end
    
    for k=1:length(z_line)
        subplot(4,5,k)
        plot(pts_CB(k).xy(:,1),pts_CB(k).xy(:,2),'r.','MarkerSize',2);
    end
end
for k=1:length(z_line)
    subplot(4,5,k)
    title(sprintf('z=%d',z_line(k)))
    plot(pts_true(k).xy(:,1),pts_true(k).xy(:,2),'r.','MarkerSize',2);
end
hold off