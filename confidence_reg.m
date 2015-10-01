function confidence_reg(file_name)
%%

name_list{1} = 'HH';
name_list{2} = 'II';
name_list{3} = 'JJ';
name_list{4} = 'KK';
name_list{5} = 'P0';
name_list{6} = 'P2';
name_list{7} = 'P3';

data = xlsread(file_name,'Sheet1');
z_index = 1.645; % 90% confidence interval

for i=1:size(name_list,2)
    subplot(2,4,i)
    temp = data(:,i);
    ix = isnan(data(:,i));
    temp(ix) = [];
    
    if (i~=1)
        ix = [];
        for j=1:size(temp,1)
            if (temp(j)>26)
                ix = [ix;j];
            end
        end
        temp(ix) = [];
    end
    
    %SEM = std(temp)/sqrt(length(temp));
    %SEM = std(temp);
    %CI = [mean(temp)-z_index*SEM, mean(temp)+z_index*SEM];
    
    CI(1) = quantile(temp,.05);
    CI(2) = quantile(temp,.95);
    [h,x] = hist(temp,100);
    h = h./sum(h);
    bar(x,h);
    title([name_list{i} '(90)']);
    hold on
    %plot([mean(temp) mean(temp)],[0 max(h)],'r-','LineWidth',2);
    plot([CI(1) CI(1)],[0 max(h)],'r:','LineWidth',2);
    plot([CI(2) CI(2)],[0 max(h)],'r:','LineWidth',2);
    hold off
    axis tight
    box on
    % --- Check confidence region
%     psum = 0;
%     for j=1:size(x,2)
%         if (x(j)>=CI(1)&&x(j)<=CI(2))
%             psum = psum + h(j);
%         end
%     end
%     fprintf('Psum of patient %s is: %.2f\n',name_list{i},psum);
end

end
