fid = fopen('./cloud_point.txt','w');
for i=1:size(predict.S_est,1)
	 fprintf(fid,[num2str(predict.S_est(i,1)) ',' ...
				  num2str(predict.S_est(i,2)) ',' ...
				  num2str(predict.S_est(i,3))]);
     fprintf(fid,'\n');
end
fclose(fid);