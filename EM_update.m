function [A,SW, logL] = EM_update(mean_T,cov_T,cov_T_,mean_T0,cov_T0,...
    cov_T_10,Sig0,data,pat_info)

T = pat_info.numScan-1; % length of available data
tmp1 = 0;
tmp2 = zeros(size(mean_T0,1),1);
tmp3 = zeros(size(cov_T0,1),size(cov_T0,2));
% --- Compute A: (t starts from 2 to avoid 1/0 error)
for t=2:T
    if (t==1)
        delta_t = data(t).scantime;
            
        tmp1 = tmp1 + delta_t^4;
        tmp2 = tmp2 + delta_t^3*(mean_T(:,t)-mean_T0);
    else
        delta_t = data(t).scantime-data(t-1).scantime;

        tmp1 = tmp1 + delta_t^4;
        tmp2 = tmp2 + delta_t^3*(mean_T(:,t)-mean_T(:,t-1));
    end
end

A = 1/tmp1*tmp2;

% --- Compute SW: (t starts from 2 to avoid 1/0 error)
for t=2:T
    if (t==1)
        delta_t = data(t).scantime;
        
        Ht = (mean_T0-mean_T(:,t))*A'+A*(mean_T0-mean_T(:,t))';
        Ct = cov_T(:,:,t) + mean_T(:,t)*mean_T(:,t)';
        Et = cov_T0 + mean_T0*mean_T0';
        Bt = cov_T_10 + mean_T(:,t)*mean_T0';
        
        tmp3 = tmp3 + delta_t^(-2)*(Ct-Bt-Bt'+Et+delta_t*Ht+delta_t^2*(A*A'));
    else
        delta_t = data(t).scantime-data(t-1).scantime;
        
        Ht = (mean_T(:,t-1)-mean_T(:,t))*A'+A*(mean_T(:,t-1)-mean_T(:,t))';
        Ct = cov_T(:,:,t) + mean_T(:,t)*mean_T(:,t)';
        Et = cov_T(:,:,t-1) + mean_T(:,t-1)*mean_T(:,t-1)';
        Bt = cov_T_(:,:,t-1) + mean_T(:,t)*mean_T(:,t-1)';
        
        tmp3 = tmp3 + delta_t^(-2)*(Ct-Bt-Bt'+Et+delta_t*Ht+delta_t^2*(A*A'));
        
    end
end

SW = 1/T*tmp3;

logL = -1/2*logdet(Sig0)-1/2*trace(Sig0^(-1)*(cov_T0+mean_T0*mean_T0'));
SWinv = SW^(-1);
% (t starts from 2 to avoid 1/0 error)
for t=2:T 
    if (t==1)
        delta_t = data(t).scantime;
        
        Ht = (mean_T0-mean_T(:,t))*A'+A*(mean_T0-mean_T(:,t))';
        Ct = cov_T(:,:,t) + mean_T(:,t)*mean_T(:,t)';
        Et = cov_T0 + mean_T0*mean_T0';
        Bt = cov_T_10 + mean_T(:,t)*mean_T0';
        
        logL = logL -1/2*logdet(delta_t^2*SW)...
            -1/2*trace(delta_t^(-2)*SWinv*(Ct-Bt-Bt'+Et+delta_t*Ht+delta_t^2*(A*A')))...
            -1/2*logdet(data(t).IS.var)...
            -1/2*trace((delta_t*data(t).IS.var)^(-1)*(cov_T(:,:,t)+...
              (data(t).IS.est-mean_T(:,t))*(data(t).IS.est-mean_T(:,t))'));
    else
        delta_t = data(t).scantime-data(t-1).scantime;
        
        Ht = (mean_T(:,t-1)-mean_T(:,t))*A'+A*(mean_T(:,t-1)-mean_T(:,t))';
        Ct = cov_T(:,:,t) + mean_T(:,t)*mean_T(:,t)';
        Et = cov_T(:,:,t-1) + mean_T(:,t-1)*mean_T(:,t-1)';
        Bt = cov_T_(:,:,t-1) + mean_T(:,t)*mean_T(:,t-1)';
        
        logL = logL -1/2*logdet(delta_t^2*SW)...
            -1/2*trace(delta_t^(-2)*SWinv*(Ct-Bt-Bt'+Et+delta_t*Ht+delta_t^2*(A*A')))...
            -1/2*logdet(data(t).IS.var)...
            -1/2*trace((delta_t*data(t).IS.var)^(-1)*(cov_T(:,:,t)+...
              (data(t).IS.est-mean_T(:,t))*(data(t).IS.est-mean_T(:,t))'));
        
    end
end

end