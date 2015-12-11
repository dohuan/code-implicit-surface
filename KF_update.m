function [mean_T,cov_T,cov_T_,mean_T0,cov_T0,cov_T_10,logL_] =...
    KF_update(A,SW,mu0,Sig0,data,pat_info,spatial_grid)
T = pat_info.numScan-1; % length of available data
ns = size(spatial_grid,1); % length of all elements of the grid
est = zeros(ns,T,T);
cov  = zeros(ns,ns,T,T);
J = zeros(ns,ns,T);

% --- Forward KF ---
for t=1:T
    if (t==1)
        delta_t = data(t).scantime;
        mean_ = mu0 + A*delta_t;    % mean(f(1)|D_0)
        cov_ = Sig0 + SW*delta_t^2; % cov(f(1)|D_0)
        K = cov_*(cov_+data(t).IS.var)^(-1);
        
        est(:,t,t) = mean_ + K*(data(t).IS.est-mean_);
        cov(:,:,t,t) = (eye(size(K,1))-K)*cov_;
    else
        delta_t = data(t).scantime-data(t-1).scantime;
        est(:,t,t-1) = est(:,t-1,t-1) + A*delta_t;
        cov(:,:,t,t-1) = cov(:,:,t-1,t-1) + SW*delta_t^2;
        K = cov(:,:,t,t-1)*(cov(:,:,t,t-1)+data(t).IS.var)^(-1);
        
        est(:,t,t) = est(:,t,t-1) + K*(data(t).IS.est-est(:,t,t-1));
        cov(:,:,t,t) = (eye(size(K,1))-K)*cov(:,:,t,t-1);
    end
end

mean_T = zeros(ns,T);
cov_T = zeros(ns,ns,T);
cov_T_ = zeros(ns,ns,T);

mean_T(:,end) = est(:,end,end);
cov_T(:,:,end) = cov(:,:,end,end);
cov_T_(:,:,T) = (eye(size(K,1))-K)*cov(:,:,T-1,T-1);

% --- Backward computation ---
for t=T:-1:1
    if (t==1)
        J0 = Sig0*cov_^(-1);
        mean_T0 = mu0 + J0*(mean_T(:,t)-mu0);    % E(f(0)|D_T) take mu0 as this one
        cov_T0  = Sig0 + J0*(cov_T(:,:,t)-cov_)*J0';
    else
        J(:,:,t-1) = cov(:,:,t-1,t-1)*(cov(:,:,t,t-1))^(-1);
        mean_T(:,t-1) = est(:,t-1,t-1)+J(:,:,t-1)*(mean_T(:,t)-est(:,t-1,t-1));
        cov_T(:,:,t-1) = cov(:,:,t-1,t-1)+...
            J(:,:,t-1)*(cov_T(:,:,t)-cov(:,:,t,t-1))*J(:,:,t-1)';
        if (t~=2)
            cov_T_(:,:,t-1) = cov(:,:,t-1,t-1)*J(:,:,t-2)'+...
                J(:,:,t-1)*(cov_T_(:,:,t)-cov(:,:,t-1,t-1))*J(:,:,t-2)';
        else
            J0 = Sig0*cov_^(-1);
            cov_T_10 = cov(:,:,t-1,t-1)*J0'+...
                J(:,:,t-1)*(cov_T_(:,:,t)-cov(:,:,t-1,t-1))*J0';
        end
    end
end

% --- Compute MODIFIED logL as in (Gupta(1974))

logL_ = 0;
for t=2:T %(t starts from 2 to avoid 1/0 error)
    delta_t = data(t).scantime-data(t-1).scantime;
    logL_ = logL_...
        -1/2*logdet(cov(:,:,t,t-1)+delta_t^2*SW)...
        -1/2*((data(t).IS.est-est(:,t,t-1))'...
        *(cov(:,:,t,t-1)+delta_t^2*SW)^(-1)*(data(t).IS.est-est(:,t,t-1)));
end

end