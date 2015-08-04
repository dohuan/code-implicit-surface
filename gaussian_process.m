function [mu,sig] = gaussian_process(hyp,covfunc,meanfunc,x,y,xs)
%% Gaussian process for discrete target
K = covariance_matrix(hyp,covfunc,x,xs);
C = covariance_matrix(hyp,covfunc,x);
Sig_0 = covariance_matrix(hyp,covfunc,xs);

Cy = C\(y - meanfunc([],x));
CK = C\K;
mu = meanfunc([],xs) + K'*Cy;
sig = Sig_0 - K'*CK;
end

function C = covariance_matrix(hyp,covfunc,x,x_)
%% x and x_ are [nt x ns]
% nt: number of variates
% ns: dimension of variates
loghyper = [];
for i=1:size(hyp,2)
    loghyper = [loghyper;hyp(i).cov];
end
if nargin<4
    C = feval(covfunc{:},loghyper,x);
else
    C = feval(covfunc{:},loghyper,x,x_);
end
end