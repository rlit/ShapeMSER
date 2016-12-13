function [evecs evals] = prepare_evec(L,n)

e = sparse([1:length(L)], [1:length(L)], ones([1,length(L)]));
n = min(size(L,1),n);
[evecs evals] = eigs(L+e*1e-8 , n ,'sm', struct('disp',0));
evals = diag(evals);
[evals,idx] = sort(evals, 'ascend');
evecs = evecs(:,idx);
