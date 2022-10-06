function [part] = make_xval_partition(n, n_folds)
% MAKE_XVAL_PARTITION - Randomly generate cross validation partition.
%
% Usage:
%
%  PART = MAKE_XVAL_PARTITION(N, N_FOLDS)
%
% Randomly generates a partitioning for N datapoints into N_FOLDS equally
% sized folds (or as close to equal as possible). PART is a 1 X N vector,
% where PART(i) is a number in (1...N_FOLDS) indicating the fold assignment
% of the i'th data point.

% YOUR CODE GOES HERE

s=mod(n,n_folds);r=n-s;%mod 求余数
p1=ceil((1:r)/ceil(r/n_folds));
p2=randperm(n_folds);p2=p2(1:s);%randperm将一列序号随机打乱
p3=[p1 p2];
part=p3(randperm(size(p3,2)));
end