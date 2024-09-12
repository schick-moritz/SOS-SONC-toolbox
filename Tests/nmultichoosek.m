function combs = nmultichoosek(vecNum, k)
% NMULTICHOOSEK computes all k-combinations with repetition (multisubsets)
% of a given set.
%
% Let S be the set containing all components of a vector 'vecNum'. This
% function computes all multisubsets of S. Order is not taken into account:
% Two multisubsets, which only differ by a permutation of their elements
% are regarded the same. The output is a matrix whose rows contain all
% k-combinations.
%
% Example: nmultichoose(1:3,2) gives the matrix
%
%       A=[1     1
%          1     2
%          1     3
%          2     2
%          2     3
%          3     3].
%
% Input:
% - vecNum: vector (row or column) whose components are the elements of the
% set of interest.
% - k: parameter speifying the k-combinations.
%
% Output
% - combs: Matrix giving all k-combinations in its rows.

n = numel(vecNum);
combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
combs = reshape(vecNum(combs),[],k);