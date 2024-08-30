function [monHalfNewF, exponentsHalfNewF]...
    = halfNewtonPolytope(f, vecVar, exponentsNewF)
% HALFNEWTONPOLYTOPE computes the lattice points in 1/2*New(f), where
% New(f) is the Newton polytope of a given polynomial 'f'.
%
%   Given a polynomial 'f', compute the exponent matrix 'expHalfNewF',
%   whose columns are all lattice point in 1/2 times the Newton polytope of
%   'f'. Further, compute the vector 'monHalfNewF' of monomials with
%   coefficient 1 such that the exponent of the i-th component of
%   'monHalfNewF' is the i-th column of 'expHalfNewF'. The variables used
%   must be YALMIP sdpvar decision variables.
%
%   Input:
%   - f: the given polynomial.
%   - vecVar: vector of variables in which 'f' is defined. Can be row or
%   column vector.
%   - exponentsNewF: optional parameter. If this is used, it should contain
%   all lattice points in the Newton polytope of f. Makes sense if the 
%   function 'exponentsFromMonomials(f, vecVar)' is used in an earlier 
%   computation.
%
%   Output:
%   - monNewF: column vector containing all monomials in 1/2 times the
%   Newton polytope of 'f'.
%   - expNewF: exponent matrix, columns are lattice points in 1/2 times the
%   Newton polytope of 'f'.


if nargin==2
    % Exponents matrix of lattice points in Newton polytope of 'f' is not
    % given as an input parameter.
    [~, exponentsNewF]=newtonPolytope(f, vecVar);
end

% First, multiply exponent matrix by 0.5. Then check, which columns have
% only integer components (i.e. are lattice points) and remove all other
% columns.
exponentsHalfNewF=0.5*exponentsNewF;
[~,numMon]=size(exponentsHalfNewF);
exponentIsLatticePoint=zeros(1,numMon);
for j=1:numMon
    if all(exponentsHalfNewF(:,j)==floor(exponentsHalfNewF(:,j)))
        exponentIsLatticePoint(1,j)=1;
    end
end

% Get exponent matrix containing all lattice points in 1/2 times the Newton
% polytope of 'f'.
exponentsHalfNewF=exponentsHalfNewF(:,logical(exponentIsLatticePoint));

% Get row vector of monomials corresponding to the columns of expHalfNewF.
numLatticePointsHalfNew=sum(exponentIsLatticePoint);
monHalfNewF=sdpvar(numLatticePointsHalfNew,1);
for i=1:numLatticePointsHalfNew
    monHalfNewF(i,1)=prod(vecVar'.^exponentsHalfNewF(:,i));
end
end
