function exponents ...
    = exponentsFromMonomials(monomials, vecVar, addOrigin)
% EXPONENTSFROMMONOMIALS compoutes the exponent matrix corresponding to a
% given set of monomials.
%
%   Given a vector of monomials (can be row or column vector), compute the
%   corresponding matrix of exponents. The i-th column of the matrix
%   'exponents' contains the exponent vector of the i-th monomial in
%   'monomials'. If wished, the origin is added in the end as first
%   exponent vector. If so, all columns of 'exponents' are shifted by one
%   to the right.The variable used must be YALMIP sdpvar decision variables.
%
%   Input:
%   - monomials: vector of monomials.
%   - vecVar: vector of variables of the underlying polynomial ring. Can be
%   row or column vector. 
%   - addOrigin: optional parameter. Valid option: 'addOrigin'. If this is
%   used, the function checks if the origin is already a column of
%   'exponents'. If not, the origin gets added in the first column-.
%
%   Output:
%   - exponents: exponent matrix.

numVar = length(vecVar);
numMon = length(monomials);

% Compute exponent matrix corresponding to monomial vector 'monomials'.
exponents = zeros(numVar, numMon);
for i = 1:numVar
    for j = 1:numMon
        [~,help] = coefficients(monomials(j), vecVar(i));
        exponents(i,j) = degree(help(1));
    end
end

% Add zero vector, if desired and if not already contained.
if nargin==3
    if ischar(addOrigin)
        if strcmp(addOrigin, 'addOrigin')
            % Check if origin is already a column. Therefore check for each
            % column, if the sum of its rows is zero. Ttrick only works for
            % nonnegative exponents. Works for polynomials, not signomials.
            if all(sum(exponents,1)~=0)
                exponents = [zeros(numVar, 1) exponents];
            end
        else
            disp(['Error: Invalid option for ''addOrigin''.',...
                'If you want to add the origin to the set of',...
                'exponents, type ''addOrigin''.']);
        end
    else
        disp(['Error: Invalid data type for ''addOrigin''',...
            '- must be character.']);
    end
end
end
