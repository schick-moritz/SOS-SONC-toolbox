function [coeffSigRep, mon, exponents, indInnerTerms] = sigRep(f, vecVar)
% SIGREP computes the signomial representative of a given polynomial.
%
%   Given a polynomial 'f' defined in the variables 'vecVar', compute its
%   signomial representative. The signomial representative is then given by
%   its coefficients 'coeffSigRep', exponent matrix 'exp' and the indices
%   of its inner terms, corresponding to columns of 'exp'. Inner terms are
%   all indices of columns of 'exp', which are not even lattice points. We
%   also give all monomials 'mon' of 'f'. The variables used must be YALMIP
%   sdpvar decision variables.
%
%   Input:
%   f: the polynomial of interest.
%   vecVar: vector of variables, in which f is defined. Can be row or
%   column vector.
%
%   Output
%   - coeffSigRep: coefficient of the signomial representative of 'f',
%   column vector.
%   - mon: monomials of 'f'.
%   - exponents: exponent matrix of 'f' respectively its signomial
%   representative.
%   - indInnerTerms: vector of indices corresponding to columns of 'f',
%   which are not even lattice points.

% Get coefficients and monomial of f.
[coeffF, mon] = coefficients(f, vecVar);
numMon = length(mon);

% Get exponent matrix corresponding to monomials 'mon' and variables
% 'vecVar'.
exponents = exponentsFromMonomials(mon,vecVar);

% Identify exponents with at least one odd component.
helpExponentsOdd = mod(exponents,2);
if size(exponents,1)>=2
    % If there are at least two variables, check all the exponents for odd
    % number.
    helpExponentsOdd = any(helpExponentsOdd);
else
    helpExponentsOdd = logical(helpExponentsOdd);
end
allIndices = 1:numMon;
indInnerTerms = allIndices(helpExponentsOdd);

% Coefficient vector of signomial representative: Start with coefficient
% vector of 'f'
coeffSigRep = coeffF;
for i=1:length(indInnerTerms)
    help = coeffSigRep(indInnerTerms(i));
    coeffSigRep(indInnerTerms(i)) = -sign(help)*help;
end
