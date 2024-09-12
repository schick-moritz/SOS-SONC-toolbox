function poly = polynomialFromExpCoeffVar(exponents,coefficients,variables)
% POLYNOMIALFROMEXPCOEFFVAR computes a YALMIP polynomial determined by the
% input parameters,
%   Input:
%   - exponents: exponent matrix of the polynomial. Exponents should be
%   columns.
%   - coefficients: coefficient vector of the polynomial.
%   - variables: vector of variables. Could be row or column vector.
%   Output:
%   - poly: the polynomial.

poly=0;
numVar=length(variables);
for r=1:size(exponents,2)
    mon=coefficients(r);
    for s=1:numVar
        mon=mon*variables(s)^exponents(s,r);
    end
    poly=poly+mon;
end
end

