function [cellExpSOS,cellCoeffSOS] ...
    = generateRandomSOS(vecNumVar,vecDeg,vecNumSquares)
% GENERATERANDOMESOS creates a collection of random SOS polynomials.
%
%   This function creates a collection of random SOS polynomials for
%   given number of variables, degree, and number of squares in an SOS
%   decomposition. The constructed polynomials all belong to the boundary
%   of the SOS cone, since they have a zero at the all ones vector.
%
%   Input:
%   - vecNumVar: vector containing the numbers of variables of SOS
%   polynomials which are constructed. Can be row or column vector.
%   - vecDeg: vector containing the degrees SOS polynomials which are 
%   constructed. Can be row or column vector.
%   - vecNumSquares: vector containing number of squares of SOS polynomials
%   which are constructed. Can be row or column vector.
%
%   Output:
%   - cellExpSONC: exponent matrices of random SOS polynomials.
%   - cellCoeffSONC: coefficient vectors of random SOS polynomials.

% Get dimensions of the cells, which are to be constructed.
sizeNumVar=length(vecNumVar);
sizeDeg=length(vecDeg);
sizeNumSquares=length(vecNumSquares);

% Initialize a vector of variables.
variables=sdpvar(max(vecNumVar),1);

% Initialize the output parameters describing the random polynomials.
% cellGramMatrix=cell(sizeNumVar,sizeDeg);
cellCoeffSOS=cell(sizeNumVar,sizeDeg);
cellExpSOS=cell(sizeNumVar,sizeDeg);

%% Generate random SOS polynomials.
for i=1:sizeNumVar
    numVar=vecNumVar(i);
    for j=1:sizeDeg
        deg=vecDeg(j);
        %% Generate random half degree polynomials for SOS decomposition.
        for k=1:sizeNumSquares
            numSquares=vecNumSquares(k);
            % Number of monomials of each square.
            numMon=min([(numVar*(numVar+1))/2 (deg/2*(deg/2+1))/2]);
            % Get exponent matrix of all half degree monomials.
            [verticesHalfNewton,innerTermsHalfNewton]...
                =latticePointsStandardSimplex(numVar,deg/2);
            exponentsHalfNewton...
                =[verticesHalfNewton(:,1:numVar) innerTermsHalfNewton];
            % Initialize the SOS polynomial, and add numSquares-many
            % random squares.
            poly=0;
            for l=1:numSquares
                % Get random set of exponents.
                indicesExp=...
                    randsample(size(exponentsHalfNewton,2),numMon-1);
                exponents=[zeros(numVar,1) exponentsHalfNewton(:,indicesExp)];
                % Get random coefficient vector. Make sure that the
                % coefficients sum up to zero, such that the polynomial has
                % a zero at the all ones vector.
                coeffHelp=-1+2*rand(1,numMon-1);
                signHelp=sign(sum(coeffHelp));
                coeff=signHelp*[sum(coeffHelp) -coeffHelp];
                polyHalfDeg=0;
                for r=1:size(exponents,2)
                    mon=coeff(r);
                    for s=1:numVar
                        mon=mon*variables(s)^exponents(s,r);
                    end
                    polyHalfDeg=polyHalfDeg+mon;
                end
                poly=poly+polyHalfDeg^2;
            end
            % Recover the coefficient vector and exponent set of the
            % generated random SOS polynomial.
            [cellCoeffSOS{i,j,k},v]=coefficients(poly);
            cellExpSOS{i,j,k}...
                =exponentsFromMonomials(v,variables(1:vecNumVar(i)));
        end
    end
end
end

