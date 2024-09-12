function [cellExpSONC,cellCoeffSONC,cellPolyValidSONC] ...
    = generateRandomSONC(vecNumVar,vecDeg,vecNumCirc,fullDim)
% GENERATERANDOMESONC creates a collection of random SOS polynomials.
%
%   This function creates a collection of random SONC polynomials for
%   given number of variables, degree, and number of PSD circuits in an 
%   SONC decomposition. 
%
%   Input:
%   - vecNumVar: vector containing the numbers of variables of SOS
%   polynomials which are to be constructed. Can be row or column vector.
%   - vecDeg: vector containing the degrees SOS polynomials which are to be
%   constructed. Can be row or column vector.
%   - vecNumSquares: vector containing number of squares of SOS polynomials
%   which are to be constructed. Can be row or column vector.
%
%   Output:
%   - cellExpSONC: exponent matrices of random SOS polynomials.
%   - cellCoeffSONC: coefficient vectors of random SOS polynomials.
% This function generates a collection of random SONC polynomials for
% different instances of number of variables, degree and number of
% circuits. In general, the polynomials constructed are non-homogeneous.
% Input:
%   - vecNumVar: vector containing number of variables of SONCs which are
%   constructed.
%   - vecDeg: vector containing degrees of SONCs which are constructed.
%   - vecNumCirc: vector containing number of circuits of SONCs which are
%   constructed.
%   - fullDim: optional boolean parameter. 1 if circuits shall have full 
%   dimensional Newton polytopes. Else 0 or leave empty.
% Output:
%   - cellExpSONC: exponent matrices of random SONCs.
%   - cellCoeffSONC: coefficient vectors of random SONCs.
%   - cellPolyValidSONC: cell containing boolean parameter for every random
%   SONC to be constructed. 1 if construction was successful, else 0.

% Get dimensions of cells, which are to be constructed.
sizeNumVar=length(vecNumVar);
sizeDeg=length(vecDeg);
sizeNumCirc=length(vecNumCirc);

% Vector of variables. Needed to get exponent matrix and coefficient
% vectors of SONCs with at least two circuits.
variables=sdpvar(max(vecNumVar),1);

% Initialize the output parameters describing the random polynomials
cellExpSONC=cell(sizeNumVar,sizeDeg,sizeNumCirc); 
cellCoeffSONC=cell(sizeNumVar,sizeDeg,sizeNumCirc);
cellPolyValidSONC=cell(sizeNumVar,sizeDeg,sizeNumCirc);

%% Get Random SONC polynomials
for i=1:sizeNumVar
    numVar=vecNumVar(i);
    for j=1:sizeDeg
        deg=vecDeg(j);
        for k=1:sizeNumCirc
            numCirc=vecNumCirc(k);
            % Get the actual random polynomial
            poly=0;
            validSONC=0;
            for l=1:numCirc
                [expCircuit,coeffCircuit,valid] ...
                    = generateRandomCircuit(numVar,deg,1000,fullDim);
                if valid
                    validSONC=1;
                    for r=1:size(expCircuit,2)
                        mon=coeffCircuit(r);
                        for s=1:numVar
                            mon=mon*variables(s)^expCircuit(s,r);
                        end
                        poly=poly+mon;
                    end
                end
            end
            cellPolyValidSONC{i,j,k}=validSONC;
            [cellCoeffSONC{i,j,k},monomialsPoly]...
                =coefficients(poly,variables(1:numVar));
            cellExpSONC{i,j,k} ...
                =exponentsFromMonomials(monomialsPoly,variables(1:numVar));
        end
    end
end
end

