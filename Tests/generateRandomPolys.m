function [cellExp,cellCoeff,cellPolyValid] ...
    = generateRandomPolys(vecNumVar,vecDeg,vecNumMon,testScenario)
% GENERATERANDOMPOLYS creates a collection of random polynomials
%
%   This function generates a collection of random polynomials for a given
%   number of variables, degree and number of monomials. Different options
%   regarding the specific Newtn polytope structure are available. In 
%   general, the polynomials constructed are not homogeneous.
%
%   Input:
%   - vecNumVar: vector containing the numbers of variables of polynomials 
%   which are constructed.
%   - vecDeg: vector containing degrees of polynomials which are 
%   constructed.
%   - vecNumMon: vector containing the numbers monomials of polynomials 
%   which are constructed.
%   - testScenario: integer parameter specifying the Newtn polytope
%   configuration of the polynomials which are constructed
%   	1 for standard simplex Newton polytope,
%   	2 for random simplex Newton polytope,
%   	3 for arbitrary Newton polytope.
% Output:
%   - cellExp: exponent matrices of random polynomials.
%   - cellCoeff: coefficient vectors of random polynomials.
%   - cellPolyValid: cell containing boolean parameter for every random
%   polynomial. 1 if construction was successful, else 0.

% Get dimensions.
sizeNumVar=length(vecNumVar);
sizeDeg=length(vecDeg);
sizeNumMon=length(vecNumMon);

% Initialize the output parameters describing the random polynomials.
cellExp=cell(sizeNumVar,sizeDeg,sizeNumMon); 
cellCoeff=cell(sizeNumVar,sizeDeg,sizeNumMon); 
cellPolyValid=cell(sizeNumVar,sizeDeg,sizeNumMon);

for i=1:sizeNumVar
    numVar = vecNumVar(i);
    for j=1:sizeDeg
        deg = vecDeg(j);
        for k=1:sizeNumMon
            numMon = vecNumMon(k);
            if testScenario==1
                %% Standard simplex Newton polytopes
                % Step 1: construct exponent matrix. First numVar+1
                % columns are vertices of the standard simplex, remaining 
                % numMon-(numVar+1) (if >=0) columns are random inner 
                % lattice points.
                [vertices,innerTermsStandardSimplex]...
                    =latticePointsStandardSimplex(numVar,deg);
                numVertices=numVar+1;
                if numMon-numVertices>=0
                    cellPolyValid{i,j,k}=1;
                    % Get random inner terms.
                    indicesInnerTerms...
                        =randsample(size(innerTermsStandardSimplex,2),...
                        numMon-numVertices);
                    innerTerms=innerTermsStandardSimplex(:,indicesInnerTerms);
                end
            elseif testScenario==2 || testScenario==3
                maxiter=200;
                if testScenario==2
                    %% Random simplex Newton polytope.
                    simplex=1;
                    fullDim=1;
                    % Needed as input for the function
                    % latticePointsRandomNewtonPolytopeOptions
                    parameterInner=[];
                else
                    %% Arbitrary Newton polytope.
                    simplex=0;
                    % When constructing arbitrary Newton polytopes, we want
                    % a lower bound on the number of innerTerms.
                    parameterInner=3; 
                    fullDim=1;
                end
                % Compute random exponent matrix.
                [cellPolyValid{i,j,k},vertices,innerTerms,numVertices]...
                    =latticePointsRandomNewPol(numVar,deg,numMon,...
                    maxiter,simplex,parameterInner,fullDim);
            end
            if cellPolyValid{i,j,k}
                % Exponents have been constructed successfully.
                exponents=[vertices innerTerms];
                cellExp{i,j,k}=exponents;
                % Compute random coefficient vector.
                coeffVertices...
                    =abs(normrnd(0,(numMon/numVar)^2,[1,numVertices]));
                coeffInnerMon=normrnd(0,1,[1,numMon-numVertices]);
                coeff=[coeffVertices,coeffInnerMon];
                cellCoeff{i,j,k}=coeff;
            end
        end
    end
end
end

