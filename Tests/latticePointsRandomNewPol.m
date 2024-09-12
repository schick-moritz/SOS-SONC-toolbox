function [valid,vertices,innerTerms,numVertices]...
    = latticePointsRandomNewPol(numVar,deg,numMon,maxiter,simplex,...
    parameterInner,fullDim)
% LATTICEPOINTSRANDOMNEWPOL creates a random set of vertices and inner
% terms of a polynomial.
%
%   This function computes the vertex set of a random polynomial and a
%   number of inner lattice points determined by the input parameters.
%   Additionally, the Newton polytope may have a particular structure
%   defined by the input parameters. It could be a simplex and/or full 
%   dimensional.
%
%   Input:
%   - numVar: number of variables of the polynomial the Newton polytope
%   corresponds to. This is the dimension of the underlying vector space.
%   - deg: degree of the polynomial the Newton polytope corresponds to.
%   - numMon: number of monomials of the polynomial. This is the total
%   number of lattice points which are computed (number of vertices plus
%   number of inner terms).
%   - maxiter: maximum number of iterations for computing an admissible
%   polytope.
%   - simplex:boolean parameter. 1 if the Newton polytope shall be a 
%   simplex, else 0.
%   - parameterInner: ineteger parameter for number of inner terms. Only
%   needed if simplex=0.
%   - fullDim: boolean parameter. 1 if the Newton polytope shall be
%   full-dimensional, else 0.
%
%   Output:
%   - valid: boolean parameter. 1 if construction was successful, else 0.
%   - vertices: vertices of the constructed Newton polytope.
%   - innerTerms: inner terms of the constructed Newton polytope.
%   - numVertices: number of vertices.

if simplex
    numInner=numMon-numVar-1;
else
    numInner=ceil(parameterInner/5*(numMon-numVar-1));
end

if numMon-numInner>=2
    % Get a random set of vertices from the scaled (factor deg/2) standard
    % simplex. Make sure that the origin is one of the vertices. Multiply 
    % them  by two to make sure they are even lattice points. Check if they 
    % satisfy the requirements specified by the input parameters.
    [verticeScaledStandardSimplex,innerTermsScaledStandardSimplex]...
        =latticePointsStandardSimplex(numVar,deg/2);
    
    % Exlude the origin from the set of vertices, from which we choose
    % numMon-numInner-1 many ones.
    verticesAll...
        =[verticeScaledStandardSimplex(:,1:numVar) ...
        innerTermsScaledStandardSimplex];
    
    % Boolean parameter, 1 if admissble set of lattice points is found.
    valid=0;
    iter=1;
    
    while ~valid && iter<=maxiter
        iter=iter+1;
        indicesVertices=randsample(size(verticesAll,2),numMon-numInner-1);
        verticesTested...
            =[zeros(numVar,1) 2*verticesAll(:,indicesVertices)];
        [vertices,innerTerms,numVertices]...
            =verticesConvexHull(verticesTested);
        numInnerUpdate=numMon-numVertices;
        % First check, if set of vertices satisfies the requirements. Then
        % attempt to find numMon-numInner many inner terms.
        tryToFindInnerTerms=0;
        if simplex
            if rank(vertices)==(numVertices-1) 
                if numVertices==numVar+1 || ~fullDim
                    % Newton polytope shall be full-dimensional
                    tryToFindInnerTerms=1;
                end
            end
        else
            if rank(vertices)==numVar || ~fullDim
                tryToFindInnerTerms=1;
            end
        end
        if tryToFindInnerTerms
            % Change the last parameter in randomInnerTermsNew to 1, if the
            % inner terms should be in the relative interior of the
            % Newton polytope.
            [innerTerms,validInnerTerms,~]...
                =randomInnerTerm(vertices,[],numInnerUpdate,0);
            if validInnerTerms
                valid=1;
            end
        end
    end
    vertices=round(vertices);
    innerTerms=round(innerTerms);
else
    valid=0;
    vertices=[];
    innerTerms=[];
    numVertices=[];
end
end