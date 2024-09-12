function [exponents,coefficients,valid] ...
    = generateRandomCircuit(numVar,deg,maxiter,fullDim)
% GENERATERANDOMCIRCUIT creates a random circuit polynomial for a given
% numer of variables and degree.
%
%   Given a number of variables and a degree, this function creates a 
%   random PSD circuit at the boundary of the SONC cone.The circuit belongs
%   to the boundary of the SONC cone, since it has a zero at the all ones
%   vector. Further, the coefficient of its inner term coincides with the
%   circuit number. Depending on the input parameters, the circuit may have
%   full-dimensional Newton polytope.
%
%   Input:
%   - numVar: number of variables.
%   - deg: degree.
%   - maxiter: maximum number of iterations to compute a suitable Newton
%   polytope (i.e. affinely independent vertices, lattice point in the
%   relative interior, possibly full-dimensional).
%   - fullDim: boolean parameter. Should either be 1 if the circuit
%   should have a full-dimensional Newton polytope, else 0.
%
%   Output:
%   - exponents: exponent matrix, columns contain exponent vectors in the
%   support of the generated circuit.
%   - coefficients: coefficient vector of the generated circuit.
%   - valid: boolean parameter. 1 if circuit has been computed 
%   successfully, else 0.

%% Get random set of vertices
% Get 'numVar'-many random lattice points from the scaled (factor deg/2) 
% standard simplex. As (numVar+1)-th point, add the origin. Compute their
% convex hull and get its set of vertices. If fullDim=1, check whether the
% polytope is full dimensional. Repeat up to 'maxiter' times until an
% admissible set of lattice points is found.

[verticesStandardSimplex,innerTermsStandardSimplex]...
    =latticePointsStandardSimplex(numVar,deg/2);

% Concatenate to a matrix containing all lattice points, except the origin.
verticesPossible...
    =[verticesStandardSimplex(:,1:numVar),innerTermsStandardSimplex];
[~,numVerticesPossible]=size(verticesPossible);

% For convenience, save the sets of indices that do not yield an affine
% independent set of lattice points.
indicesTested=[];

% Boolean parameter, 1 if admissble set of lattice points is found.
valid=0;
% Counter for the numnber of iterations.
iter=1;

while ~valid && iter<=maxiter
    iter=iter+1;
    indicesVertices=sort(randsample(numVerticesPossible,numVar));
    % Add the origin and multiply all other exponents by 2 to have even
    % lattice points.
    verticesInitial...
        =[zeros(numVar,1) 2*verticesPossible(:,indicesVertices)];
    [vertices,~,numVertices] = verticesConvexHull(verticesInitial);
    if size(indicesTested,2)==0 ...
            || ~any(all(indicesTested-indicesVertices==0,1))
        % Proceed only if the set of vertices has not failed before.
        if (~fullDim || numVertices==numVar+1) ...
                && rank(vertices)==(numVertices-1)
            % We wish to find a full dimensional Newton polytope.
            %% Get random inner Term. 
            numInner=1;
            [innerTerm,valid,convexComb]...
                =randomInnerTerm(vertices,[],numInner,1);
            if ~valid
                % There is no lattice point in the relative interior of the
                % convex hull of vertices. Thus 'indicesVertices' fails to
                % deliver an admissible vertex set
                indicesTested=[indicesTested indicesVertices];
            end
        else
            indicesTested=[indicesTested indicesVertices];
        end
    end
end

if valid
    %% Remove all vertices, which are not needed for the convex comb.
    vertexNeeded=~convexComb==0;
    indicesVertices=1:numVertices;
    indicesVertices=indicesVertices(vertexNeeded);
    vertices=vertices(:,indicesVertices);
    % Update number of vertices and get exponent set.
    numVertices=length(vertices);
    exponents=[vertices innerTerm];
    %% Get coefficients.
    % Get coefficients of the vertices.
    coeffVertices=ones(1,numVertices);
    for i=2:numVertices
        coeffVertices(i)=convexComb(i)/convexComb(1);
    end
    % Compute the circuit number
    circuitNumber=1;
    for i=1:numVertices
        circuitNumber=circuitNumber...
            *(coeffVertices(i)/convexComb(i))^convexComb(i);
    end
    %% Get coefficient of the inner term.
    % Coincides with minus the circuit number.
    coeffInnerTerm=-circuitNumber;
    coefficients=[coeffVertices coeffInnerTerm];
    % Multiply all coefficients with a random number
    coefficients=rand(1)*coefficients;
else
    exponents=[];
    coefficients=[];
end

end

