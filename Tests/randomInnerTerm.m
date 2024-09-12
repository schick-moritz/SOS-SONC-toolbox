function [innerTerms,valid,convexComb] ...
    = randomInnerTerm(vertices,innerTerms,numInner,properInner)
% RANDOMINNERTERM computes a lattices point in the relative interior of a
% given polytope (or in one of its faces).
%
%   Consider a polytope, which is given by its vertex set 'vertices'. This
%   function computes a lattice point in the relative interior of the
%   polytope. Some optional input parameters can be used. If a set of
%   relative interior points is given in 'innerTerms', the function makes
%   sure that the new point does not coincide with any of the given points.
%   By 'numInner', you can tell the function how many inner terms should be
%   computed. The parameter 'properInner' can be used if the new inner term
%   should be in the relative interior of the polytope.
% 
% Input:
%   - vertices: matrix containing the vertices of the given polytope in its
%   columns.
%   - innerTerms: matrix containing given inner terms in its columns. The
%   function computes a new inner term, which is distinct from these.
%   - numInner: number of inner terms which are to be computed. 
%   - relInt: optional boolean parameter. 1 if the computed lattice points
%   should be in the relative interior of the polytope, 0 if they are also 
%   allowed to be in the relative interior of on of its faces.
%
% Output:
%   - innerTerms: concatenation of the input parameter 'innerTerms' and the 
%   'numInner'-many newly computed inner terms.
%   - valid: boolean parameter. 1 if inner terms have been computed 
%   successfully, else 0.
%   - convexComb: matrix containing coefficient vectors of the convex
%   combinations of the newly computed inner terms in its rows. 

[numVar,numVertices]=size(vertices);

valid=1;
convexComb=[];

for i=1:numInner 
    %% Optimization
    % Solve an optimization problem to find an inner term.
    exponentTest=sdpvar(numVar,1); % This will be the new inner Term.
    coeffConvexComb=sdpvar(numVertices,1);
    % This parameter is used to make sure that the new inner term is not a
    % vertex and in the relative interior, if necessary.
    posConstr=double(lcm(sym(1:max(vertices,[],'all'))));
    constraints=[vertices*coeffConvexComb==exponentTest, ...
            integer(exponentTest), sum(coeffConvexComb)==1,];
    if nargin==3 || (nargin==4 && ~properInner)
        % New inner term must not be in the relative interior but cannot be
        % a vertex. 
        constraints=[constraints, coeffConvexComb>=0, ...
            coeffConvexComb<=(1-1/posConstr)];
    elseif nargin==4 && properInner
        % Inner term must be in the relative interior and cannot be a
        % vertex. Relative interior shows that the convex combination
        % coefficients must be strictly positive. However, as strict
        % inequalities do not make sense in a numerical program, we bound
        % them from below by 1/posConstr.
        constraints = [constraints, coeffConvexComb>=1/posConstr];
    end
    % If there are already inner terms, make sure we compute a new one.
    if ~size(innerTerms,2)==0
        constraints=[constraints,sum(...
            abs(innerTerms-exponentTest*ones(1,size(innerTerms,2))))>=1];
    end
    % Random cost functional in order to get random inner terms.
    objective = (rand(numVar,1)-0.5)'*exponentTest;
    solverOptions = sdpsettings('solver', '', 'verbose', 0);
    diagnostics = optimize(constraints, objective, solverOptions);
    if diagnostics.problem == 0
        % New inner term was found.
        innerTerms=[innerTerms value(exponentTest)];
        convexComb=[convexComb value(coeffConvexComb)];
    else
        valid=0;
        break;
    end
end

