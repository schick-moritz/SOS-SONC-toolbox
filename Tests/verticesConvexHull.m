function [vertices,innerTerms,numVertices]...
    = verticesConvexHull(exponentsGiven)
% VERTICESCONVEXHULL computes all sets of vertices and lattice points in
% the relatve interior of the convex hull of a given set of exponents.
%
%   This function considers the polytope determined by the convex hull of a
%   given set 'exponentsGiven' of exponents. It computes the the set of
%   vertices and lattice points in the relative interior of this polytope.
%
%   Input:
%   - exponentsGiven: given set of lattice points. 
%
%   Output:
%   - vertices: matrix describing the vertex set of the polytope. Vertices 
%   are columns.
%   - innerTerms: matrix containing the lattice points in the relative 
%   interior of the polytope in its columns.
%   - numVertices: the number of vertices.

%% From the given exponents, check which are vertices.
vertices = exponentsGiven;
innerTerms = [];
solverOptions = sdpsettings('solver', '', 'verbose', 0);
if size(exponentsGiven,2) > 1
    % Only compute the convex hull if there are at least two given
    % exponents, else it is trivial.
    i=1;
    while i<=size(vertices,2)
        % Stop if there is only one vertex left.
        exponentsConvexHull = vertices(:,[1:i-1 i+1:size(vertices,2)]);
        % Check if exponentsConvexHull contain only the origin only. Note:
        % only works for matrices with nonnegative integer components.
        if all(sum(exponentsConvexHull)==0)
            i=i+1;
        else
            %% Optimization
            coeffConvexComb = sdpvar(size(vertices,2)-1,1);
            % Test if the i-th column of 'vertices' is a vertex.
            exponentTest = vertices(:,i);
            % Assemble the optimization problem determining if we have a
            % vertex.
            constraints...
                =[exponentsConvexHull*coeffConvexComb==exponentTest,...
                sum(coeffConvexComb)==1, coeffConvexComb>=0];
            objective=0;
            diagnostics=optimize(constraints, objective, solverOptions);
            if diagnostics.problem==0
                % Feasible: 'testExponent is no vertex'. Remove it from the
                % matrix containing the vertices and add it to the matrix
                % containing the inner Terms.
                vertices=exponentsConvexHull;
                innerTerms=[innerTerms exponentTest];
            elseif diagnostics.problem == 1
                % Infeasible: testExponent is vertex
                i=i+1;
            else
                disp(['Error occurred while attempting to solve',...
                    'optimization problem in verticesConvexHull.m']);
            end
        end
    end
end
numVertices=size(vertices,2);
end