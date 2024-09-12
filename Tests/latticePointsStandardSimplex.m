function [vertices,innerTerms] = latticePointsStandardSimplex(numVar,deg)
% LATTICEPOINTSSTANDARDSIMPLEX computes all lattice points contained in the
% scaled (factor 'deg') standard 'numVar'-dimensional standard simplex. 
%
%   The function computes an exponent matrix 'vertices' whose columns are
%   the vertices of the scaled standard simples and an exponent matrix
%   'innerTerms' whose columns are all lattice points in the relative
%   interior.
%
%   Input:
%   - numVar: dimensional of the simplex which is to be constructed.
%   - deg: factor of the scaled standard simplex.
%
%   Output:
%   - vertices: exponent matrix containing all vertices of the scaled
%   standard simplex.
%   - innerTerms: exponent matrix containing all lattice points in the
%   interior of the scaled standard simplex.

indicesVar=1:(numVar+1);
% Identify the monomials in 'numVar' variables of degree at most 'deg' with
% multisubsets of 'indicesVar'. A monomial x_1^{a_1}*...*x_n^{a_n} of
% degree k<='deg' is identified with the multisubset containing a_1 times
% 1, a_2 times 2,..., a_n times n.
% Rows of 'helpExponents' represent the multisets.
helpExponents=nmultichoosek(indicesVar,deg);
vertices=[];
innerTerms=[];
for i=1:size(helpExponents,1)
    % Column vector which will be the next exponent.
    exponent=zeros(numVar,1);
    for j=1:deg
        indexVar=helpExponents(i,j);
        if indexVar<=numVar
            % indexVar==numVar+1 is possible to ensure we also have
            % monomials of smaller degree than 'deg'.
            exponent(indexVar)=exponent(indexVar)+1;
        end
    end
    % Check if 'exponent' is vertex or lattice point in the interior.
    if any(exponent==deg) || all(exponent==0)
        vertices=[vertices exponent];
    else
        innerTerms=[innerTerms exponent];
    end
end
vertices=round(vertices);
innerTerms=round(innerTerms);
end