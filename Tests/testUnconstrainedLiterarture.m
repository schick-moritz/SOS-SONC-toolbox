close all; clc; clear;
yalmip('clear');

disp(today('datetime'));
fprintf('Moritz Schick, University of Konstanz\n');
fprintf('Testing membership in the SOS, SONC and SOS+SONC cones \n');

%% Summary
% In this MATLAB script, the SOS, SONC and SOS+SONC relaxations are applied
% to a list of explicit polynomials that can be found throughout the
% thesis.

%% Explanation
% YALMIP decision variables must be used.
% We construct a cell 'test' containing the dara describing the polynomials
% which are to be tested. The first column of 'test' contains the
% polynomials and the second column the corresponding row vectors of
% variables. For example test{1,1} is the Motzkin polynomial and test{1,2}
% the variable vector [x y].

sdpvar x y;
vecTwoVar=[x y];
test={};
i=1;

%% Literature: 
% Dressler, Iliman, de Wolff: An approach to constrained polynomial 
% optimization via nonnegative circuit polynomials and geometric 
% programming (2019).

% Example 5.4.
% Optimal values in the literature: SONC: 3.269 or 3.572 (depending on 
% triangulation). SOS: 3.8673.
test{i,1}=6+x^2*y^6+2*x^4*y^6+x^8*y^2-1.2*x^2*y^3-...
    0.85*x^3*y^5-0.9*x^4*y^3-0.73*x^5*y^2-1.14*x^7*y^2;
test{i,2}=vecTwoVar;
i=i+1;

% Example 5.5.
% Optimal values in the literature: SONC: 0.5732 or 0.6583 (depending on 
% distribution of coefficients). SOS: 0.8383. 
% In Seidler, de Wolff: SONC: 0.693158.
test{i,1}=1+3*x^2*y^6+2*x^6*y^2+6*x^2*y^2-x*y^2-2*x^2*y-3*x^3*y^3;
test{i,2}=vecTwoVar;
i=i+1;

% Example 5.6.
% Optimal values in the literature: None. (GP approach not applicable).
test{i,1}=1+x^4+y^2+x^2*y^4+x^4*y^4-x*y-x*y^2-x^2*y^3-x^3*y^3;
test{i,2}=vecTwoVar;

%% Run the test
[numPolynomials,~]=size(test);
numRelaxations=3; 

optVal=-42*ones(numPolynomials,numRelaxations);
problem=-42*ones(numPolynomials,numRelaxations);
runTime=-42*ones(numPolynomials,numRelaxations);

for i=1:numPolynomials
    fprintf(['\n', 'Polynomial number ', num2str(i), '\n']);
    f=test{i,1};
    vecVar=test{i,2};
    %% SOS relaxation.
    tic;
    [optVal(i,1),problem(i,1)]=globalMinSOS(f,vecVar);
    runTime(i,1)=toc;
    
    %% SONC relaxation.
    tic;
    [optVal(i,2),problem(i,2)]=globalMinSONC(f,vecVar);
    runTime(i,2)=toc;

    %% SOS+SONC relaxation.
    tic;
    [optVal(i,3),problem(i,3)]=globalMinSOSpSONC(f,vecVar);
    runTime(i,3)=toc;
end

% Collect all the test data.
% Without certificates.
data=[1:numPolynomials;problem';optVal';runTime']';
% Symbolic array.
dataSym=[sym(1:numPolynomials);sym(problem)';...
    sym(compose('%8.4f',optVal))';sym(compose('%8.2f',runTime))']';

% Get LaTex code
latexData=latex(dataSym);
 
% Save the data
save(strcat("data/unconstrainedLiterature ",...
    string(datetime(datetime,'InputFormat',...
    'yyyy-MM-dd HH:mm:ss.SSS'))));

