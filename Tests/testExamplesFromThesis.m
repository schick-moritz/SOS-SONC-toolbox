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

sdpvar x y z;
vecTwoVar=[x y];
vecThreeVar=[x y z];
test={};
i=1;
%% PSD/ SONC polynomial that are not SOS or not SOS+SONC.
% 1: Motzkin polynomial.
M=x^4*y^2+x^2*y^4+1-3*x^2*y^2; 
test{i,1}=M; 
test{i,2}=vecTwoVar;
i=i+1;
% 2: Choi-Lam polynomial.
Q=x^2*y^2+x^2*z^2+y^2*z^2+1-4*x*y*z; 
test{i,1}=M; 
test{i,2}=vecThreeVar;
i=i+1;
% 3: Second Choi-Lam polynomial.
S=x^4*y^2+y^4+x^2-3*x^2*y^2; 
test{i,1}=S; 
test{i,2}=vecTwoVar;
i=i+1;
% 4: Modified Motzkin polynomial.
Mhat=1-x^2*y^2+x^4*y^2+x^2*y^4; 
test{i,1}=Mhat; 
test{i,2}=vecTwoVar;
i=i+1;
% 5: Robinson polynomial.
R=x^6+y^6+1-(x^4*y^2+x^4+x^2*y^4+x^2+y^4+y^2)+3*x^2*y^2; 
test{i,1}=R; 
test{i,2}=vecTwoVar;
i=i+1;
% 6: Second Robinson polynomial.
Rhat=x^2*(x-1)^2+y^2*(y-1)^2+z^2*(z-1)^2+2*x*y*z*(x+y+z-2); 
test{i,1}=Rhat; 
test{i,2}=vecThreeVar;
i=i+1;
% 7: Schmüdgen polynomial.
P=200*((x^3-4*x)^2+(y^3-4*y)^2)+(y^2-x^2)*x*(x+2)*(x*(x-2)+2*(y^2-4)); 
test{i,1}=P; 
test{i,2}=vecTwoVar;
i=i+1;

%% SOS polynomials that are not SONC.
% 8: Iliman and de Wolff Example.
test{i,1}=(x-1)^2*(x-2)^2; 
test{i,2}=x;
i=i+1;
% 9: Dressler Example.
test{i,1}=(x+y+1)^2; 
test{i,2}=vecTwoVar;
i=i+1;
% 10: Dressler Example.
test{i,1}=(x^2+2*x+1)^2; 
test{i,2}=x;
i=i+1;
% 11: Dressler Example.
test{i,1}=(x*y+x+y)^2; 
test{i,2}=vecTwoVar;
i=i+1;
% 12: Dressler Example.
test{i,1}=(x+1)^4; 
test{i,2}=x;
i=i+1;

%% Powers of PSD/ circuit polynomials.
% 13: Stengle polynomial.
T=x^3+(y^2-x^3-x)^2;
test{i,1}=T; 
test{i,2}=vecTwoVar;
i=i+1;
% 14: Motzkin square.
test{i,1}=M^2; 
test{i,2}=vecTwoVar;
i=i+1;
% 15: Motzkin cube.
test{i,1}=M^3; 
test{i,2}=vecTwoVar;
i=i+1;
% 16: Choi-Lam square.
test{i,1}=S^2; 
test{i,2}=vecTwoVar;
i=i+1;
% 17: Choi-Lam cube.
test{i,1}=S^3; 
test{i,2}=vecTwoVar;
i=i+1;

%% Intersection of SOS and SONC cones. Polynomials are no SOBS.
% 18: Motzkin plus SOBS.
test{i,1}=M+(x^2-1)^2+(y^2-1)^2; 
test{i,2}=vecTwoVar;
i=i+1;
% 19: Choi-Lam plus SOBS.
test{i,1}=Q+(1-x^2)^2+(1-y^2)^2+(1-z^2)^2; 
test{i,2}=vecThreeVar;
i=i+1;

%% The SOS+SONC cone is not closed under multiplication.
% 20: Quadratic times Motzkin.
test{i,1}=(x+y+1)^2*M; 
test{i,2}=vecTwoVar;
i=i+1;
% 21: Quadratic times Choi-Lam.
test{i,1}=(x+y+z+1)^2*Q; 
test{i,2}=vecThreeVar;
i=i+1;

%% The SOS+SONC cone is not closed under affine variable transformations.
% 22: Motzkin variable transformation.
test{i,1}=replace(M,[x y],[x-1 y-1]); 
test{i,2}=vecTwoVar;
i=i+1;
% 23: Choi-Lam variable transformation.
test{i,1}=replace(Q,[x y z],[x-1 y-1 z-1]); 
test{i,2}=vecThreeVar;
i=i+1;

%% SOS+SONC polynomials that are neither SOS not SONC.
% 24: SOS plus Motzkin.
test{i,1}=1/2*(1+2*x*y+x^2*y)^2+M; 
test{i,2}=vecTwoVar;
i=i+1;
% 25: SOS plus Choi-Lam.
test{i,1}=(x*y+x*z+y*z)^2+1+Q; 
test{i,2}=vecThreeVar;
i=i+1;

%% SOS+SONC polynmials with no support-invariant SOS+SONC decomposition.
% 26: SOS plus Motzkin;
test{i,1}=2*(x*y^2+x*y)^2+(x^2*y+2*x*y+1)^2+2*M; 
test{i,2}=vecTwoVar;
i=i+1;
% 27: SOS plus Choi-Lam.
test{i,1}=2*(x*y+x+y*z)^2+Q; 
test{i,2}=vecThreeVar;
i=i+1;

%% SOS+SONC polynomials with several SOS+SONC decompositions.
% 28: SOS polynomial, not SONC.
test{i,1}=3*(x*y^2+y)^2+3*(x^2*y+x)^2+6*(x^2*y+x*y^2)^2+x^4*y^2+x^2*y^4+4; 
test{i,2}=vecTwoVar;
i=i+1;
% 29: Both SOS polynomial and SONC polynomial.
test{i,1}=2*(x*y+y*z)^2+(y*z+z)^2+(x*y+x)^2+x^2*z^2+1 ; 
test{i,2}=vecThreeVar;
i=i+1;
% 30: SONC polynomial
test{i,1}=(x*y^2+y+1)^2+(x^2*y+x+1)^2+6*(x^2*y-x)^2+6*(x*y^2-y)^2+...
    2*(x^4*y^2+x^2*y^4+4-2*x*y); 
test{i,2}=vecTwoVar;
i=i+1;
% 31: SOS polynomial, not SONC.
test{i,1}=(x^2*y+x+y+2)^2+(x*y^2+x+y+2)^2+2*(x^2*y-y)^2+2*x^2*y^4+2*x^2+6; 
test{i,2}=vecTwoVar;
% i=i+1;

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
data=[1:numPolynomials;problem';optVal';runTime']';
% Symbolic array.
dataSym=[sym(1:numPolynomials);sym(problem)';...
    sym(compose('%8.4g',optVal))';sym(compose('%8.4g',runTime))']';

% Get LaTex code
latexData=latex(vpa(dataSym,3));
 
% Save the data
save(strcat("data/polynomialsThesis ",...
    string(datetime(datetime,'InputFormat',...
    'yyyy-MM-dd HH:mm:ss.SSS'))));

