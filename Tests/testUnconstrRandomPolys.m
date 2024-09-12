close all; clc; clear;
yalmip('clear');

disp(today('datetime'));
fprintf('Moritz Schick, University of Konstanz\n');
fprintf(['Unconstrained polynomial optimization and \n',...
    'testing membership in the SOS, SONC and SOS+SONC cones \n']);

%% Summary
% In this script, we first construct a collection of random polynomials
% with either (1) standard simplex, (2) simplex or (3) arbitrary Newton
% polytopes. We then test our MATLAB functions for computing lower bounds
% on their global minima using an SOS, SONC, SOS+SONC and - if desired - an
% SOS-before-SONC or SONC-before-SOS relaxation.

% Numbers of variables tested.
vecNumVar=[2 3 4];

% Degrees tested.
vecDeg=[6 8];

% Number of monomials tested.
vecNumMon=12;

sizeNumVar=length(vecNumVar);
sizeDeg=length(vecDeg);
sizeNumMon=length(vecNumMon);

% Choose the type of Newton polytopes you wish to be tested.
testScenario=1; % Standard Simplex Newton polytopes.
% testScenario=2; % Simplex Newton polytopes.
% testScenario=3; % Arbitrary Newton polytopes.

%% Get exponents and coefficients of random polynomials.
[cellExp1,cellCoeff1,cellPolyValid1]...
    =generateRandomPolys(vecNumVar,vecDeg,vecNumMon,1);
[cellExp2,cellCoeff2,cellPolyValid2]...
    =generateRandomPolys(vecNumVar,vecDeg,vecNumMon,2);
[cellExp3,cellCoeff3,cellPolyValid3]...
    =generateRandomPolys(vecNumVar,vecDeg,vecNumMon,3);

% From the exponent matrices and the coefficient vectors, get the actual
% monomials.
variables=sdpvar(max(vecNumVar),1);
cellPoly1=cell(sizeNumVar,sizeDeg,sizeNumMon);
cellPoly2=cell(sizeNumVar,sizeDeg,sizeNumMon);
cellPoly3=cell(sizeNumVar,sizeDeg,sizeNumMon);

for i=1:sizeNumVar
    for j=1:sizeDeg
        for k=1:sizeNumMon
            if cellPolyValid1{i,j,k}
                cellPoly1{i,j,k}=polynomialFromExpCoeffVar(...
                    cellExp1{i,j,k},cellCoeff1{i,j,k},...
                    variables(1:vecNumVar(i)));
                cellPoly2{i,j,k}=polynomialFromExpCoeffVar(...
                    cellExp2{i,j,k},cellCoeff2{i,j,k},...
                    variables(1:vecNumVar(i)));
                cellPoly3{i,j,k}=polynomialFromExpCoeffVar(...
                    cellExp3{i,j,k},cellCoeff3{i,j,k},...
                    variables(1:vecNumVar(i)));
            end
        end
    end
end


%% Run the test
numRelaxations=3;

optVal1=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);
runTime1=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);
problem1=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);

optVal2=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);
runTime2=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);
problem2=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);

optVal3=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);
runTime3=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);
problem3=-42*ones(sizeNumVar,sizeDeg,sizeNumMon,numRelaxations);

for i=1:sizeNumVar
    numVar = vecNumVar(i);
    for j=1:sizeDeg
        for k=1:sizeNumMon
            %% Polynomials with standard simplex Newton polytopes
            if cellPolyValid1{i,j,k}
                fprintf(['\n', 'Test a random polynomial with ',...
                    'standard simplex Newton polytope, ', '\n',...
                    'in ',num2str(vecNumVar(i)), ' variables, degree ',...
                    num2str(vecDeg(j)), ' and ', num2str(vecNumMon(k)),...
                    ' monomials.','\n']);
                
                f1=cellPoly1{i,j,k};
                list=variables(1:numVar)';
                %% SOS relaxation
                tic;
                [optVal1(i,j,k,1),problem1(i,j,k,1)]=globalMinSOS(f1,list);
                runTime1(i,j,k,1)=toc;
                
                %% Compute global minimum using the SONC cone.
                tic;
                [optVal1(i,j,k,2),problem1(i,j,k,2)]=globalMinSONC(f1,list);
                runTime1(i,j,k,2)=toc;
                
                %% Compute global minimum using the SOS+SONC cone.
                tic;
                [optVal1(i,j,k,3),problem1(i,j,k,3)]=globalMinSOSpSONC(f1,list);
                runTime1(i,j,k,3)=toc;
            end
            %% Polynomials with arbitrary simplex Newton polytopes.
            if cellPolyValid2{i,j,k}
                fprintf(['\n', 'Test a random polynomial with ',...
                    'simplex Newton polytope, ', '\n',...
                    'in ',num2str(vecNumVar(i)), ' variables, degree ',...
                    num2str(vecDeg(j)), ' and ', num2str(vecNumMon(k)),...
                    ' monomials.','\n']);
                
                f2=cellPoly2{i,j,k};
                list=variables(1:numVar)';
                %% SOS relaxation
                tic;
                [optVal2(i,j,k,1),problem2(i,j,k,1)]=globalMinSOS(f2,list);
                runTime2(i,j,k,1)=toc;
                
                %% Compute global minimum using the SONC cone.
                tic;
                [optVal2(i,j,k,2),problem2(i,j,k,2)]=globalMinSONC(f2,list);
                runTime2(i,j,k,2)=toc;
                
                %% Compute global minimum using the SOS+SONC cone.
                tic;
                [optVal2(i,j,k,3),problem2(i,j,k,3)]=globalMinSOSpSONC(f2,list);
                runTime2(i,j,k,3)=toc;
            end
            %% Polynomials with arbitrary Newton polytopes.
            if cellPolyValid3{i,j,k}
                fprintf(['\n', 'Test a random polynomial with ',...
                    'random Newton polytope, ', '\n',...
                    'in ',num2str(vecNumVar(i)), ' variables, degree ',...
                    num2str(vecDeg(j)), ' and ', num2str(vecNumMon(k)),...
                    ' monomials.','\n']);
                
                f3=cellPoly3{i,j,k};
                list=variables(1:numVar)';
                %% SOS relaxation
                tic;
                [optVal3(i,j,k,1),problem3(i,j,k,1)]=globalMinSOS(f3,list);
                runTime3(i,j,k,1)=toc;
                
                %% Compute global minimum using the SONC cone.
                tic;
                [optVal3(i,j,k,2),problem3(i,j,k,2)]=globalMinSONC(f3,list);
                runTime3(i,j,k,2)=toc;
                
                %% Compute global minimum using the SOS+SONC cone.
                tic;
                [optVal3(i,j,k,3),problem3(i,j,k,3)]=globalMinSOSpSONC(f3,list);
                runTime3(i,j,k,3)=toc;
            end
        end
    end
end

% Get nice arrays containing the test results.
N=sizeNumVar*sizeDeg*sizeNumMon;
dataOptVal1=reshape(permute(optVal1,[3,2,1,4]),N,numRelaxations);
dataRunTime1=reshape(permute(runTime1,[3,2,1,4]),N,numRelaxations);
dataProblem1=reshape(permute(problem1,[3,2,1,4]),N,numRelaxations);

dataOptVal2=reshape(permute(optVal2,[3,2,1,4]),N,numRelaxations);
dataRunTime2=reshape(permute(runTime2,[3,2,1,4]),N,numRelaxations);
dataProblem2=reshape(permute(problem2,[3,2,1,4]),N,numRelaxations);

dataOptVal3=reshape(permute(optVal3,[3,2,1,4]),N,numRelaxations);
dataRunTime3=reshape(permute(runTime3,[3,2,1,4]),N,numRelaxations);
dataProblem3=reshape(permute(problem3,[3,2,1,4]),N,numRelaxations);

% Symbolic arrays.
dataSym1=[sym(dataProblem1)';...
    sym(compose('%8.4g',dataOptVal1))';...
    sym(compose('%8.4g',dataRunTime1))']';
dataSym2=[sym(dataProblem2)';...
    sym(compose('%8.4g',dataOptVal2))';...
    sym(compose('%8.4g',dataRunTime2))']';
dataSym3=[sym(dataProblem3)';...
    sym(compose('%8.4g',dataOptVal3))';...
    sym(compose('%8.4g',dataRunTime3))']';

% LaTex array
latexData1=latex(vpa(dataSym1,4));
latexData2=latex(vpa(dataSym2,4));
latexData3=latex(vpa(dataSym3,4));

% Save the data
save(strcat("data/randomPolynomialsDifferentPolytopes ",...
    string(datetime(datetime,'InputFormat',...
    'yyyy-MM-dd HH:mm:ss.SSS'))));
% % Get LaTex code containing the results.
% latexDataCert=latex(sym(dataCert));
% latexDataRunTime=latex(sym(compose('%8.2f',dataRunTime)));
% latexDataOptVal=latex(sym(compose('%8.6f',dataOptVal)));
% latexDataProblem=latex(sym(dataProblem));

% Save the data
% save(strcat("dataTestsUnconstrained/testArbitraryNewton: ",...
%     string(datetime(datetime,'InputFormat',...
%     'yyyy-MM-dd HH:mm:ss.SSS'))),...
%     "vecNumVar","vecDeg","vecNumMon",...
%     "cellExp","cellCoeff","cellPolyValid",...
%     "cert","optVal","runTime","problem",...
% "dataCert","dataOptVal","dataRunTime","dataProblem",...
% "latexDataCert","latexDataOptVal","latexDataRunTime","latexDataProblem");

