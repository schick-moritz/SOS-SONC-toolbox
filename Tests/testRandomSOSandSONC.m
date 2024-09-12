close all; clc; clear all;
yalmip('clear');

disp(today('datetime'));
fprintf('Moritz Schick, University of Konstanz\n');
fprintf(['Unconstrained polynomial optimization and \n',...
    'testing membership in the SOS, SONC and SOS+SONC cones \n']);

%% Summary
% In this script, we first construct a collection of SOS and SONC
% polynomials for different instances of number of variables and degrees.
% We also define SOS+SONC polynomials by taking their sum. We then apply
% the SOS, SONC and SOS+SONC relaxation to all of these polynomials.

% Numbers of variables tested.
vecNumVar=[2 3 4 6];
% vecNumVar=[4 5 6];
% vecNumVar=[2];
% vecNumVar=[2];

% Degrees tested.
vecDeg=[6 8 10];
% vecDeg=[6 8];
% vecDeg=[6];
% vecDeg=[10 12];

% Number of squares of SOS polynomial tested.
vecNumSquares=1;

% Number of circuits of SONC polynomials tested.
vecNumCirc=1;
fullDim=1; % Circuits should have full-dimensional Newton polytopes.

sizeNumVar=length(vecNumVar);
sizeDeg=length(vecDeg);
sizeNumSquares=length(vecNumSquares);
sizeNumCirc=length(vecNumCirc);

%% Get exponent matrices and coefficients of the random polynomials.
[cellExpSOS,cellCoeffSOS] ...
    = generateRandomSOS(vecNumVar,vecDeg,vecNumSquares);
[cellExpSONC,cellCoeffSONC,cellPolyValidSONC] ...
    = generateRandomSONC(vecNumVar,vecDeg,vecNumCirc,fullDim);

%% From the exponent matrices and coefficient vectors, get the monomials.
variables=sdpvar(max(vecNumVar),1);
cellPolySOS=cell(sizeNumVar,sizeDeg,sizeNumSquares);
cellPolySONC=cell(sizeNumVar,sizeDeg,sizeNumCirc);
cellPolySOSpSONC=cell(sizeNumVar,sizeDeg,sizeNumSquares,sizeNumCirc);

for i=1:sizeNumVar
    for j=1:sizeDeg
        for k=1:sizeNumSquares
            cellPolySOS{i,j,k}=polynomialFromExpCoeffVar(...
                cellExpSOS{i,j,k},cellCoeffSOS{i,j,k},...
                variables(1:vecNumVar(i)));
        end
        for l=1:sizeNumCirc
            if cellPolyValidSONC{i,j,l}
                cellPolySONC{i,j,l}=polynomialFromExpCoeffVar(...
                    cellExpSONC{i,j,l},cellCoeffSONC{i,j,l},...
                    variables(1:vecNumVar(i)));
                for k=1:sizeNumSquares
                    cellPolySOSpSONC{i,j,k,l}=cellPolySOS{i,j,k}+...
                        cellPolySONC{i,j,l};
                end
            end
        end
    end
end


%% Run the test
numRelaxations = 3;

optValSOS=-42*ones(sizeNumVar,sizeDeg,sizeNumSquares,numRelaxations);
runTimeSOS=-42*ones(sizeNumVar,sizeDeg,sizeNumSquares,numRelaxations);
problemSOS=-42*ones(sizeNumVar,sizeDeg,sizeNumSquares,numRelaxations);

optValSONC=-42*ones(sizeNumVar,sizeDeg,sizeNumCirc,numRelaxations);
runTimeSONC=-42*ones(sizeNumVar,sizeDeg,sizeNumCirc,numRelaxations);
problemSONC=-42*ones(sizeNumVar,sizeDeg,sizeNumCirc,numRelaxations);

optValSOSpSONC=-42*...
    ones(sizeNumVar,sizeDeg,sizeNumSquares,sizeNumCirc,numRelaxations);
runTimeSOSpSONC=-42*...
    ones(sizeNumVar,sizeDeg,sizeNumSquares,sizeNumCirc,numRelaxations);
problemSOSpSONC=-42*...
    ones(sizeNumVar,sizeDeg,sizeNumSquares,sizeNumCirc,numRelaxations);

for i=1:sizeNumVar
    numVar = vecNumVar(i);
    for j=1:sizeDeg
        %% Test SOS polynomials
        fprintf(['\n', 'Test a random SOS polynomial in ', ...
            num2str(vecNumVar(i)), ' variables and degree ', ...
            num2str(vecDeg(j)) ,'.\n']);
        for k=1:sizeNumSquares
            fSOS = cellPolySOS{i,j,k};
            list = variables(1:numVar)';
            %% SOS relaxation
            tic;
            [optValSOS(i,j,k,1),problemSOS(i,j,k,1)]...
                =globalMinSOS(fSOS,list);
            runTimeSOS(i,j,k,1)=toc;
            
            %% SONC relaxation
            tic;
            [optValSOS(i,j,k,2),problemSOS(i,j,k,2)]...
                =globalMinSONC(fSOS,list);
            runTimeSOS(i,j,k,2)=toc;
            
            %% SOS+SONC relaxation
            tic;
            [optValSOS(i,j,k,3),...
                problemSOS(i,j,k,3)]...
                =globalMinSOSpSONC(fSOS,list);
            runTimeSOS(i,j,k,3)=toc;
        end
        
        %% Test SONC polynomials
        fprintf(['\n', 'Test a random SONC polynomial in ', ...
            num2str(vecNumVar(i)), ' variables and degree ', ...
            num2str(vecDeg(j)) ,'.\n']);
        for l=1:sizeNumCirc
            if cellPolyValidSONC{i,j,l}
                fSONC = cellPolySONC{i,j,l};
                list = variables(1:numVar)';
                %% SOS relaxation
                tic;
                [optValSONC(i,j,l,1), problemSONC(i,j,l,1)]...
                    =globalMinSOS(fSONC,list);
                runTimeSONC(i,j,l,1)=toc;
                
                %% SONC relaxation
                tic;
                [optValSONC(i,j,l,2), problemSONC(i,j,l,2)]...
                    =globalMinSONC(fSONC,list);
                runTimeSONC(i,j,l,2)=toc;
                
                %% SOS+SONC relaxation
                tic;
                [optValSONC(i,j,l,3), problemSONC(i,j,l,3)]...
                    =globalMinSOSpSONC(fSONC,list);
                runTimeSONC(i,j,l,3)=toc;
                
                %% Test SOS+SONC polynomials
                fprintf(['\n', 'Test a random SOS+SONC polynomial in ', ...
                    num2str(vecNumVar(i)), ' variables and degree ', ...
                    num2str(vecDeg(j)) ,'.\n']);
                for k=1:sizeNumSquares
                    fSOS = cellPolySOS{i,j,k};
                    f=fSOS+fSONC;
                    list = variables(1:numVar)';
                    %% SOS relaxation
                    tic;
                    [optValSOSpSONC(i,j,k,l,1),...
                        problemSOSpSONC(i,j,k,l,1)]...
                        =globalMinSOS(f,list);
                    runTimeSOSpSONC(i,j,k,1)=toc;
                    
                    %% SONC relaxation
                    tic;
                    [optValSOSpSONC(i,j,k,l,2),...
                        problemSOSpSONC(i,j,k,l,2)]...
                        =globalMinSONC(f,list);
                    runTimeSOSpSONC(i,j,k,2)=toc;
                    
                    %% SOS+SONC relaxation
                    tic;
                    [optValSOSpSONC(i,j,k,l,3),...
                        problemSOSpSONC(i,j,k,l,3)]...
                        =globalMinSOSpSONC(f,list);
                    runTimeSOSpSONC(i,j,k,3)=toc;
                end
            end
        end
    end
end

% Get nice arrays containing the test results.
NSOS=sizeNumVar*sizeDeg*sizeNumSquares;
dataOptValSOS=reshape(permute(optValSOS,[3,2,1,4]),NSOS,numRelaxations);
dataRunTimeSOS...
    =reshape(permute(runTimeSOS,[3,2,1,4]),NSOS,numRelaxations);
dataProblemSOS...
    =reshape(permute(problemSOS,[3,2,1,4]),NSOS,numRelaxations);

NSONC=sizeNumVar*sizeDeg*sizeNumCirc;
dataOptValSONC...
    =reshape(permute(optValSONC,[3,2,1,4]),NSONC,numRelaxations);
dataRunTimeSONC...
    =reshape(permute(runTimeSONC,[3,2,1,4]),NSONC,numRelaxations);
dataProblemSONC...
    =reshape(permute(problemSONC,[3,2,1,4]),NSOS,numRelaxations);

NSOSpSONC=sizeNumVar*sizeDeg*sizeNumSquares*sizeNumCirc;
dataOptValSOSpSONC=reshape(permute(optValSOSpSONC,[4,3,2,1,5]),...
    NSOSpSONC,numRelaxations);
dataRunTimeSOSpSONC=reshape(permute(runTimeSOSpSONC,[4,3,2,1,5]),...
    NSOSpSONC,numRelaxations);
dataProblemSOSpSONC...
    =reshape(permute(problemSOSpSONC,[4,3,2,1,5]),NSOS,numRelaxations);


% Symbolic arrays.
dataSOSSym=[sym(dataProblemSOS)';...
    sym(compose('%8.4g',dataOptValSOS))';...
    sym(compose('%8.4g',dataRunTimeSOS))']';
dataSONCSym=[sym(dataProblemSONC)';...
    sym(compose('%8.3g',dataOptValSONC))';...
    sym(compose('%8.3g',dataRunTimeSONC))']';
dataSOSpSONCSym=[sym(dataProblemSOSpSONC)';...
    sym(compose('%8.3g',dataOptValSOSpSONC))';...
    sym(compose('%8.3g',dataRunTimeSOSpSONC))']';

% LaTexArrays
latexDataSOS=latex(vpa(dataSOSSym,4));
latexDataSONC=latex(vpa(dataSONCSym,4));
latexDataSOSpSONC=latex(vpa(dataSOSpSONCSym,4));

% Save the data
save(strcat("data/randomSOS_SONC_and_SOSpSONC_polynomials ",...
    string(datetime(datetime,'InputFormat',...
    'yyyy-MM-dd HH:mm:ss.SSS'))));

% Save the data.
% save(strcat("dataTestsUnconstrained/testSOSandSONCandSOSpSONC: ",...
%     string(datetime(datetime,'InputFormat',...
%     'yyyy-MM-dd HH:mm:ss.SSS'))),...
%     "vecNumVar","vecDeg","vecNumSquares","vecNumCirc",...
%     "cellExpSOS","cellCoeffSOS",...
%     "cellExpSONC","cellCoeffSONC","cellPolyValidSONC",...
%     "certSOS","optValSOS","runTimeSOS","problemSOS",...
%     "certSONC","optValSONC","runTimeSONC","problemSONC",...
%     "certSOSpSONC","optValSOSpSONC","runTimeSOSpSONC","problemSOSpSONC",...
%     "dataCertSOS","dataOptValSOS","dataRunTimeSOS","dataProblemSOS",...
%     "dataCertSONC","dataOptValSONC","dataRunTimeSONC","dataProblemSONC",...
%     "dataCertSOSpSONC","dataOptValSOSpSONC","dataRunTimeSOSpSONC",...
%     "dataProblemSOSpSONC");
