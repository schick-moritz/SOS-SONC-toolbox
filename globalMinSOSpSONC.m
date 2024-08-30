function [cert, val, problem, fSOS, fSONC, hSOS, decompVec, hSONC] = ...
    globalMinSOSpSONC(f, vecVar)
% GLOBALMINSOSpSONCNC computes a lower bound for the global minimum of a 
% given polynomial using the sums of squares plus sums of nonnegative
% circuit polynomial cone.
%
%   Given a real polynomial 'f' defined in the variables given in 'vecVar',
%   we compute a lower bound for the global minimum of f. Therefore, we
%   solve the (SOS+SONC)-relaxation
%       lambda* = sup lambda s.t. f-lambda is SOS+SONC
%   of the unconstrained polynomial optimization problem
%       inf f = sup lambda s.t. f-lambda is PSD.
%
%   This function can also be used to check membership in the SOS+SONC 
%   cone: If lambda* is nonnegative, then f is SOS+SONC.
%   The variables used must be YALMIP sdpvar decision variables.
%
%   Input:
%   - f: the given polynomial.
%   - vecVar: vector of variables in which 'f' is defined. Can be row or
%   column vector.
%
%   Output:
%   - cert: 1 if (SOS+SONC)-relaxation is solved successfully, else 0.
%   - val: val=lambda* is the optimal value of the (SOS+SONC)-relaxation.
%   - problem: error code from the YALMIP 'optimize' function.
%   - fSOS: optional parameter. SOS summand of f-lambda*, if
%   (SOS+SONC)-relaxation is solved successfully, else -1.
%   - fSONC: optional parameter. SONC summand of f-lambda*, if
%   (SOS+SONC)-relaxation is solved successfully, else -1.
%   - decompVec: optional parameter. Vector 
%   c(1);...;c(m); nu(1);...;nu(m)] of relative entropy constraint 
%   variables used to decide membership in the polynomial SAGE/ SONC cone. 
%   'm' is the number of exponents/ monomials of 'f'. This output parameter
%   comes from a call to 'constSAGE'.
%   - hSOS: optional parameter. Gives a vector containing all (at most)
%   half degree polynomials yielding an SOS decomposition of fSOS, if
%   SOS-relaxation is solved successfully.
%   - hSONC: optional parameter. Gives a vector containing all polynomial
%   AGE functions in a polynomial SAGE decomposition of 'f'.

sdpvar objVar;
% Constraint polynomial, which is to be SOS+SONC.
fConst = f - objVar;
% Flip cost function, since 'optimize' computes a minimum.
obj = -objVar;

% Initialize the output.
cert = 0;
val = -1;

%% Generate SOS constraints
% We use a Gram matrix approach. We compute all monomials and exponents in
% the Newton polytope of 'f', as we use it anyways for the SONC 
% constraints. 
[monNew, exponentsNew] = newtonPolytope(fConst, vecVar);
numMonNew = length(monNew);
[monHalfNew, ~] = halfNewtonPolytope([], vecVar, exponentsNew);
Q = sdpvar(length(monHalfNew));
fSOS = monHalfNew' * Q * monHalfNew;

%% Generate SONC constraints.
helpExponentsOdd = mod(exponentsNew,2);
helpExponentsOdd = any(helpExponentsOdd);
allIndices = 1:numMonNew;
indInnerTerms = allIndices(helpExponentsOdd);
coeffSONC = sdpvar(numMonNew, 1);
fSONC = coeffSONC' * monNew;
[constSONC, decompVec] = constSAGE(coeffSONC, exponentsNew, indInnerTerms);

%% Assemble the full constraint set
const = [...
    (coefficients(fConst - fSOS - fSONC, vecVar) == 0): 'decomposition',...
    (constSONC): 'SAGE constraints',...
    (Q>=0): 'SOS constraint'];

%% Solve the optimization problem
% Let YALMIP choose a solver on its own.
solverOptions = sdpsettings('solver', '', 'verbose', 1);
% % Use SCS
% solverOptions = sdpsettings('solver','scs','verbose',1,'SCS.max_iters',...
%     100000,'scs.rho_x',1e-5,'scs.eps',1e-5);
% % Use MOSEK.
% solverOptions = sdpsettings('solver', 'mosek', 'verbose', 1);
diagnostics = optimize(const, obj, solverOptions);

if diagnostics.problem == 0
    disp('YALMIP: Successfully solved.')
    val = value(objVar);
    cert = 1;
    if nargout >= 4
        % Get SOS and SONC summands of 'f' numerically.
        [cSOS, vSOS] = coefficients(fSOS, vecVar);
        fSOS = value(cSOS)' * vSOS;
        
        [cSONC, vSONC] = coefficients(fSONC, vecVar);
        fSONC = value(cSONC)' * vSONC;
    end
    if nargout >=6
        % Compute explicit SOS decomposition.
        hSOS = sosd(fSOS);
        % Compute explicit polynomial SAGE decomposition.
        decompVec = value(decompVec);
        hSONC = sdpvar(numMonNew, 1);
        for i = 1:numMonNew
            hSONC(i) = decompVec((numMonNew * (i - 1))...
                + 1:(numMonNew * i))' * monNew;
        end
    end
elseif diagnostics.problem == 1
    disp('YALMIP: Infeasible problem.')
else
    disp('YALMIP: Problem occurred:')
    disp(['Problem: ',num2str(diagnostics.problem)]);
end
problem=diagnostics.problem;

if cert == 0
    fSOS = -1;
    fSONC = -1;
    hSOS = [];
    decompVec = [];
    hSONC = [];
end
end
