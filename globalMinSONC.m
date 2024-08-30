function [cert, val, problem, fSONC, decompVec, hSONC] = ...
    globalMinSONC(f, vecVar)
% GLOBALMINSONC computes a lower bound for the global minimum of a given
% polynomial using the sums of nonnegative circuit polynomials cone.
%
%   Given a real polynomial 'f' defined in the variables given in 'vecVar',
%   we compute a lower bound for the global minimum of f. Therefore, we
%   solve the SONC-relaxation
%       lambda* = sup lambda s.t. f-lambda is SONC
%   of the unconstrained polynomial optimization problem
%       inf f = sup lambda s.t. f-lambda is PSD.
%
%   This function can also be used to check membership in the SONC cone: If
%   lambda* is nonnegative, then f is SONC.
%   The variables used must be YALMIP sdpvar decision variables.
%
%   Input:
%   - f: the given polynomial.
%   - vecVar: vector of variables in which 'f' is defined. Can be row or
%   column vector.
%
%   Output:
%   - cert: 1 if SONC-relaxation is solved successfully, else 0.
%   - val: val=lambda* is the optimal value of the SONC-relaxation.
%   - problem: error code from the YALMIP 'optimize' function.
%   - fSONC: optional parameter. Gives the SONC polynomial f-lambda*, if
%   SOS-relaxation is solved successfully, else -1.
%   - decompVec: optional parameter. Vector 
%   c(1);...;c(m); nu(1);...;nu(m)] of relative entropy constraint 
%   variables used to decide membership in the polynomial SAGE/ SONC cone. 
%   'm' is the number of exponents/ monomials of 'f'. This output parameter
%   comes from a call to 'constSAGE'.
%   - hSONC: optional parameter. Gives a vector containing all polynomial
%   AGE functions in a polynomial SAGE decomposition of 'f'.

sdpvar objVar;
% Polynomial in the constraint of the SOS-relaxation, which is to be SOS.
fConst = f - objVar;
% Flip cost function, since 'optimize' computes a minimum.
obj = -objVar;

% Initialize the output.
cert = 0;
val = -1;

% Compute signomial representative.
[~, mon, exponents, indInnerTerms] = sigRep(fConst, vecVar);

% Generate SONC constraints.
numMon = length(mon);
coeff = sdpvar(numMon, 1);
fSONC = coeff' * mon;
[const, decompVec] = constSAGE(coeff, exponents, indInnerTerms);

% Solve optimization problem
const = [coefficients(fConst - fSONC, vecVar) == 0, const];

%% Solve the optimization problem
% Let YALMIP choose a solver on its own.
solverOptions = sdpsettings('solver', '', 'verbose', 1);
% % Use ECOS.
% solverOptions = sdpsettings('solver','ecos','verbose',1,...
%     'ECOS.maxit',1000,'SCS.max_iters',1000);
% % Use MOSEK.
% solverOptions = sdpsettings('solver', 'mosek', 'verbose', 1);
diagnostics = optimize(const, obj, solverOptions);

%% Get the output parameters, if problem is successfully solved.
if diagnostics.problem == 0
    disp('YALMIP: Successfully solved.')
    val = value(objVar);
    cert = 1;
    if nargout>=4
        % Get coefficients of fSONC numerically.
        [coeff, monSONC] = coefficients(fSONC, vecVar);
        coeff = value(coeff);
        fSONC = coeff' * monSONC;
    end
    if nargout>=5
        % Compute explicit polynomial SAGE decomposition.
        decompVec = value(decompVec);
        hSONC = sdpvar(numMon, 1);
        for i = 1:numMon
            hSONC(i) = decompVec((numMon * (i - 1)) + 1 : (numMon * i))' * mon;
        end
    end
elseif diagnostics.problem == 1
    disp('YALMIP: Infeasible problem.')
else
    disp('YALMIP: Problem occurred:')
    disp(['Problem: ',num2str(diagnostics.problem)]);
end
problem=diagnostics.problem;

if cert==0
    fSONC = -1;
    decompVec = [];
    hSONC = [];
end
end
