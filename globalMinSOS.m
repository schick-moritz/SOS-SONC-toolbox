function [cert, val, problem, fSOS, hSOS] = globalMinSOS(f, vecVar)
% GLOBALMINSOS computes a lower bound for the global minimum of a given
% polynomial using the sums of squares cone.
%
%   Given a real polynomial 'f' defined in the variables given in 'vecVar',
%   we compute a lower bound for the global minimum of f. Therefore, we 
%   solve the SOS-relaxation
%       lambda* = sup lambda s.t. f-lambda is SOS
%   of the unconstrained polynomial optimization problem
%       inf f = sup lambda s.t. f-lambda is PSD.
%  
%   This function can also be used to check membership in the SOS cone: If 
%   lambda* is nonnegative, then f is SOS. 
%   The variables used must be YALMIP sdpvar decision variables.
%
%   Input:
%   - f: the given polynomial.
%   - vecVar: vector of variables in which 'f' is defined. Can be row or
%   column vector.
%
%   Output:
%   - cert: 1 if SOS-relaxation is solved successfully, else 0.
%   - val: val=lambda* is the optimal value of the SOS-relaxation.
%   - problem: error code from the YALMIP 'optimize' function.
%   - fSOS: optional parameter. Gives the SOS polynomial f-lambda*, if
%   SOS-relaxation is solved successfully, else -1.
%   - hSOS: optional parameter. Gives a vector containing all (at most)
%   half degree polynomials yielding an SOS decomposition of fSOS, if
%   SOS-relaxation is solved successfully.

sdpvar objVar;
% Polynomial in the constraint of the SOS-relaxation, which is to be SOS.
fConst = f - objVar;
% Flip cost function, since 'optimize' computes a minimum.
obj = -objVar;

% Initialize the output.
cert = 0;
val = -1;
hSOS = [];

% SOS constraints are given by a Gram matrix approach, using the function
% 'halfNewtonPolytope'.

[monHalfNew, ~] = halfNewtonPolytope(fConst, vecVar);
Q = sdpvar(length(monHalfNew));

% Gram matrix decomposition of constraint polynomial.
fSOS = monHalfNew' * Q * monHalfNew;
% Assemble the constraint set.

const = [coefficients(fConst - fSOS, vecVar) == 0, Q >= 0];

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
        % Get coefficients of fSOS numerically.
        [coefffSOS, monomialsfSOS] = coefficients(fSOS, vecVar);
        coefffSOS = value(coefffSOS);
        fSOS = coefffSOS' * monomialsfSOS;
    end
    if nargout>=5
        % Compute explicit SOS decomposition.
        hSOS = sosd(fSOS);
    end
elseif diagnostics.problem == 1
    disp('YALMIP: Infeasible problem.')
else
    disp('YALMIP: Problem occurred:')
    disp(['Problem: ',num2str(diagnostics.problem)]);
end
problem=diagnostics.problem;

if cert==0
    fSOS = -1;
    hSOS = [];
end
end