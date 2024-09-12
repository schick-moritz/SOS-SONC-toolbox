close all; clear all; clc;
yalmip('clear');

disp(today('datetime'));
fprintf('Moritz Schick, University of Konstanz\n');

% Initialize variables.
sdpvar x y
vecVar = [x y];

% Define polynomials.
M = x^4*y^2+x^2*y^4+1-3*x^2*y^2; % Motzkin.
% This is SOS+SONC but not SOS and not SONC.
f = 1/2*(1+2*x*y+x^2*y)^2+M;

% Compute the Newton polytope and one half times the Newton polytope.
[monNewF, expNewF] = newtonPolytope(f,vecVar);
[monHalfNewF, expHalfNewF] = halfNewtonPolytope(M, vecVar, expNewF);

% Solve the SOS, SONC, and SOS+SONC relaxations
[valSOS, problSOS, fSOS1, hSOS1] = globalMinSOS(f, vecVar);
[valSONC, problSONC, fSONC1, decompVec1, hSONC1] ...
    = globalMinSONC(f, vecVar);
[valSOSpSONC, problSOSpSONC, fSOS2, fSONC2, hSOS2,decompVec2, hSONC2] ...
    = globalMinSOSpSONC(f, vecVar);