%%%%%%%%% GdR MACS 2023 Spring School %%%%%%%%%
%%%%%%%%%%% Solution to Exercise 1 %%%%%%%%%%%%
% Code by Matteo Tacchi & Corbinian Schlosser %

%% Plotting f

fplot(@(x) 8*x^5 - 9*x^3 + 2*x + 14*x^6 - 20*x^4 + 7*x^2,[-1 1])

%% Problem formulation

mset clear % Delete all existing GloptiPoly variables from the workspace

mpol x; % Declare decision variable

f = 8*x^5 - 9*x^3 + 2*x + 14*x^6 - 20*x^4 + 7*x^2; % Declare objective function

X = [x^2 <= 1]; % Declare feasible set

P = msdp(min(f), X); % Declare (POP) instance

%% Solving the (POP)

mset('yalmip',true);

[status,obj] = msol(P) % Solve (POP) instance

meas % Check the structure of the computed moment distribution

double(mvec(meas)) % Display the computed moments

double(x) % Extract minimizers
