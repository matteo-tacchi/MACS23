%%%%%%%%% GdR MACS 2023 Spring School %%%%%%%%%
%%%%%%% Bolza's optimal control problem %%%%%%%
% Code by Matteo Tacchi & Corbinian Schlosser %

%% Problem formulation

mset clear % Delete all existing GloptiPoly variables from the workspace

% Choose a maximum degree for the moments
k = % Enter an integer; the degree will be twice that value (even)
d = 2*k; 

% Declare random variables
mpol t0; % initial time
mpol x0; % initial condition
mpol tT; % terminal time
mpol xT; % terminal condition
mpol t; % running time
mpol x; % running state
mpol u; % running input

% Assign variables to probability distributions
mu0 = meas([t0;x0]); % initial
muT = meas([tT;xT]); % terminal
nu = meas([t;x;u]); % occupation

% Declare dynamics and running cost

% Declare feasible set: polynomial inequality constraints /!\ Do not forget equality constraints on x0, xT, t0, tT!

% Declare Liouville's PDE
v0 = mmon([t0;x0],d);
vT = mmon([tT;xT],d);
v = mmon([t;x],d); % v is the list of all monomials of degree at most d
momcon = []; % Test Liouville PDE against v, v0 := v(0,·), vT := v(T,·) (use mom and diff functions)

% Declare that you are working with probability distributions
momcon = [momcon; mass(mu0) == 1];

% Declare (OCP) instance
P = % Use msdp as for POP, with optimization objective, moment and support constraints (see user's guide)

%% Solving the (OCP)
                                                                                        
mset('yalmip',true);

[status,cost] = % Use msol as for POP
