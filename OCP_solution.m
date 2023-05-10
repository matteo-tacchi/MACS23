%%%%%%%%% GdR MACS 2023 Spring School %%%%%%%%%
%%%%%%%%%%% Solution to Exercise 1 %%%%%%%%%%%%
% Code by Matteo Tacchi & Corbinian Schlosser %

%% Problem formulation

mset clear % Delete all existing GloptiPoly variables from the workspace

% Choose a maximum degree for the moments
k = 1; 
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

% Declare dynamics
f = u; 

% Declare running cost
L = x^4 + (u^2 - 1)^2;

% Declare feasible set
X = [x^2 <= 1; u^2 <=1; (1-t)*t >= 0; x0 == 0; xT == 0; t0 == 0; tT == 1]; 

% Declare Liouville's PDE
v0 = mmon([t0;x0],d);
vT = mmon([tT;xT],d);
v = mmon([t;x],d);
Lf = [mom( diff(v,t) + diff(v,x)*f ) - mom(vT) + mom(v0) == 0];

% Declare that you are working with probability distributions
momcon = [Lf; mass(mu0) == 1];

% Declare (OCP) instance
P = msdp(min(L), momcon, X);

%% Solving the (OCP)

[status,cost] = msol(P);

double(mvec(mu0)) % Display the moments of mu0
double(mvec(muT)) % Display the moments of muT
double(mvec(nu)) % Display the moments of nu
