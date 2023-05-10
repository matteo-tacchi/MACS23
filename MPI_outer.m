% *************************************************************************
% Problem:  Van der Pol oscillator outer MPI approximation - dual SOS problem
% Requires: Gloptipoly 3, Yalmip + SDP solver (tested with MOSEK and SeDuMi) 
% Authors:  Milan Korda, Didier Henrion, 2013
% *************************************************************************


% SDP solver
SDPsolver = 'sedumi'; % or mosek

% degree of decision variables
k = 5; 
d=2*k;

% scaling parameter
alpha = 1.02;

% Variables
x = sdpvar(2,1);


% Define polynomials w(x) and v(x)
[w,cw] = polynomial(x,d);
[v,cv] = polynomial(x,d);


% Discount factor
b = 1 ;

% Dynamics
f = -[ 2*x(2) ; -0.8*x(1) - 10*(alpha^2*x(1)^2-0.20)*x(2) ] ;


% [0,T] x X
gx = 1 - x'*x ;


% SOS multipliers
[q1,c1] = polynomial(x,d-2);
[q2,c2] = polynomial(x,d-2);
[q3,c3] = polynomial(x,d-2);


% Lebesgue moments on X
l = 2*pi;
for deg = 1 : d
    for y_ind = 0 : deg
        x_ind = deg - y_ind ;
        
        if  mod(x_ind,2)==0 && mod(y_ind,2)==0
            
            l=[l; 2*beta((y_ind+1)/2,(x_ind+1)/2)/(x_ind+y_ind+2)  ];
            
        else
            l=[l; 0] ;
        end
        
    end
    
end


Lv = jacobian(v,x)*f;

% Constraints 
con = [ sos( b*v-Lv - q1*gx) ; sos(q1) ];    
con = [ con ; sos(w - v - 1 - q2*gx); sos(q2) ]; % w >= v + 1 on {0} x X
con = [ con ; sos( w - q3*gx) ; sos(q3)];   % w >= 0 on X

% Objective
obj = cw'*l; % minimization of int w d_lambda

% Solver parameteres
options = sdpsettings('solver',SDPsolver,'verbose',1);

% Solve
solvesos(con,obj,options,[cw;cv;c1;c2;c3]);

% Retrieve coefficients of w and v
cw = double(cw);
cv = double(cv);

%% Plots

% Level set contour plot of {x ; w(x) = 0}
figure
X = sdpvar(1,1); Y = sdpvar(1,1);
vv = monolist([X;Y],d);
p = vectorize(sdisplay(cv'*vv + 1)); % v(x) + 1
[X,Y] = meshgrid(-1:0.005:1,-1:0.005:1);
Z = eval(p);
contour(X,Y,Z, [1 1], '-b', 'linewidth',2); hold on

% Simulate trajectory with reversed time to get the boundary of the true MPI
f_vpo =  @(t,x)([ 2*x(2) ; -0.8*x(1) - 10*(alpha^2*x(1)^2-0.20)*x(2) ]);
[~, xv] = ode45(f_vpo,0:0.01:100,[0.1;0.1]);

plot(xv(2000:2574,1),xv(2000:2574,2),'-r','LineWidth',2)
xlabel('x_1'); ylabel('x_2');
%legend('Outer', 'True')
title('Van der Pol oscillator MPI set')





