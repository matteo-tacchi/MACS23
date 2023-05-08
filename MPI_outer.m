% *************************************************************************
% Problem:  Van der Pol oscillator ROA - dual SOS problem
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


% Define polynomials w(x) and v(t,x)
[w,cw] = polynomial(x,d);
[v,cv] = polynomial(x,d);


% Discount factor
b = 1 ;

% Dynamics
f = -[ 2*x(2) ; -0.8*x(1) - 10*(alpha^2*x(1)^2-0.20)*x(2) ] ;


% [0,T] x X
gx = 1 - x'*x ;


% SOS multipliers
[q,cq] = polynomial(x,d-2);
[p,cp] = polynomial(x,d-2);
[s,cs] = polynomial(x,d-2);


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
con = [ sos( b*v-Lv - q*gx) ; sos(q) ];    
con = [ con ; sos(w - v - 1 - p*gx); sos(p) ]; % w >= v + 1 on {0} x X
con = [ con ; sos( w - s*gx) ; sos(s)];   % w >= 0 on X

% Objective
obj = cw'*l; % minimization of int w d_lambda

% Solver parameteres
options = sdpsettings('solver',SDPsolver,'verbose',1);

% Solve
solvesos(con,obj,options,[cw;cv;cq;cp;cs]);

% Retrieve coefficients of w and V
cw = double(cw);
cv = double(cv);

%% Plots

% Level set contour plot of {x ; w(x) = 0}
figure
X = sdpvar(1,1); Y = sdpvar(1,1);
vv = monolist([X;Y],d);
p = vectorize(sdisplay(cv'*vv + 1)); % v(0,x) + 1
[X,Y] = meshgrid(-1:0.005:1,-1:0.005:1);
Z = eval(p);
contour(X,Y,Z, [1 1], '-b', 'linewidth',2); hold on

% Simulate trajectory with reversed time to get the boundary of the true MPI
f_vpo =  @(t,x)([ 2*x(2) ; -0.8*x(1) - 10*(alpha^2*x(1)^2-0.20)*x(2) ]);
[~, xv] = ode45(f_vpo,0:0.01:100,[0.1;0.1]);

plot(xv(2000:2574,1),xv(2000:2574,2),'-r','LineWidth',2)
xlabel('x_1'); ylabel('x_2');
%legend('Outer', 'True')
title('Van der Pol oscillator ROA')





