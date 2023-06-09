%%%%%%%%% GdR MACS 2023 Spring School %%%%%%%%%
%%%%%%%%%%% Solution to Exercise 3 %%%%%%%%%%%%
% Code by A. Oustry, M. Tacchi & C. Schlosser %

% SDP solver
SDPsolver = 'sedumi'; % or mosek

% degree of decision variables
k = 5; 
d=2*k;

% scaling parameter
alpha = 1.02;

% Variables
x = sdpvar(2,1);

% bound on occupation measure mass
r = 10;

% Define polynomials w(x) and v(t,x)
[w,cw] = polynomial(x,d);
[v,cv] = polynomial(x,d);
s = sdpvar(1);

% Dynamics
f = -[ 2*x(2) ; -0.8*x(1) - 10*(alpha^2*x(1)^2-0.20)*x(2) ] ;

% [0,T] x X
gx = 1 - x'*x ;

% SOS multipliers
[p,cp] = polynomial(x,d-2);
[q1,c1] = polynomial(x,d-2);
[q2,c2] = polynomial(x,d-2);
[q3,c3] = polynomial(x,d-2);

% Lebesgue moments on X
l = [2*pi];
for deg = 1 : d
    for y_ind = 0 : deg
        x_ind = deg - y_ind ;
        
        if  mod(x_ind,2)==0 && mod(y_ind,2)==0
            
            l=[l; 2*(1^(2+x_ind+y_ind))*beta((y_ind+1)/2,(x_ind+1)/2)/(x_ind+y_ind+2)  ];
            
        else
            l=[l; 0] ;
        end
        
    end
    
end


Lv = jacobian(v,x)*f;

% Constraints 
con = [ sos(v - p*gx)];
con = [ con ; sos( w - q1*gx); sos(q1)];
con = [ con ; sos(w - 1 - v - q2*gx); sos(q2)];
con = [ con ; sos(s - Lv - q3* gx); sos(q3)];
con = [ con ; s>=0];

% Objective
obj = cw'*l + r * s * l(1); % minimization of int w d_lambda

% Solver parameteres
options = sdpsettings('solver',SDPsolver,'verbose',1);

% Solve
solvesos(con,obj,options,[s;cw;cv;cp;c1;c2;c3]);

% Retrieve coefficients of w and V
cw = double(cw);
cv = double(cv);
s=double(s)

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
                       
% Plot the state constraint set
theta = linspace(0, 2*pi, 1000);
x = cos(theta);
y = sin(theta);
plot(x, y, 'k', 'lineWidth',2); 
