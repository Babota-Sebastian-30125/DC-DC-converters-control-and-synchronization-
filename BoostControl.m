%% BOOST
gama_boost = 0.3; % valoarea random
L1=2.57e-3;
rL1 = 130e-3;
rDs1 = 0.01;
C1 = 4.7e-6;
rC1 = 270e-3;
Cin = 3.57e-6;
Vf1 = 0.2;
L2 = 1.71e-3;
rL2 = 110e-3;
rDs2 = 80e-3;
C2 = 3.57e-6;
rC2 = 350e-3;
rCin = 270e-3;
Vf2 = 0.62;
R = 27;
E = 100;

a011_boost = (-1/L1)*(rL1 +(R*rC1)/(R+rC1)+rDs1);
a012_boost = (-1/L1)*(R/(R+rC1));
a021_boost = (1/C1)*(R/(R+rC1));
a022_boost = -1/(C1*(R+rC1));

a111_boost = (R*rC1)/(L1*(R+rC1));
a112_boost = R/(L1*(rC1+R));
a121_boost = -R/(C1*(rC1+R));
a122_boost = 0;

A0_boost = [a011_boost, a012_boost; a021_boost, a022_boost];
A1_boost = [a111_boost, a112_boost; a121_boost, a122_boost];

N1_boost = [a111_boost; a112_boost];
N2_boost = [a121_boost; 0];

B1_boost = [0; 0];    
B0_boost = [E/L1; 0]; 

num_vertices = 4;
iL_limit_min_boost = -5;
iL_limit_max_boost = 5;
vC_limit_min_boost = -50;     
vC_limit_max_boost = 50;

ak_boost = zeros(2, 1, 4); 
ak_boost(:,:,1) = [ 1/iL_limit_max_boost;  0 ];   
ak_boost(:,:,2) = [ 1/iL_limit_min_boost;  0 ];  
ak_boost(:,:,3) = [ 0;  1/vC_limit_max_boost ];    
ak_boost(:,:,4) = [ 0; 1/vC_limit_min_boost ]; 

X_boost = [
  iL_limit_min_boost, vC_limit_min_boost;   
  iL_limit_min_boost, vC_limit_max_boost;  
  iL_limit_max_boost, vC_limit_min_boost;   
  iL_limit_max_boost, vC_limit_max_boost    
];

% Parameter Structure
params_boost = struct();
params_boost.ak = ak_boost;
params_boost.A0 = A0_boost;
params_boost.A1 = A1_boost; 
params_boost.B1 = B1_boost;
params_boost.B0 = B0_boost; 
params_boost.N1 = N1_boost;
params_boost.N2 = N2_boost;
params_boost.X = X_boost;
params_boost.gama = gama_boost;

opts = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','active-set',...
    'MaxIterations',200,...
    'FiniteDifferenceType','forward',...
    'MaxFunctionEvaluations',1e5,...
    'ConstraintTolerance', 1e-10,...
    'FunctionTolerance',1e-6,...
    'StepTolerance',1e-6,...
    'OptimalityTolerance',1e-10',...
    'PlotFcn',{'optimplotx','optimplotfunccount',...  % ,'optimplotfval'
    'optimplotfvalconstr','optimplotconstrviolation',...
    'optimplotstepsize','optimplotfirstorderopt'},...
    'UseParallel',false, ...
    'SpecifyObjectiveGradient',false,...
    'SpecifyConstraintGradient',false, ...
    'DiffMaxChange',1e-2 ... % Inf
    );

%x0 = -ones(6,1);
x0_boost = [1; 0; 1; 0; 0; 0];

[xopt_boost,fval_boost,EXITFLAG_boost,OUTPUT_boost,LAMBDA_boost,GRAD_boost,HESSIAN_boost] = fmincon(...
   @(x) x(end),x0_boost,[],[],[],[],[],[],@(x)constraints_LMI(x,params_boost),opts);

Popt_boost = [xopt_boost(1) xopt_boost(2); xopt_boost(2) xopt_boost(3)]
Wopt_boost = [xopt_boost(4) xopt_boost(5)]
alpha_opt_boost = xopt_boost(6)
kopt_boost = Wopt_boost*inv(Popt_boost)
%%
kopt_boost = [-0.0068964,-0.0011605];