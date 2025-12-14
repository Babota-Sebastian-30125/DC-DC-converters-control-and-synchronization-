
%% BUCK
gama_buck = 0.5;%valoarea random
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
E_buck =300;
a011_buck = (-1/L1)*(rL1 +(R*rC1)/(R+rC1)-rDs1);
a012_buck = (-1/L1)*(R/(R+rC1));
a021_buck = (1/C1)*(R/(R+rC1));
a022_buck = -1/(C1*(R+rC1));
a111_buck = -2/L1*rDs1;
b111_buck = E_buck/L1;
A0_buck = [a011_buck, a012_buck; a021_buck, a022_buck];
N1_buck = [a111_buck; 0];
N2_buck = [0; 0];
B1_buck = [b111_buck; 0];
num_vertices = 4;

iL_limit_min_buck = -20;
iL_limit_max_buck = 20;
vC_limit_buck = 70;
vC_limit_min_buck = -50;

ak_buck = zeros(2, 1, 4); 
ak_buck(:,:,1) = [ 1/iL_limit_max_buck;  0 ];   
ak_buck(:,:,2) = [ 1/iL_limit_min_buck;  0 ];   
ak_buck(:,:,3) = [ 0;  1/vC_limit_buck ];    
ak_buck(:,:,4) = [ 0; 1/vC_limit_min_buck ];   

X_buck = [
 iL_limit_min_buck, vC_limit_min_buck;   
 iL_limit_min_buck,  vC_limit_buck;  
  iL_limit_max_buck, vC_limit_min_buck;   
  iL_limit_max_buck,  vC_limit_buck    
];
params = struct();
params.ak = ak_buck;
%params.num_vertices = num_vertices;
params.A0 = A0_buck;
params.B1 = B1_buck;
params.N1 = N1_buck;
params.N2 = N2_buck;
params.X = X_buck;
%params.alpha = alpha;
params.gama = gama_buck;

opts = optimoptions(@fmincon,'Display','iter',...u
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
x0_buck = [1; 0; 1; 0; 0; 0];
[xopt,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(...
   @(x) x(end),x0_buck,[],[],[],[],[],[],@(x)constraints_LMI(x,params),opts)

Popt_buck = [xopt(1) xopt(2); xopt(2) xopt(3)]
Wopt_buck = [xopt(4) xopt(5)]
alpha_opt_buck = xopt(6)
kopt_buck = Wopt_buck*inv(Popt_buck)
