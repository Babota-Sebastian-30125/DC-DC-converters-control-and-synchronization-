function [c,ceq] = constraints_LMI(x,params)

ak = params.ak;
%num_vertices = params.num_vertices;
A0 = params.A0;
B1 = params.B1;
N1 = params.N1;
N2 = params.N2;
%N3 = params.N3;
X =  params.X;
alpha = x(end);
gama = params.gama;


ceq = [];

P = [x(1) x(2); x(2) x(3)];
W = [x(4) x(5)]
% P > 0
c = [-eig(P);-alpha];

%6d
for i = 1:length(ak) 
    a = ak(:,:,i);
    M = [1 gama*a'*P; P*a*gama P];
    c = [c; -eig(M)];
end
%6e
for i = 1:length(X)
    xi = X(i,:);
    M = [1 xi; xi' P];
    c = [c; -eig(M)];
end
%6f
for i = 1:length(X)
    xi=X(i,:);
    XN = [ xi*N1; xi*N2];
    NX = [N1'*xi' N2'*xi'];
    M  = 2*alpha*P + gama*(A0*P + P*A0') + gama*(B1*W + W'*B1') + (XN*W + W'*NX);
    c = [c; eig(M)];
end