function [PTP, PF, PR, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np, Q);
[m1, ~] = size(Cp);
[n1, n_in] = size(Bp);
A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ap;
A_e(n1+1:n1+m1,1:n1) = Cp*Ap;
B_e = zeros(n1+m1,n_in);
B_e(1:n1,:) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;
C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);

n = n1+m1;
h = zeros(1, n);
F = zeros(1, n);
h(1:n1,:) = C_e;
F(1:n1,:) = C_e*A_e;
% for kk = 2:Np
%     h(kk,:) = h(kk-1,:)*A_e;
%     F(kk,:) = F(kk-1,:)*A_e;
% end
for kk = n1+1:n1:Np*n1
    h(kk:kk+n1-1,:) = h(kk-n1:kk-1,:)*A_e;
    F(kk:kk+n1-1,:) = F(kk-n1:kk-1,:)*A_e;
end

v = h*B_e;
%size(v)

Phi = zeros(Np*n1,Nc);
Phi(:,1) = v;
for i = 2:Nc
    Phi(:,i) = [zeros((i-1)*n1,1);v(1:n1*(Np-i+1),1)];
end
%size(Phi)

bQ = zeros(Np*n1);
for j = 1:n1:Np*n1
    bQ(j:j+n1-1,j:j+n1-1) = Q;
end
%size(bQ)

Ref = ones(Np*n1,1);
%Ref(end-6:end) = 2;
PTP = Phi'*Phi;
PF = Phi'*F;
PR = Phi'*bQ*Ref;

end