clear; close all;
%mpc regulation of segway. continuous model:
Ac = [0 -18.74 0;
    0 0 1;
    0 21.99 0];
Bc = [0.09975 0 -0.0844]';
Cc = eye(3);
Dc = zeros(3,1);
dt = 0.05;

%obtain discrete state space
sys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sys,dt);
A = sysd.A;
B = sysd.B;
C = sysd.C;

%control and prediction horizons
Np = 20;
Nc = Np;

Q = [50 0 0;
    0 1 0;
    0 0 1];
tic;
[PTP, PF, PR, Ae, Be, Ce] = mpcgain(A,B,C,Nc,Np,Q);
toc;
%augmented state size
n = length(Be);

%initial augmented feedback
Xf = zeros(n,1);

%initial plant state
xm = [-1 0 0]';

%simulation time vector
simTime = 20;
tvec = 0:dt:simTime;
nSim = length(tvec);

%setpoint, initial u and y
r = zeros(nSim,1);
y = xm; u = 0;


%control penalty, inverted hessian
R = 0.01;
cPen = R*eye(Nc);
ihes = inv(PTP+cPen);
Kmpc = ihes*PF;
u1 = zeros(1, nSim);
y1 = zeros(3, nSim);
x = zeros(3, nSim);

disp(eig(Ae-Be*Kmpc(1,:)));

for k = 1:nSim
    
    %get change in u
    deltaU = -Kmpc*Xf;
    du = deltaU(1,1);
    
    %find control input, store 
    u = u+du;
    u1(k) = u;
    y1(:,k) = y;
    
    %store prev plant state, run plant
    x(:,k) = xm;
    xPrev = xm;
    xm = A*xm + B*u;
    y = C*xm;
    
    %augmented system states
    Xf = [xm-xPrev; y];
end

figure;
subplot(311);hold on;
plot(tvec,y1(1,:));
plot(tvec,r,'--');ylabel('[m]');
legend('x position','reference');
subplot(312);hold on;
plot(tvec,y1(2,:));
plot(tvec,r,'--');ylabel('[rad]');
legend('rod angle','reference');
subplot(313);hold on;
plot(tvec,y1(3,:));
plot(tvec,r,'--');ylabel('[rad/s]');
legend('rod rotation speed','reference');
xlabel('Time [s]');

figure();
stairs(tvec,u1,'r');
legend('control input');
ylabel('Torque [Nm]');
xlabel('Time [s]');

