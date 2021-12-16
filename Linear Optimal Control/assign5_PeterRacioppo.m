
% MAE 270B
% Peter Racioppo

% ------------------------------
% HW 5

% Problem definition of the state space model
% x1' = x3
% x2' = x4
% x3' = u
% x4' = (6*9.8)*x2 + 6*u
% X' = [0,0,1,0; 0,0,0,1; 0,0,0,0; 0,6*9.8,0,0]*X + [0;0;1;6]*u;

A = [0,0,1,0; 0,0,0,1; 0,0,0,0; 0,6*9.8,0,0];
B = [0;0;1;6];

q1 = 1;
q2 = 100;

% q1 = 4;
% q2 = 1000;

R = 1;
Q = zeros(4,4);
Q(1,1) = q1; Q(2,2) = q2;
[K,~,~] = lqr(A,B,Q,R);

C = eye(4);
D = zeros(4, 1);
sys1 = ss(A-B*K, B, C, D);

T = [0:0.01:10]';
v = zeros(size(T,1),1);
x0 = [1 0 0 0]';
x = lsim(sys1,v,T,x0);
u = -x*K';

idx = find(T==4);
xa4 = x(idx:end,1);
x1max = max(xa4) % Max value of x1 after 4 seconds
x2max = max(x(:,2)) % Max value of x2 at any time

figure;
subplot(211);
plot(T,x(:,1:2));
% xlim([4,10])
grid, xlabel('time (sec)');
legend('horizontal position of mass center (m)','angle from vertical (rad)');
subplot(212);
plot(T,u);
% xlim([4,10])
grid, xlabel('time (sec)');
legend('control command');
