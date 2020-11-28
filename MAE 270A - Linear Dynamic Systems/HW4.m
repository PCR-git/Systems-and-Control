
% Peter Racioppo
% Linear Systems - HW4

% -----------------------------------------------
% QS 1A

% Perform Gram-Schmidt on columns of A

A = [1,1,1;...
     1,1,1;....
     1,2,3;...
     1,2,4];
 
v1 = A(:,1);
v2 = A(:,2);
v3 = A(:,3);

u1 = v1/norm(v1);

u2_tild = v2 - dot(v2,u1)*u1;
u2 = u2_tild/norm(u2_tild);

u3_tild = v3 - dot(v3,u1)*u1 - dot(v3,u2)*u2;
u3 = u3_tild/norm(u3_tild);

U = [u1,u2,u3];
R = transpose(U)*A;
Q = U;

R_inv = inv(R);
y = [2;1;1;2];

x_LS = R_inv*transpose(Q)*y;

% -----------------------------------------------
% -----------------------------------------------
% QS 5

syms t
A1 = [0,1;-10,-0.8];
A2 = [-1,3;0,-1];

[T1,D1] = eig(A1); % V = matrix of eigenvectors; D = eigenvalues on diagonals
[T2,D2] = eig(A2); % V = matrix of eigenvectors; D = eigenvalues on diagonals
% D1,D2 have negative real parts => should be stable
% T1 has rank 2. T2 has rank 1.

% -----------------------------------------------
% Computing Matrix Exponentials

% To find the first matrix exponential, using a similarity transformation.
Q1 = (T1\A1)*T1; % Similar diagonal matrix
%Q1 = D1;
e_Q1 = [exp(Q1(1,1)*t),0;0,exp(Q1(2,2)*t)];
e_A1T = T1*e_Q1/T1; % Matrix Exponential 1

% To find the second matrix exponential, decompose A2 into 2 matrices
% Let A2 = -I2 + B
B = [0,3;0,0];
% e_A2T := e^At = (e^-I2*t)(e^Bt)
% eL := e^-I2*t
eL = exp(-t)*eye(2);
% B^2 = 0 => e_BT := e^Bt = I + Bt
e_BT = [1,3*t;0,1];
e_A2T = eL*e_BT; % Matrix Exponential 2

% -----------------------------------------------
% Computing Trajectories

% Initial Conditions
x1 = [0;1];
x2 = [1;0];
x3 = [-1;-1];
x4 = [-1;1/2];

% Linspace of time
t1 = linspace(0,10,500);

% Trajectories, A1
Traj1A = e_A1T*x1; % Trajectory 1
Traj1A_sub = real(double(subs(Traj1A,t,t1)));
Traj1B = e_A1T*x2; % Trajectory 2
Traj1B_sub = real(double(subs(Traj1B,t,t1)));
Traj1C = e_A1T*x3; % Trajectory 3
Traj1C_sub = real(double(subs(Traj1C,t,t1)));
Traj1D = e_A1T*x4; % Trajectory 4
Traj1D_sub = real(double(subs(Traj1D,t,t1)));

% Trajectories, A2
Traj2A = e_A2T*x1; % Trajectory 1
Traj2A_sub = real(double(subs(Traj2A,t,t1)));
Traj2B = e_A2T*x2; % Trajectory 2
Traj2B_sub = real(double(subs(Traj2B,t,t1)));
Traj2C = e_A2T*x3; % Trajectory 3
Traj2C_sub = real(double(subs(Traj2C,t,t1)));
Traj2D = e_A2T*x4; % Trajectory 4
Traj2D_sub = real(double(subs(Traj2D,t,t1)));

% -----------------------------------
% Computing Vectors for Use in Drawing Arrows

% Scale factors for the arrows
scale1A = 1;
scale1B = 0.8;
scale1C = 0.8;
scale1D = 0.8;

% A1:
Arrow1Ax = scale1A*([Traj1A_sub(1,:),0]-[0,Traj1A_sub(1,:)]);
Arrow1Ay = scale1A*([Traj1A_sub(2,:),0]-[0,Traj1A_sub(2,:)]);
Arrow1Ax(end) = [];
Arrow1Ay(end) = [];
Arrow1Ax(1) = 0;
Arrow1Ay(1) = 0;

Arrow1Bx = scale1B*([Traj1B_sub(1,:),0]-[0,Traj1B_sub(1,:)]);
Arrow1By = scale1B*([Traj1B_sub(2,:),0]-[0,Traj1B_sub(2,:)]);
Arrow1Bx(end) = [];
Arrow1By(end) = [];
Arrow1Bx(1) = 0;
Arrow1By(1) = 0;

Arrow1Cx = scale1C*([Traj1C_sub(1,:),0]-[0,Traj1C_sub(1,:)]);
Arrow1Cy = scale1C*([Traj1C_sub(2,:),0]-[0,Traj1C_sub(2,:)]);
Arrow1Cx(end) = [];
Arrow1Cy(end) = [];
Arrow1Cx(1) = 0;
Arrow1Cy(1) = 0;

Arrow1Dx = scale1D*([Traj1D_sub(1,:),0]-[0,Traj1D_sub(1,:)]);
Arrow1Dy = scale1D*([Traj1D_sub(2,:),0]-[0,Traj1D_sub(2,:)]);
Arrow1Dx(end) = [];
Arrow1Dy(end) = [];
Arrow1Dx(1) = 0;
Arrow1Dy(1) = 0;

scale2A = 6;
scale2B = 3;
scale2C = 4;
scale2D = 2;

% A2:
Arrow2Ax = scale2A*([Traj2A_sub(1,:),0]-[0,Traj2A_sub(1,:)]);
Arrow2Ay = scale2A*([Traj2A_sub(2,:),0]-[0,Traj2A_sub(2,:)]);
Arrow2Ax(end) = [];
Arrow2Ay(end) = [];
Arrow2Ax(1) = 0;
Arrow2Ay(1) = 0;

Arrow2Bx = scale2B*([Traj2B_sub(1,:),0]-[0,Traj2B_sub(1,:)]);
Arrow2By = scale2B*([Traj2B_sub(2,:),0]-[0,Traj2B_sub(2,:)]);
Arrow2Bx(end) = [];
Arrow2By(end) = [];
Arrow2Bx(1) = 0;
Arrow2By(1) = 0;

Arrow2Cx = scale2C*([Traj2C_sub(1,:),0]-[0,Traj2C_sub(1,:)]);
Arrow2Cy = scale2C*([Traj2C_sub(2,:),0]-[0,Traj2C_sub(2,:)]);
Arrow2Cx(end) = [];
Arrow2Cy(end) = [];
Arrow2Cx(1) = 0;
Arrow2Cy(1) = 0;

Arrow2Dx = scale2D*([Traj2D_sub(1,:),0]-[0,Traj2D_sub(1,:)]);
Arrow2Dy = scale2D*([Traj2D_sub(2,:),0]-[0,Traj2D_sub(2,:)]);
Arrow2Dx(end) = [];
Arrow2Dy(end) = [];
Arrow2Dx(1) = 0;
Arrow2Dy(1) = 0;

% Parameterizing Circle
n = 1000;
t_circ_1 = linspace(0,2*pi,n);
x_circ_1 = cos(t_circ_1);
y_circ_1 = sin(t_circ_1);

% ------------------------------------
% Deriving Delta Neighborhood using Singular Value Bounds

r1 = 0.8;

% Matrix 1
%e_Q1_conj = conj(e_Q1);
e_Q1_conj = [exp(-t*((246^(1/2)*1i)/5 + 2/5)), 0;...
          0, exp(t*((246^(1/2)*1i)/5 - 2/5))];
e_Q1_star = transpose(e_Q1_conj);
e_Q1_norm = simplify(e_Q1_star*e_Q1);
e_Q1_norm1 = [exp(-(4*t)/5), 0;...
              0, exp(-(4*t)/5)];
sigma_exp = sqrt(max(max(subs(e_Q1_norm1,t,0))));
% exp(-(4*t)/5) it's maximized at t = 0
% => max eig = 1  =>  max singular value is 1

sigma_T1 = sqrt(max(eig(transpose(conj(T1))*T1)));
sigma_T1inv = sqrt(max(eig(transpose(conj(inv(T1)))*inv(T1))));

sigma_product = double(sigma_T1*sigma_exp*sigma_T1inv);
delta1 = r1/sigma_product;
r1B = delta1;

% exp_norm1 = (transpose(conj(T1*e_Q1/T1)))*(T1*e_Q1/T1);
% simplify(eig(exp_norm1))
% exp_norm_sub = simplify(subs(exp_norm1,t,0))

% Matrix 2
% exp_norm2 = transpose(conj(e_A2T))*e_A2T
% eig(exp_norm2)
eig(exp(-2*t)*t^2*transpose(conj(B))*B);
delta2 = r1/2.0621;

% ------------------------------------
% Plots

% Figure 1
figure;
% Plotting Trajectories
plot(Traj1A_sub(1,:),Traj1A_sub(2,:),'b-');
hold on;
plot(Traj1B_sub(1,:),Traj1B_sub(2,:),'r-');
plot(Traj1C_sub(1,:),Traj1C_sub(2,:),'g-');
% plot(Traj1D_sub(1,:),Traj1D_sub(2,:),'m-');

% Plotting Arrrows
quiver(Traj1A_sub(1,:),Traj1A_sub(2,:),Arrow1Ax,Arrow1Ay,'b-');
quiver(Traj1B_sub(1,:),Traj1B_sub(2,:),Arrow1Bx,Arrow1By,'r-','AutoScale','off');
quiver(Traj1C_sub(1,:),Traj1C_sub(2,:),Arrow1Cx,Arrow1Cy,'g-','AutoScale','off');
% quiver(Traj1D_sub(1,:),Traj1D_sub(2,:),Arrow1Dx,Arrow1Dy,'m-','AutoScale','off');

% Plotting Circle
plot(r1*x_circ_1,r1*y_circ_1,'k--');
plot(r1B*x_circ_1,r1B*y_circ_1,'k.-');

grid on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);
xlabel('x','FontSize',20,'FontAngle', 'italic');
ylabel('y','FontSize',20,'FontAngle', 'italic');
title('Phase Plane Portrait','FontSize',20);
set(gca, 'FontSize',20,'FontName', 'Times New Roman')

% Figure 2
figure;
% Plotting Trajectories
plot(Traj2A_sub(1,:),Traj2A_sub(2,:),'b-');
hold on;
plot(Traj2B_sub(1,:),Traj2B_sub(2,:),'r-');
plot(Traj2C_sub(1,:),Traj2C_sub(2,:),'g-');
plot(Traj2D_sub(1,:),Traj2D_sub(2,:),'m-');

% Plotting Arrrows
quiver(Traj2A_sub(1,:),Traj2A_sub(2,:),Arrow2Ax,Arrow2Ay,'b-');
quiver(Traj2B_sub(1,:),Traj2B_sub(2,:),Arrow2Bx,Arrow2By,'r-','AutoScale','off');
quiver(Traj2C_sub(1,:),Traj2C_sub(2,:),Arrow2Cx,Arrow2Cy,'g-','AutoScale','off');
quiver(Traj2D_sub(1,:),Traj2D_sub(2,:),Arrow2Dx,Arrow2Dy,'m-','AutoScale','off');

% Plotting Circle
plot(r1*x_circ_1,r1*y_circ_1,'k--');
plot(r2B*x_circ_1,r2B*y_circ_1,'k.-');

grid on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);
xlabel('x','FontSize',20,'FontAngle', 'italic');
ylabel('y','FontSize',20,'FontAngle', 'italic');
title('Phase Plane Portrait','FontSize',20);
set(gca, 'FontSize',20,'FontName', 'Times New Roman')

% ------------------------------------
