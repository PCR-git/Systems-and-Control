% Peter Racioppo
% 103953689
% HW 5

syms t

M1 = [-4,0;0,1];
M2 = [0,3;-3,0];

% To find the matrix exponential, use a similarity transformation.
[T,D] = eig(M2);
Q = (T\M2)*T; % Similar diagonal matrix
e_Q1 = [exp(Q(1,1)*t),0;0,exp(Q(2,2)*t)];
e_M2 = simplify(T*e_Q1/T); % Matrix Exponential

e_M2 = [cos(3*t), sin(3*t);...
       -sin(3*t), cos(3*t)];

transpose(e_M2)*e_M2;

e_M2_p = -3*[sin(3*t), -cos(3*t);...
             cos(3*t),  sin(3*t)];

A_tild = simplify(e_M2*(M1*transpose(e_M2)-transpose(e_M2_p)))
EigA = simplify(eig(A_tild))
