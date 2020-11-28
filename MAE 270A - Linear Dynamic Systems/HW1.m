
% Linear Systems - Homework 1
% Peter Racioppo
% 103953689


close all
clear all
clc

% ---------------------------------
% Verifying that A is a normal matrix

A = [1,1,0;...
     0,1,1;...
     1,0,1];
A_star = transpose(A); % A is real, so it equals it's complex conjugate
M = A*A_star - A_star*A; % M is a 3x3 zero matrix iff A is normal

% ---------------------------------
% Showing that v1 is an eigenvector of A

v1 = transpose([1 1 1]);
A*v1; % This yields transpose[[2 2 2]], implying that transpose([1 1 1])
% is an eigenvector with eigenvalue 2

% ---------------------------------
% Algorithm to derive the eigenvectors of the 3x3 matrix

syms lambda
I3 = eye(3,3); % 3x3 identity matrix
LHS = A - lambda*I3; % LHS of characteristic equation

% Check Det3 function
%expand(Det3(LHS))
%det(LHS)

char_poly = Det3(LHS); % Compute the characteristic polynomial

lambda_vec = zeros(3,1); % Preallocate vec of eigenvalues
lambda_vec = lambda_vec + solve(char_poly == 0); % Vec of eigenvalues
% Note: we could use Cardan's formula for 3rd order polynomials here,
% instead of "solve"

Mat1 = A - lambda_vec(1)*I3; % Substitute eigenvalue into LHS
Mat2 = rref(Mat1); % Reduce to reduced row echelon form
%Mat3 = Mat2(1,:); % Take a row and find eigenvector

% Thus, v = transpose([1 1 1]) is an eigenvector

% ---------------------------------
% Gram-Schmidt Orthnormalization for the 3x3 matrix A

u1 = normalize(v1); % Normalize the eigenvector

v2 = [rand;rand;rand]; % Generate (2nd) random vector
u2_s = v2 - (transpose(conj(u1))*v2)*u1; % Subtract out part parallel to u1
u2 = normalize(u2_s); % Normalize

v3 = [rand;rand;rand]; % Generate (3rd) random vector
u3_s = v3 - (transpose(conj(u1))*v3)*u1 - (transpose(conj(u2))*v3)*u2; % Subtract out parts parallel to u1 & u2
u3 = normalize(u3_s); % Normalize

T3 = [u1 u2 u3];
rank(T3); % Verify that T3 is full rank
transpose(conj(T3))*T3; % Verify that T3 is unitary

L = transpose(conj(T3))*A*T3;
A2 = L([2 3], [2 3]);
% A([a b], [c d]): Extract the ath and bth rows, and the cth and dth columns

% ---------------------------------
% Algorithm to derive the eigenvectors of the 2x2 matrix

syms lambda2
I2 = eye(2,2); % 2x2 identity matrix
LHS2 = A2 - lambda2*I2; % LHS of characteristic equation

% Check Det3 function
%expand(Det3(LHS))
%det(LHS)

char_poly2 = Det2(LHS2); % Compute the characteristic polynomial

lambda_vec2 = zeros(2,1); % Preallocate vec of eigenvalues
lambda_vec2 = lambda_vec2 + solve(char_poly2 == 0); % Vec of eigenvalues

Mat1 = A2 - lambda_vec2(1)*I2; % Substitute eigenvalue into LHS
Mat2 = rref(Mat1); % Reduce to reduced row echelon form
%Mat3 = Mat2(1,:) % Take a row and find eigenvector

% By inspection, [i; -1] is an eigenvector

% ---------------------------------
% Gram-Schmidt Orthnormalization for the 2x2 matrix A2

v1 = transpose([1i -1]);
u1 = normalize(v1); % Normalize the eigenvector

v2 = [rand;rand]; % Generate (2nd) random vector
u2_s = v2 - (transpose(conj(u1))*v2)*u1; % Subtract out part parallel to u1
u2 = normalize(u2_s); % Normalize

T2 = [u1 u2]; % T2, the matrix that diagonalizes A2
rank(T2); % Verify that T2 is full rank
transpose(conj(T2))*T2; % Verify that T2 is unitary
L2 = transpose(conj(T2))*A2*T2; % Verify that T2 diagonalizes A2

% ---------------------------------
% Constructing U, the unitary matrix that diagonalizes A

RM = zeros(3,3);
RM(1,1) = 1;
RM([2 3], [2 3]) = T2;
U = T3*RM

transpose(conj(U))*U; % Verifying that U is unitary
A_diag = transpose(conj(U))*A*U % Verifying that U diagonalizes A

% ---------------------------------
% Solution

% U = [0.5774 + 0.0000i  -0.1852 + 0.5468i   0.5749 - 0.0531i;...
%      0.5774 + 0.0000i   0.5662 - 0.1130i  -0.2414 + 0.5244i;...
%      0.5774 + 0.0000i  -0.3810 - 0.4338i  -0.3335 - 0.4713i];
% 
% A_diag = [2.0000        0                      0         ;...
%           0      0.5000 - 0.8660i              0         ;...
%           0             0               0.5000 + 0.8660i];
% A_diag := transpose(conj(U))*A*U

% U is not unique, since it depends upon the orthnormal basis used to
% construct it, of which there are infinitely many. This is to verify
% by running the above code twice.

% ---------------------------------
% ---------------------------------

% QS 2

A = [1,1,0;0,1,1;1,0,1];
B = [0;0;1];
C = [1,0,0];

T = [0.5774 + 0.0000i  -0.1852 + 0.5468i   0.5749 - 0.0531i;...
     0.5774 + 0.0000i   0.5662 - 0.1130i  -0.2414 + 0.5244i;...
     0.5774 + 0.0000i  -0.3810 - 0.4338i  -0.3335 - 0.4713i];

Tinv = transpose(conj(T));
Lambda = Tinv*A*T;
lambda1 = Lambda(1,1);
lambda2 = Lambda(2,2);
lambda3 = Lambda(3,3);

M1 = C*T;
M2 = Tinv*B;

M1(1)*M2(1);
M1(2)*M2(2);
M1(3)*M2(3);

% ---------------------------------
% ---------------------------------

% QS 5

A3 = [6,7,15;-1,-1,-3;-1,-1,-2];

% Find eigenvector of A3
% Pick an eigenvector and normalize it
% Use Gram-Schmidt to construct two orthnormal vectors
% Form a unitary matrix T
% T* A T
% Repeat process for 2x2 in lower-right corner
% Form TN = T*[1, TN-1]

% ---------------------------------
% Algorithm to derive the eigenvectors of the 3x3 matrix

syms lambda3
I3 = eye(3,3); % 3x3 identity matrix
LHS3 = A3 - lambda3*I3; % LHS of characteristic equation

char_poly3 = Det3(LHS3); % Compute the characteristic polynomial

lambda_vec3 = zeros(3,1); % Preallocate vec of eigenvalues
lambda_vec3 = lambda_vec3 + solve(char_poly3 == 0); % Vec of eigenvalues

Mata = A3 - lambda_vec3(1)*I3; % Substitute eigenvalue into LHS
Matb = rref(Mata); % Reduce to reduced row echelon form

% By inspection, [3; 0; -1] is an eigenvector

v1 = [3;0;-1];

% ---------------------------------

% Gram-Schmidt Orthnormalization for the 3x3 matrix

u1 = normalize(v1); % Normalize the eigenvector

v2 = [rand;rand;rand]; % Generate (2nd) random vector
u2_s = v2 - (transpose(conj(u1)))*v2*u1; % Subtract out part parallel to u1
u2 = normalize(u2_s); % Normalize

v3 = [rand;rand;rand]; % Generate (3rd) random vector
u3_s = v3 - (transpose(conj(u1)))*v3*u1 - (transpose(conj(u2)))*v3*u2; % Subtract out parts parallel to u1 & u2
u3 = normalize(u3_s); % Normalize

T3 = [u1 u2 u3];
rank(T3); % Verify that T3 is full rank
transpose(conj(T3))*T3; % Verify that T3 is unitary
L = transpose(conj(T3))*A3*T3; % Verify that bottom left element is zero

A2 = L([2 3], [2 3]);
% A([a b], [c d]): Extract the ath and bth rows, and the cth and dth columns

% ---------------------------------
% Algorithm to derive the eigenvectors of the 2x2 matrix

syms lambda2
I2 = eye(2,2); % 2x2 identity matrix
LHS2 = A2 - lambda2*I2; % LHS of characteristic equation

char_poly2 = Det2(LHS2); % Compute the characteristic polynomial

lambda_vec2 = zeros(2,1); % Preallocate vec of eigenvalues
lambda_vec2 = lambda_vec2 + solve(char_poly2 == 0); % Vec of eigenvalues

Mat1 = A2 - lambda_vec2(1)*I2; % Substitute eigenvalue into LHS
Mat2 = double(rref(Mat1)); % Reduce to reduced row echelon form

% ---------------------------------
% Gram-Schmidt Orthnormalization for the 2x2 matrix A2

v1 = [Mat2(1,2);-1];
u1 = normalize(v1); % Normalize the eigenvector

v2 = [rand;rand]; % Generate (2nd) random vector
u2_s = v2 - (transpose(conj(u1))*v2)*u1; % Subtract out part parallel to u1
u2 = normalize(u2_s); % Normalize

T2 = [u1 u2]; % T2, the matrix that diagonalizes A2
rank(T2); % Verify that T2 is full rank
transpose(conj(T2))*T2; % Verify that T2 is unitary
L2 = transpose(conj(T2))*A2*T2; % Verify that T2 diagonalizes A2

% ---------------------------------
% Constructing TN, the unitary matrix that diagonalizes A

RM = zeros(3,3);
RM(1,1) = 1;
RM([2 3], [2 3]) = T2;
T = T3*RM

transpose(conj(T))*T; % Verifying that TN is unitary
A_diag = transpose(conj(T))*A3*T % Verifying that TN makes A3 upper-diagonal

% ---------------------------------
% Solution

% TN = [0.95   0.17   0.27;...
%       0.00  -0.85   0.53;...
%      -0.32   0.51   0.80];

% A_diag = [1.00   2.67  17.24;...
%           0.00   1.00   4.43;...
%           0.00   0.00   1.00];
  
% ---------------------------------
% ---------------------------------
% Functions

function Det = Det3(X) % Calculates 3x3 Determinant

%A([a b], [c d]): Extract the ath and bth rows, and the cth and dth columns
X1 = X([2 3], [2 3]);
X2 = X([2 3], [1 3]);
X3 = X([2 3], [1 2]);

Det = X(1,1)*Det2(X1) - X(1,2)*Det2(X2) + X(1,3)*Det2(X3);
% 3x3 determinant in terms of 2x2 determinants

end

function Det = Det2(X) % Calculates 2x2 Determinant

Det = X(1,1)*X(2,2) - X(1,2)*X(2,1); % Leibniz Formula for 2x2 determinant

end

function x_norm = normalize(x) % Normalizes the input vector

x_norm = x/sqrt(transpose(conj(x))*x);

end

% ---------------------------------
% ---------------------------------


