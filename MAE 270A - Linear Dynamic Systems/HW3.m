
% Peter Racioppo
% 103953689
% Linear Systems - HW3

clear all;
clc;

% Defining Matrices
M1 = [1,1;2,1];
M2 = [1,1;1,2];
M3 = [-1,1;-1,-2];
M4 = [1,1;1,1];

% Singular Value Decompositions
[U1,S1,V1] = svd(M1);
[U2,S2,V2] = svd(M2);
[U3,S3,V3] = svd(M3);
[U4,S4,V4] = svd(M4);

% Plotting the Ellipses
f_ellipse(S1,V1,U1)
f_ellipse(S2,V2,U2)
f_ellipse(S3,V3,U3)
f_ellipse(S4,V4,U4)

% ----------------------------------------------------------

% Computes and Plots Ellipses and Input/Output Vectors
function f_ellipse(S,V,U)

a1 = max(diag(S)); % Major axis length
b1 = min(diag(S)); % Minor axis length
high_gain_in_vec1 = V(:,1); % High Gain Input Vector
low_gain_in_vec1  = V(:,2); % Low Gain Input Vector
high_gain_out_vec1 = U(:,1)*a1; % High Gain Output Vector
low_gain_out_vec1  = U(:,2)*b1; % Low Gain Output Vector

% Angle of major axis with x-axis
ang1 = -atan2(high_gain_out_vec1(2),high_gain_out_vec1(1));
% Rotation Matrix
Rot1 = [cos(ang1), sin(ang1);...
        -sin(ang1),cos(ang1)];

% Parameterizing Circle
n = 1000;
t_circ_1 = linspace(0,2*pi,n);
x_circ_1 = cos(t_circ_1);
y_circ_1 = sin(t_circ_1);

% Parameterizing Ellipse
t_ell_1 = linspace(0,2*pi,n);
x_ell_1 = a1*cos(t_ell_1);
y_ell_1 = b1*sin(t_ell_1);
ell_vec_1 = [x_ell_1;y_ell_1];
% Rotated Ellipse
ell_rot_1 = Rot1*ell_vec_1;
x_ell_rot_1 = ell_rot_1(1,:);
y_ell_rot_1 = ell_rot_1(2,:);

% Plot
figure;

plot(x_circ_1,y_circ_1,'b-'); % Plot Circle
hold on;
plot(x_ell_rot_1,y_ell_rot_1,'r-'); % Plot Ellipse

% Inputs
plot([0;high_gain_in_vec1(1)],[0;high_gain_in_vec1(2)],'r--');
plot(high_gain_in_vec1(1),high_gain_in_vec1(2),'ro','MarkerFaceColor', 'r');
plot([0;low_gain_in_vec1(1)],[0;low_gain_in_vec1(2)],'k--');
plot(low_gain_in_vec1(1),low_gain_in_vec1(2),'ko','MarkerFaceColor', 'k');

% Outputs
plot([0;high_gain_out_vec1(1)],[0;high_gain_out_vec1(2)],'r-');
plot(high_gain_out_vec1(1),high_gain_out_vec1(2),'ro','MarkerFaceColor', 'r');
plot([0;low_gain_out_vec1(1)],[0;low_gain_out_vec1(2)],'k-');
plot(low_gain_out_vec1(1),low_gain_out_vec1(2),'ko','MarkerFaceColor', 'k');

grid on;
axis equal;

end

% ----------------------------------------------------------

% The input and output vectors of the Hermitian matrices are aligned.
% The final matrix has rank 1 and its output ellipse collapses to a line.
