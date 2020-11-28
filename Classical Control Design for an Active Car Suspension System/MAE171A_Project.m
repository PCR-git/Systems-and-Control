% MAE 171A - Controller Design for an Active Suspension System
% Peter Racioppo

clear all;
close all;
clc;

%% System Parameters
% Masses [kg]: m1=2500, m2=320
% Spring constants [N/m]: k1=80,000, k2=500,000
% Damping coefficients [M/m/s]: b1=350, b2=15,000

m1 = 2500;
m2 = 320;
k1 = 80000;
k2 = 500000;
b1 = 350;
b2 = 15000;

g = 9.8; % Gravitational acceleration at Earth's surface

% Equilibrium positions
r1 = (((k1+k2)/(k1*k2))*m1 + m2/k2)*g;
r2 = (m1+m2)*g/k2;

%% Symbolic Calculations 1 (W = 0)

syms k1 k2 m1 m2 b1 b2 s F

A = s^2 + (k1+b1*s)*(1/m1+1/m2);
B = (k2+b2*s)/m2;
C = (s^2+(k2+b2*s)/m2);
D = (k1+b1*s)/m2;
E = F*(1/m1+1/m2);
Num = -B/m2 + C*(1/m1+1/m2);
Denom = A-B*D/C;
G = simplify(Num/Denom);

Num2 = (k2 + b2*s + m1*s^2 + m2*s^2);
Denom2 = (m1*m2*(s^2 - ((k1 + b1*s)*(k2 + b2*s))/(m2*(m2*s^2 + b2*s + k2)) + ((m1 + m2)*(k1 + b1*s))/(m1*m2))); 
collect(Num2,s);
collect(Denom2,s);

Num3 = ((m1 + m2)*s^2 + b2*s + k2)*(m2*s^2 + b2*s + k2);
collect(Num3,s);

clear all;

%% Symbolic Calculations 2 (F = 0)

syms k1 k2 m1 m2 b1 b2 s F

A = s^2 + (k1+b1*s)*(1/m1+1/m2);
B = (k2+b2*s)/m2;
C = (s^2+(k2+b2*s)/m2);
D = (k1+b1*s)/m2;

Num = B*(B-C);
Denom = A*C-B*D;
G = simplify(Num/Denom);

Num2 = -(s^2*(k2 + b2*s));
Denom2 = (m2*(((k2 + b2*s)/m2 + s^2)*((1/m1 + 1/m2)*(k1 + b1*s) + s^2) - ((k1 + b1*s)*(k2 + b2*s))/m2^2));

collect(Num2,s);
collect(Denom2,s);

Num_f = -b2*s^3 - k2*s^2;
Denom_f = m2*s^4 + m2*(b2/m2 + b1*(1/m1 + 1/m2))*s^3 + m2*(k2/m2 + k1*(1/m1 + 1/m2) - (b1*b2)/m2^2 + (b1*b2*(1/m1 + 1/m2))/m2)*s^2 - m2*((b1*k2)/m2^2 + (b2*k1)/m2^2 - (b1*k2*(1/m1 + 1/m2))/m2 - (b2*k1*(1/m1 + 1/m2))/m2)*s - m2*((k1*k2)/m2^2 - (k1*k2*(1/m1 + 1/m2))/m2);

clear all;

%% G1(s) (Input to Output Transfer Function)

nvec = [m2*(m1 + m2),(b2*m2 + b2*(m1 + m2)),(b2^2 + k2*(m1 + m2) + k2*m2),2*b2*k2,k2^2];
dvec = [m1*m2^2,(b1*m2^2 + b1*m1*m2 + b2*m1*m2),(k1*m2^2 + b1*b2*m2 + k1*m1*m2 + k2*m1*m2),(b1*k2*m2 + b2*k1*m2),k1*k2*m2];

Roots1 = roots(dvec);

sys1 = tf(nvec,dvec);

% Bode Plot
figure;
opts = bodeoptions('cstprefs');
opts.XLabel.FontSize = 30;
opts.YLabel.FontSize = 30;
opts.Title.String = '';
opts.TickLabel.FontSize = 30;
bode(sys1, opts, 'k');
grid on;

% Nyquist Plot
figure;
opts = nyquistoptions('cstprefs');
opts.XLabel.FontSize = 30;
opts.YLabel.FontSize = 30;
opts.Title.String = '';
opts.TickLabel.FontSize = 25;
nyquist(sys1, opts, 'k');
grid on;

% Impulse Response
figure;
impulse(sys1,'k');
title('', 'FontSize', 30);
xlabel('Time', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
set(gca,'FontSize',30);
grid on;

% Step Response
figure;
step(sys1,'k');
title('', 'FontSize', 30);
xlabel('Time', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
set(gca,'FontSize',30);
grid on;

% Root Locus
figure;
rlocus(sys1);
title('', 'FontSize', 30);
xlabel('Real Axis', 'FontSize', 30);
ylabel('Imaginary Axis', 'FontSize', 30);
set(gca,'FontSize',20);
grid on;

%% G2(s) (Noise to Output Transfer Function)

nvec = [-b2,-k2,0,0];
dvec = [m2,m2*(b2/m2 + b1*(1/m1 + 1/m2)),m2*(k2/m2 + k1*(1/m1 + 1/m2) - (b1*b2)/m2^2 + (b1*b2*(1/m1 + 1/m2))/m2),-m2*((b1*k2)/m2^2 + (b2*k1)/m2^2 - (b1*k2*(1/m1 + 1/m2))/m2 - (b2*k1*(1/m1 + 1/m2))/m2),-m2*((k1*k2)/m2^2 - (k1*k2*(1/m1 + 1/m2))/m2)];

Roots2 = roots(dvec);

sys2 = tf(nvec,dvec);

% Bode Plot
figure;
opts = bodeoptions('cstprefs');
opts.XLabel.FontSize = 30;
opts.YLabel.FontSize = 30;
opts.Title.String = '';
opts.TickLabel.FontSize = 30;
bode(sys2, opts, 'k');
grid on;

% Nyquist Plot
figure;
opts = nyquistoptions('cstprefs');
opts.XLabel.FontSize = 30;
opts.YLabel.FontSize = 30;
opts.Title.String = '';
opts.TickLabel.FontSize = 20;
nyquist(sys2, opts, 'k');
grid on;

% Nyquist Plot (one sided)
figure;
opts = nyquistoptions('cstprefs');
opts.XLabel.FontSize = 30;
opts.YLabel.FontSize = 30;
% opts.Xlim = [-2 0];
% opts.Ylim = [-5 5];
opts.Title.String = '';
opts.TickLabel.FontSize = 20;
opts.ShowFullContour = 'off';
nyquist(sys2, opts, 'k');
grid on;

% Impulse Response
figure;
impulse(sys2, 'k');
title('', 'FontSize', 30);
xlabel('Time', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
set(gca,'FontSize',30);
grid on;

% Step Response
figure;
step(sys2, 'k');
title('', 'FontSize', 30);
xlabel('Time', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
set(gca,'FontSize',30);
grid on;

%% Controller Design
% Specifications: we want a maximum overshoot of 5% and a 1 seconds
% settling time within 1% of the final value. Find the corresponding
% sector in the complex domain where the poles should be placed. Apply
% the formulas used for a second order system. They are not exact for
% this case, but they will work reasonably well since thereare dominant
% poles.

% Choose -5 +/- 2i
% Denominator: s^4 + 57.8892*s^3 + 2320.86*s^2 + 19518.5*s + 52576.1

% syms s
% num_t_sym = (m1+m2)*s^2/(m1*m2) + b2*s/(m1*m2)+ k2/(m1*m2);
% denom_t_sym = s^4 + 57.8892*s^3 + 2320.86*s^2 + 19518.5*s + 52576.1;
% vpa(simplify(num_t_sym/denom_t_sym),2);
% % (3.5e-3*s^2+0.019*s+0.62)/(s^4+59.0*s^3+2.3e3*s^2+2.0e4*s+5.3e4)

num_t = [(m1+m2)/(m1*m2), b2/(m1*m2), k2/(m1*m2)];
denom_t = [1,57.8892,2320.86,19518.5,52576.1];
sys_t = tf(num_t,denom_t);

% Step Response
figure;
step(sys_t, 'k');
title('', 'FontSize', 30);
xlabel('Time', 'FontSize', 30);
ylabel('Amplitude', 'FontSize', 30);
set(gca,'FontSize',30);
grid on;

% Sisotool / controlSystemDesigner(tf())

% controlSystemDesigner(sys1)
% sisotool(sys1)

% ----------------------------

% syms s
% % C = 1000*(s+1.8-6i)*(s+1.8+6i)/((s+5.5)*(s+14.5))
% % Cnum = 1000 s^2 + 3600 s + 39240
% % Cdenom = s^2 + 20 s + 79.75
% 
% Cnum = [1000,3600,39240];
% Cdnm = [1,20,79.75];
% sysC = tf(Cnum,Cdnm);
% 
% Cnum2 = [1000*0.16^2,91,1000];
% Cdnm2 = [0.01242,0.249,1];
% sysC2 = tf(Cnum/1000,Cdnm2);
% 
% sysL = sys1*sysC
% sysL2 = sys1*sysC2
% 
% % % Step Response
% % figure;
% % step(sysL2, 'k');
% % title('', 'FontSize', 30);
% % xlabel('Time', 'FontSize', 30);
% % ylabel('Amplitude', 'FontSize', 30);
% % set(gca,'FontSize',30);
% % grid on;

