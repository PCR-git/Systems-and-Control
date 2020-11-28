
load('Assign_9_2018_test.mat')
% load('Assignment_9_Example_2.mat')

[sysKP, sysKPerr, sysKF, sysKFerr] = assign9_PeterRacioppo(Plant, W);

figure;
bode(sysKP);
figure;
bode(sysKPerr);
figure;
bode(sysKF);
figure;
bode(sysKFerr);

clc
max(max(sysKF.A-sysKF_a9sg.A))
max(max(sysKFerr.A-sysKFerr_a9sg.A))
max(max(sysKP.A-sysKP_a9sg.A))
max(max(sysKPerr.A-sysKPerr_a9sg.A))
max(max(sysKF.B-sysKF_a9sg.B))
max(max(sysKFerr.B-sysKFerr_a9sg.B))
max(max(sysKP.B-sysKP_a9sg.B))
max(max(sysKPerr.B-sysKPerr_a9sg.B))
max(max(sysKF.C-sysKF_a9sg.C))
max(max(sysKFerr.C-sysKFerr_a9sg.C))
max(max(sysKP.C-sysKP_a9sg.C))
max(max(sysKPerr.C-sysKPerr_a9sg.C))
max(max(sysKF.D-sysKF_a9sg.D))
max(max(sysKFerr.D-sysKFerr_a9sg.D))
max(max(sysKP.D-sysKP_a9sg.D))
max(max(sysKPerr.D-sysKPerr_a9sg.D))
