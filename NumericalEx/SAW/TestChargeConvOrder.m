clear
% 测试表面电荷的收敛阶
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
omega = 1.5*1e9;
N_HCT=8;
N=18;
Nx1 = 4*N; Ny1 = 2; Nz1 = 2*N;
Nx2 = 2*N; Nz2 = N;
Nx_pml=N;
type="quadratic";
[Q,SizeSchur] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
q=sum(Q);
% save('result','Q','SizeSchur')
load('D:\Code\M\Mortar_FEM_Wavelet\NumericalEx\SAW\Data\N32ChargeIDT256.mat')
q_ref=sum(Q);
err=norm(q-q_ref,inf)/norm(q_ref)
