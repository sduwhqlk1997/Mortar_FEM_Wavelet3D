clear
% clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
omega = 4.39822971502571*1e9;
N_HCT=5;
N=8;
Nx1 = 2*N+1; Ny1 = 2; Nz1 = 2*N+1;
Nx2 = N+1; Nz2 = N+1;
Nx_pml=N+1;
type="quadratic";
Q = SAWHCTUnbounded(omega,Nx1,Ny1,Nz1,Nx2,Nz2,Nx_pml,N_HCT,type);