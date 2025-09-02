% 使用外推法计算表面电荷，s=4,m=2，四个自由维度分别为：
%   (1)基底、电极、PML的x方向
%   (2)基底z方向
%   (3)电极z方向
%   (4)PML深度方向
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
omega = 1.5*1e9;
N_HCT=8;
N=8;
Nx1 = 4*N; Ny1 = 2; Nz1 = 2*N;
Nx2 = 2*N; Nz2 = N;
Nx_pml=N;
type="quadratic";
type_Extra=2;
% 计算粗网格解
switch type_Extra
    case 1
        % (0,0,0,0)
        Q0000 = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        disp("(0,0,0,0)已计算完成")
        % (1,0,0,0)
        Q1000 = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q1000 = Q1000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q1000 = Q1000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        Q1000 = Q1000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        % (1,1,0,0)
        Q1100 = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,2*Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,2*Nx_pml+1,N_HCT,type);
        disp("(1,1,0,0)已计算完成")
        % (2,0,0,0)
        Q2000 = SAWHCTUnbounded(omega,4*Nx1+1,Ny1,Nz1+1,4*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q2000 = Q2000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,4*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q2000 = Q2000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,4*Nz2+1,Nx_pml+1,N_HCT,type);
        Q2000 = Q2000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,4*Nx_pml+1,N_HCT,type);
        disp("(2,0,0,0)已计算完成")
    case 2
        % (0,0,0,0)
        Q0000 = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        disp("(0,0,0,0)已计算完成")
        % (1,0,0,0)
        Q1000 = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q1000 = Q1000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q1000 = Q1000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        Q1000 = Q1000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        % (1,1,0,0)
        Q1100 = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,2*Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        Q1100 = Q1100 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,2*Nx_pml+1,N_HCT,type);
        disp("(1,1,0,0)已计算完成")
        % (2,0,0,0)
        Q2000 = SAWHCTUnbounded(omega,3*Nx1+1,Ny1,Nz1+1,3*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q2000 = Q2000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,3*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        Q2000 = Q2000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,3*Nz2+1,Nx_pml+1,N_HCT,type);
        Q2000 = Q2000 + SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,3*Nx_pml+1,N_HCT,type);
        disp("(2,0,0,0)已计算完成")
end
u_extra = ComputeExtraSol([Q0000,Q1000,Q1100,Q2000],4,2,type_Extra);
load('D:\Code\M\Mortar_FEM_Wavelet\NumericalEx\SAW\Data\N32ChargeIDT256.mat')
q_ref=sum(Q);
err=norm(u_extra-q_ref,inf)/norm(q_ref)
