% 逐步加密网格，计算固定频率的电荷
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
Charge=[];
SizeSchur=[];
for N=2:2:16
    [Q,SchurSize]=ComputeChargeExtra(N);
    Charge=[Charge;Q];
    SizeSchur=[SizeSchur;SchurSize];
    save('ChargeHCTNIDT8ConvOrderExtra.mat','Charge','SizeSchur')
    disp(['N=',num2str(N),'已完成'])
end

function [Q,SchurSize]=ComputeChargeExtra(N)
SchurSize=[];
omega=1.5*1e9;
N_HCT=8;
Nx1 = 4*N; Ny1 = 2; Nz1 = 2*N;
Nx2 = 2*N; Nz2 = N;
Nx_pml=N;
type="quadratic";
type_Extra=2;
% 计算粗网格解
switch type_Extra
    case 1
        % (0,0,0,0)
        [Q0000,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        disp("(0,0,0,0)已计算完成")
        % (1,0,0,0)
        [Q1000,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1000=Q1000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1000=Q1000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1000=Q1000+Q;
        % (1,1,0,0)
        [Q1100,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,2*Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        disp("(1,1,0,0)已计算完成")
        % (2,0,0,0)
        [Q2000,SizeSchurCopy] = SAWHCTUnbounded(omega,4*Nx1+1,Ny1,Nz1+1,4*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,4*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q2000=Q2000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,4*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q2000=Q2000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,4*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q2000=Q2000+Q;
        disp("(2,0,0,0)已计算完成")
    case 2
        % (0,0,0,0)
        [Q0000,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        disp("(0,0,0,0)已计算完成")
        % (1,0,0,0)
        [Q1000,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1000=Q1000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1000=Q1000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1000=Q1000+Q;
        % (1,1,0,0)
        [Q1100,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,2*Nz1+1,2*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,2*Nx1+1,Ny1,Nz1+1,2*Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,2*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,2*Nz1+1,Nx2+1,Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,2*Nz2+1,2*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q1100=Q1100+Q;
        disp("(1,1,0,0)已计算完成")
        % (2,0,0,0)
        [Q2000,SizeSchurCopy] = SAWHCTUnbounded(omega,3*Nx1+1,Ny1,Nz1+1,3*Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,3*Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q2000=Q2000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,3*Nz2+1,Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q2000=Q2000+Q;
        [Q,SizeSchurCopy] = SAWHCTUnbounded(omega,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,3*Nx_pml+1,N_HCT,type);
        SchurSize=[SchurSize,sum(SizeSchurCopy)];
        Q2000=Q2000+Q;
        disp("(2,0,0,0)已计算完成")
end
Q = ComputeExtraSol([Q0000,Q1000,Q1100,Q2000],4,2,type_Extra);
Q = sum(Q);
end