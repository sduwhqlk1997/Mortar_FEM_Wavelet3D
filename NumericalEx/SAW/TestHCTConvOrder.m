clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
N_HCT=8;
Charge=[];
SizeSchur=[];
type="quadratic";
for N=10:2:16
    Nx1 = 4*N; Ny1 = 2; Nz1 = 2*N;
    Nx2 = 2*N; Nz2 = N;
    Nx_pml=N;
    [Q,SizeSchurCopy]=SAWHCTUnbounded(1.5*1e9,Nx1+1,Ny1,Nz1+1,Nx2+1,Nz2+1,Nx_pml+1,N_HCT,type);
    Charge=[Charge;sum(Q)];
    SizeSchur=[SizeSchur;sum(SizeSchurCopy)];
    save('ChargeHCTNIDT8ConvOrder.mat','Charge','SizeSchur')
    disp(['N=',num2str(N),'已完成'])
end
