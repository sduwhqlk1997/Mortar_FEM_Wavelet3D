%扫频计算表面电荷
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
N_HCT=8;
N=2;
Nx1 = 4*N; Ny1 = 2; Nz1 = 2*N;
Nx2 = 2*N; Nz2 = N;
Nx_pml=N;
type="quadratic";
Charge=zeros(2^8,1500-1400);
stat=1;
for i=1400:1500
    Charge(:,stat)=SAWHCTUnbounded(i*1e6,Nx1,Ny1,Nz1,Nx2,Nz2,Nx_pml,N_HCT,type);
    save('ChargeHCTOnlyNIDT8N2.mat','Charge')
    disp(['已完成',num2str(stat),'个频点'])
    stat=stat+1;
end