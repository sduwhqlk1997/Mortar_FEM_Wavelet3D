% 有精确解的特征值问题算例
clear
% clc
% currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
currentPath = 'D:\code\Mortar_FEM_Wavelet3D';
addpath(genpath(currentPath));
% 材料参数
% ========== 压电材料参数 (PZT-5A 量级) ==========
% 所有参数使用 Voigt 表示法

% 密度 [kg/m³]
rho = 7500;

% 弹性矩阵 (6×6, 单位: Pa)
% Voigt 顺序: 11, 22, 33, 23, 13, 12
c_LN = 1e9 * [ % 乘以 10^9 转换为 Pa
    120,  75,  75,   0,   0,   0;
     75, 120,  75,   0,   0,   0;
     75,  75, 120,   0,   0,   0;
      0,   0,   0, 22.5,  0,   0;
      0,   0,   0,   0, 22.5,  0;
      0,   0,   0,   0,   0, 22.5
];

% 压电矩阵 (3×6, 单位: C/m²)
% 行: 电场方向 (1=x, 2=y, 3=z)
% 列: Voigt 力学索引
e_LN = [
     0,   0,   0,   0, 12.3,  0;
     0,   0,   0, 12.3,   0,  0;
  -5.4,-5.4,15.8,   0,   0,  0
];

% 介电矩阵 (3×3, 单位: F/m)
epcl_LN = 1e-9 * [ % 乘以 10^{-9} 
    8.85, 0,    0;
    0,    8.85, 0;
    0,    0,    8.85
];

% ========== 最小特征值精确解参数 ==========
% 特征频率 [rad/s]
omega_exact = pi * sqrt(c_LN(5,5) / rho); % ≈ 5442 rad/s

% 位移场函数 (x,y,z 为坐标向量)
u_x = @(x,y,z) sin(pi*z);
u_y = @(x,y,z) 0;
u_z = @(x,y,z) 0;

% 电势场函数
phi = @(x,y,z) 0;

% 几何参数
ori = [0,0,0];
a = 1; b = 1; c = 1;
kappa_bar_sq=1;
% 离散
type="quadratic";
N = 2;
Nx=N+1;Ny=N+1;Nz=N+1;
[K,M,~,Dof_Index] =...
    AssemblePiezMatFEM(c_LN,e_LN,epcl_LN,rho,kappa_bar_sq,...
    ori,a,b,c,Nx,Ny,Nz,...
    [],[],[],type);
% 边界条件
bound_vanish=find(Dof_Index(:,4)==ori(3)|Dof_Index(:,4)==ori(3)+c);
bound_xL = find(Dof_Index(:,2)==ori(1));
bound_xR = find(Dof_Index(:,2)==ori(1)+a);
bound_yL = find(Dof_Index(:,3)==ori(2));
bound_yR = find(Dof_Index(:,3)==ori(2)+b);

K(bound_xL,:)=K(bound_xL,:)+K(bound_xR,:);
K(bound_yL,:)=K(bound_yL,:)+K(bound_yR,:);
K(:,bound_xL)=K(:,bound_xL)+K(:,bound_xR);
K(:,bound_yL)=K(:,bound_yL)+K(:,bound_yR);

M(bound_xL,:)=M(bound_xL,:)+M(bound_xR,:);
M(bound_yL,:)=M(bound_yL,:)+M(bound_yR,:);
M(:,bound_xL)=M(:,bound_xL)+M(:,bound_xR);
M(:,bound_yL)=M(:,bound_yL)+M(:,bound_yR);

bound_vanish=[bound_vanish;bound_xR;bound_yR];

K(bound_vanish,:)=[];
K(:,bound_vanish)=[];
M(bound_vanish,:)=[];
M(:,bound_vanish)=[];
Dof_Index(bound_vanish,:)=[];

%% 计算特征值
[V_saddle_eigs,lambda_saddle_eigs]=eigs(K,M,1,'sm');
err = abs(lambda_saddle_eigs-omega_exact^2)/omega_exact^2;