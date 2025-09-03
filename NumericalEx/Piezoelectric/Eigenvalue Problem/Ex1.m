clear
% clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 材料参数
ModelCoeff = 'D:\Code\M\Mortar_FEM_Wavelet\Piezoelectric\Data\ModelCoef2.mat';
%% 数值求解
% equ_type = "saddle";
equ_type = "schur";
% 导入材料参数
load(ModelCoeff)
c_LN = cell2mat(materials(2));
e_LN = cell2mat(materials(5));
epcl0 = cell2mat(materials(6));
epcl_LN = cell2mat(materials(7));
epcl_LN = epcl_LN*epcl0; 
rho = cell2mat(materials(8));
% 几何参数
ori = [0,0,0];
a = 1; b = 1; c = 1;
kappa_bar_sq=1;
% 无量纲化
[c_LN, e_LN, epcl_LN, L_bar, L0, vs, Phi0, kappa_bar_sq] = ...
         piezo_dimensionless(c_LN, e_LN, epcl_LN, rho, a, b, c); rho=1;
% 离散
type="quadratic";
N = 8;
Nx=N+1;Ny=N+1;Nz=N+1;
[K,M,~,Dof_Index] =...
    AssemblePiezMatFEM(c_LN,e_LN,epcl_LN,rho,kappa_bar_sq,...
    ori,a,b,c,Nx,Ny,Nz,...
    [],[],[],type);
% 处理Dirichlete边界条件
boundary_u=find(Dof_Index(:,1)==1&...
    (Dof_Index(:,2)==ori(1)|Dof_Index(:,2)==ori(1)+a|...
    Dof_Index(:,3)==ori(2)|Dof_Index(:,3)==ori(2)+b|...
    Dof_Index(:,4)==ori(3)|Dof_Index(:,4)==ori(3)+c));
boundary_v=find(Dof_Index(:,1)==2&...
    (Dof_Index(:,2)==ori(1)|Dof_Index(:,2)==ori(1)+a|...
    Dof_Index(:,3)==ori(2)|Dof_Index(:,3)==ori(2)+b|...
    Dof_Index(:,4)==ori(3)|Dof_Index(:,4)==ori(3)+c));
boundary_w=find(Dof_Index(:,1)==3&...
    (Dof_Index(:,2)==ori(1)|Dof_Index(:,2)==ori(1)+a|...
    Dof_Index(:,3)==ori(2)|Dof_Index(:,3)==ori(2)+b|...
    Dof_Index(:,4)==ori(3)|Dof_Index(:,4)==ori(3)+c));
boundary_phi=find(Dof_Index(:,1)==4&...
    (Dof_Index(:,2)==ori(1)|Dof_Index(:,2)==ori(1)+a|...
    Dof_Index(:,3)==ori(2)|Dof_Index(:,3)==ori(2)+b|...
    Dof_Index(:,4)==ori(3)|Dof_Index(:,4)==ori(3)+c));
K([boundary_u;boundary_v;boundary_w;boundary_phi],:) = [];
K(:,[boundary_u;boundary_v;boundary_w;boundary_phi]) = [];
M([boundary_u;boundary_v;boundary_w;boundary_phi],:) = [];
M(:,[boundary_u;boundary_v;boundary_w;boundary_phi]) = [];
Dof_Index([boundary_u;boundary_v;boundary_w;boundary_phi],:)=[];
%% 计算特征值
if equ_type == "saddle"
% 直接计算鞍点系统
% [V_saddle,lambda_saddle]=eigs(K,M,10,'sm');
[lambda_saddle,V_saddle]=inverse_iteration_GPT(K, M, 0,1E-10,1000);
[lambda_saddle, V_saddle, phi] = piezo_restore_dimension(lambda_saddle, V_saddle, L0, vs, Phi0);
elseif equ_type == "schur"
% 计算Schur补系统
index_Pot = Dof_Index(:,1)==4;
index_disp = ~index_Pot;
S = K(index_Pot,index_Pot)\K(index_Pot,index_disp);
S = K(index_disp,index_disp)-K(index_disp,index_Pot)*S;
M_S = M(index_disp,index_disp);
[V_schur,lambda_schur]=eigs(S,M_S,10,'sm');
[lambda_schur, V_schur, phi] = piezo_restore_dimension(lambda_schur, V_schur, L0, vs, Phi0);
end





function [lambda, v] = inverse_iteration_GPT(A, B, sigma, tol, maxit)
% 求解广义特征值问题 A u = lambda B u
% 使用反幂(shift-invert)法, 收敛到最靠近 sigma 的特征值
%
% 输入:
%   A, B   - 矩阵 (n x n), B SPD
%   sigma  - 移位参数
%   tol    - 收敛容忍度 (如 1e-10)
%   maxit  - 最大迭代次数
%
% 输出:
%   lambda - 近似特征值
%   v      - 对应特征向量 (B-单位化)

    if nargin < 5, maxit = 200; end
    if nargin < 4, tol = 1e-10; end
    n = size(A,1);

    % 初始向量
    v = randn(n,1);
    %v = v / sqrt(v'*B*v);  % B-归一化
    v = v / sqrt(v'*v);

    % 因式分解 A - sigma B
    [L,U,P,Q] = lu(A - sigma*B);

    lambda_old = inf;
    for k = 1:maxit
        % 1) 解线性方程 (A - sigma B) w = B v
        rhs = B*v;
        w = Q * (U \ (L \ (P * rhs)));   % 使用LU分解求解
        
        % 2) B-单位化
        % lambda=max(w);
        v=w/max(abs(w));
        % v = w / sqrt(v' * B * w);
        % v = w / sqrt(v' * w);

        % 3) Rayleigh 商 (近似特征值)
        lambda = (v' * A * v) / (v' * B * v);

        % 4) 收敛判据
        res = norm(A*v - lambda*B*v);
        if res < tol %|| abs(lambda - lambda_old) < tol*(1+abs(lambda))
            fprintf('Converged in %d iterations, residual = %.2e\n', k, res);
            return;
        end
        lambda_old = lambda;
    end

    warning('未收敛，最大迭代数 = %d, 最后残差 = %.2e', maxit, res);
end
