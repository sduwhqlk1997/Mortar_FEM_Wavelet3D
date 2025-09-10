clear
% clc
% currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
currentPath = 'D:\code\Mortar_FEM_Wavelet3D';
addpath(genpath(currentPath));
% 材料参数
% ModelCoeff = 'D:\Code\M\Mortar_FEM_Wavelet\Piezoelectric\Data\ModelCoef2.mat';
ModelCoeff = 'ModelCoef2.mat';
%% 数值求解
% 求解的方程形式
equ_type = "saddle";
% equ_type = "schur";
% 使用的数值方法
method_type="power";
% method_type="J-D";
% 无量纲化？
dimensionless=true;
% dimensionless=false;
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
if dimensionless==true
[c_LN, e_LN, epcl_LN, L_bar, L0, vs, Phi0, kappa_bar_sq] = ...
         piezo_dimensionless(c_LN, e_LN, epcl_LN, rho, a, b, c); rho=1;
end
% 离散
type="quadratic";
N = 2;
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
[V_saddle_eigs,lambda_saddle_eigs]=eigs(K,M,1,'sm');
if method_type == "power"
[lambda_saddle,V_saddle]=inverse_iteration_GPT(K, M, 0,1E-10,1000);
elseif method_type =="J-D"
[lambda_saddle,V_saddle]=JD_iteration_Gen(K,[],M,1E-10,200);
end
if dimensionless == true
    [lambda_saddle_eigs, ~, ~] = piezo_restore_dimension(lambda_saddle_eigs, V_saddle_eigs, L0, vs, Phi0);
    [lambda_saddle, ~, ~] = piezo_restore_dimension(lambda_saddle, V_saddle, L0, vs, Phi0);
end
rel_err=abs(lambda_saddle-lambda_saddle_eigs)/abs(lambda_saddle_eigs)
elseif equ_type == "schur"
% 计算Schur补系统
index_Pot = Dof_Index(:,1)==4;
index_disp = ~index_Pot;
S = K(index_Pot,index_Pot)\K(index_Pot,index_disp);
S = K(index_disp,index_disp)-K(index_disp,index_Pot)*S;
M_S = M(index_disp,index_disp);
[V_saddle_eigs,lambda_saddle_eigs]=eigs(S,M_S,1,'sm');
if method_type == "power"
elseif method_type =="J-D"
    [lambda_saddle,V_saddle,u]=JD_iteration(S,M_S,1E-10,1000);
end
if dimensionless == true
    [lambda_saddle,~, ~] = piezo_restore_dimension(lambda_saddle,ones(size(K,1),1), L0, vs, Phi0);
    [lambda_saddle_eigs, ~, ~] = piezo_restore_dimension(lambda_saddle_eigs,ones(size(K,1),1), L0, vs, Phi0);
end
rel_err=abs(lambda_saddle-lambda_saddle_eigs)/abs(lambda_saddle_eigs)
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
    % v = randn(n,1);
    v=ones(n,1);
    v = v / sqrt(v'*B*v);  % B-归一化
    % v = v / sqrt(v'*v);

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
        if res < tol || abs(lambda - lambda_old) < tol*(1+abs(lambda))
            fprintf('Converged in %d iterations,\n residual = %.2e\n, distance between two adjoint step is %.2e', k, res,abs(lambda - lambda_old)/(1+abs(lambda)));
            return;
        end
        lambda_old = lambda;
    end

    warning('未收敛，最大迭代数 = %d, 最后残差 = %.2e', maxit, res);
end

function [lambda,vec]=JD_iteration_Gen(A,PA,B,tol, maxit)
if ~exist('PA','var')||isempty(PA), PA=A; end
n = size(A,1);
% v = randn(n,1);
v = ones(n,1);
v = v/norm(v);
w=A*v;
w_tilde=B*v;
vAv=v'*w; % Rayleigh商
vBv=v'*w_tilde;
u=v; theta = vAv/vBv;
% r = w-theta*w_tilde;
for k=1:maxit
    M = PA-theta*B;
    % M = spdiags(diag(M), 0, size(M,1), size(M,2));
    % PM=ilu(M);
    Mu=B*u; Mu = M\Mu;
    % Mu=gmres(M,Mu,[],1e-4,30,PM);
    epcl=(u'*u)/(u'*Mu);
    t=-u+epcl*Mu;
    t = modifiedGramSchmidt(t, v);
    t = t/norm(t);
    v=[v,t];
    w=[w,A*t];
    vAv=v'*A*v;
    vBv=v'*B*v;
    % h=[h;t'*w(:,1:k)];
    % h=[h,v'*w(:,k+1)];
    % h(end,1:end-1)=[h;h(1:)]
    [s,theta]=eig(vAv,vBv,'vector');
    [~, idx] = min((theta)); % 以最小模特征值为例
    theta = theta(idx);        % 当前近似特征值
    s = s(:, idx);
    s = s/norm(s);
    u=v*s;
    uhat=w*s;
    r=uhat-theta*B*u;
    % r=A*u-theta*u;
    res=norm(r);
    if res<tol
        lambda=theta;
        vec=u;
        fprintf('Converged in %d iterations, residual = %.2e\n', k, res);
        return;
    end
end
fprintf('Reached the maximum number of %d iterations, residual = %.2e\n', maxit, res);
lambda=theta;
vec=u;
function t_orth = modifiedGramSchmidt(t, V)
% 函数：modifiedGramSchmidt
% 输入：
%   t - 待正交化的向量（列向量）
%   V - 标准正交基组成的矩阵（每列是一个基向量）
% 输出：
%   t_orth - 正交化后的向量，满足与V的每一列正交

% 检查输入维度
[m, n] = size(V);
if size(t, 1) ~= m
    error('向量t的维度和矩阵V的行数不一致');
end
if size(t, 2) ~= 1
    error('向量t必须是列向量');
end

% 初始化正交化结果
t_orth = t;

% Modified Gram-Schmidt过程
for j = 1:n
    vj = V(:, j);          % 取出第j个标准正交基
    proj_coeff = vj' * t_orth; % 计算投影系数
    t_orth = t_orth - proj_coeff * vj; % 减去投影分量
end
end
end