function [K00,K11,K12,K22,K13,K23,K33,...
    Dof_index,P,T,Block_interval] = ...
    CombineSomeWaveDomain3D(Block,...
    K00_ref,K11_ref,K12_ref,K22_ref,K13_ref,K23_ref,K33_ref,...
    DofIndex_sparse,BaseFunIndex_ref)
%CombineSomeWaveDomain3D 按Block的分块组装整个区域上的基本矩阵，
% 在整个区域上要求方程系数为常数（注：目前只能在同一方向进行区域分解）
%   Block:区域分块，第一个元素表示x轴的剖分，第二个元素表示y轴的剖分,为一个cell
%   Kij_ref：参考刚度矩阵，在[0,1]^3参考单元上组装的基本刚度矩阵
%   DofIndex_sparse:稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号
%   BaseFunIndex_ref：[0,1]区间上的小波基函数的信息，行号为基函数编号，1列为type2，2~5列分别表示基函数的j,k,支集，6列为坐标索引

% 根据Block进行区域剖分
X=cell2mat(Block(1)); Y=cell2mat(Block(2)); Z=cell2mat(Block(3));
Nx=length(X);Ny=length(Y);Nz=length(Z);
[X_mesh,Y_mesh,Z_mesh] = meshgrid(X,Y,Z);
P = [X_mesh(:),Y_mesh(:),Z_mesh(:)];
index = 1:Nx*Ny*Nz;

index = reshape(index,Ny,Nx,Nz);
T = [reshape(index(1:Ny-1,1:Nx-1,1:Nz-1),[],1),...
    reshape(index(2:Ny,1:Nx-1,1:Nz-1),[],1),...
    reshape(index(1:Ny-1,2:Nx,1:Nz-1),[],1),...
    reshape(index(2:Ny,2:Nx,1:Nz-1),[],1),...
    reshape(index(1:Ny-1,1:Nx-1,2:Nz),[],1),...
    reshape(index(2:Ny,1:Nx-1,2:Nz),[],1),...
    reshape(index(1:Ny-1,2:Nx,2:Nz),[],1),...
    reshape(index(2:Ny,2:Nx,2:Nz),[],1)];
% Block_interval每行表示一个区域，1~2个元素表示该区域x取值范围，3~4表示y取值范围
Block_interval=[min(reshape(P(T,1),size(T,1),size(T,2)),[],2),max(reshape(P(T,1),size(T,1),size(T,2)),[],2),...
    min(reshape(P(T,2),size(T,1),size(T,2)),[],2),max(reshape(P(T,2),size(T,1),size(T,2)),[],2),...
    min(reshape(P(T,3),size(T,1),size(T,2)),[],2),max(reshape(P(T,3),size(T,1),size(T,2)),[],2)];

Jacobi=(Block_interval(:,2)-Block_interval(:,1)).*...
    (Block_interval(:,4)-Block_interval(:,3)).*...
    (Block_interval(:,6)-Block_interval(:,5));
Jacobi_inv=[1./(Block_interval(:,2)-Block_interval(:,1)),...
    1./(Block_interval(:,4)-Block_interval(:,3)),...
    1./(Block_interval(:,6)-Block_interval(:,5))];
% 第一个区域对应地矩阵
K00=Jacobi(1)*K00_ref;
K11=Jacobi(1)*Jacobi_inv(1,1)^2*K11_ref;
K12=Jacobi(1)*Jacobi_inv(1,1)*Jacobi_inv(1,2)*K12_ref;
K22=Jacobi(1)*Jacobi_inv(1,2)^2*K22_ref;
K13=Jacobi(1)*Jacobi_inv(1,1)*Jacobi_inv(1,3)*K13_ref;
K23=Jacobi(1)*Jacobi_inv(1,2)*Jacobi_inv(1,3)*K23_ref;
K33=Jacobi(1)*Jacobi_inv(1,3)^2*K33_ref;

BaseFunIndex_x=BaseFunIndex_ref;
BaseFunIndex_y=BaseFunIndex_ref;
BaseFunIndex_z=BaseFunIndex_ref;

BaseFunIndex_x(:,4:6)=BaseFunIndex_x(:,4:6)*(Block_interval(1,2)-Block_interval(1,1))+Block_interval(1,1);
BaseFunIndex_y(:,4:6)=BaseFunIndex_y(:,4:6)*(Block_interval(1,4)-Block_interval(1,3))+Block_interval(1,3);
BaseFunIndex_z(:,4:6)=BaseFunIndex_z(:,4:6)*(Block_interval(1,6)-Block_interval(1,5))+Block_interval(1,5);

Dof_index = GenWEMBaseIndex3D([Block_interval(1,1:2);Block_interval(1,3:4);Block_interval(1,5:6)],...
    DofIndex_sparse,BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z);
if size(T,1)>1 % 若多于一个单元
    for i=2:size(T,1) % 逐单元Matching
        BaseFunIndex_x=BaseFunIndex_ref;
        BaseFunIndex_y=BaseFunIndex_ref;
        BaseFunIndex_z=BaseFunIndex_ref;
        BaseFunIndex_x(:,4:6)=BaseFunIndex_x(:,4:6)*(Block_interval(i,2)-Block_interval(i,1))+Block_interval(i,1);
        BaseFunIndex_y(:,4:6)=BaseFunIndex_y(:,4:6)*(Block_interval(i,4)-Block_interval(i,3))+Block_interval(i,3);
        BaseFunIndex_z(:,4:6)=BaseFunIndex_z(:,4:6)*(Block_interval(i,6)-Block_interval(i,5))+Block_interval(i,5);
        Dof_index2 = GenWaveletBaseIndex3D([Block_interval(i,1:2);Block_interval(i,3:4);Block_interval(i,5:6)],...
            DofIndex_sparse,BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z);

        [K00,Dof_index1] = MatchingTwoCubeWaveDomain(K00,Jacobi(i)*K00_ref,Dof_index,Dof_index2,1);
        [K11,~] = MatchingTwoCubeWaveDomain(K11,Jacobi(i)*Jacobi_inv(i,1)^2*K11_ref,Dof_index,Dof_index2,1);
        [K12,~] = MatchingTwoCubeWaveDomain(K12,Jacobi(i)*Jacobi_inv(i,1)*Jacobi_inv(i,2)*K12_ref,Dof_index,Dof_index2,1);
        [K22,~] = MatchingTwoCubeWaveDomain(K22,Jacobi(i)*Jacobi_inv(i,2)^2*K22_ref,Dof_index,Dof_index2,1);
        [K13,~] = MatchingTwoCubeWaveDomain(K13,Jacobi(i)*Jacobi_inv(i,1)*Jacobi_inv(i,3)*K13_ref,Dof_index,Dof_index2,1);
        [K23,~] = MatchingTwoCubeWaveDomain(K23,Jacobi(i)*Jacobi_inv(i,2)*Jacobi_inv(i,3)*K23_ref,Dof_index,Dof_index2,1);
        [K33,~] = MatchingTwoCubeWaveDomain(K33,Jacobi(i)*Jacobi_inv(i,3)^2*K33_ref,Dof_index,Dof_index2,1);
        Dof_index=Dof_index1;
    end
end
end

