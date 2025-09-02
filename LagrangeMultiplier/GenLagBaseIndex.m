function [LagBaseIndex,LagMeshIndex] = GenLagBaseIndex(baseType,x_step,y_step)
%GENLAGBASEINDEX 生成Lagrange乘子基函数信息表
%   interval：区域，为2*2的矩阵，每一行一个维度
%   baseType：有限元网格
%   Nunknows：未知量个数
%   x_step,y_step：x，y方向的网格。注：不是有限元网格
%   LagBaseIndex：Lagrange乘子基函数信息表，每一行表示一个基函数，各列存储的信息如下所示：
%       1~4：所属单元，1~2第一个维度，3~4第二个维度
%       5~6：两个维度对应的一维基函数编号
%       7~8：标记两个维度是否为边界单元，若为0则为内部单元，若为1则为左侧单元，若为2则为右侧单元
%   LagMeshIndex：Lagrange乘子网格索引，每一行表示一个网格
%       1~2：网格第一个维度的区间，3~4：网格第二个维度的区间
%       剩下元素为对的基函数编号，即对应于LagBaseIndex的某些行,若为0则不表示任何基函数

switch baseType
    case "linear"
        T_x=[x_step(1:end-1),x_step(2:end)];
        T_y=[y_step(1:end-1),y_step(2:end)];
        Nx=size(T_x,1);
        Ny=size(T_y,1);
        %{
        LagBaseIndex1D_x：一维基函数信息表，每一行表示一个基函数，每一列存储的信息为：
        1~2：所属单元的左右端点；3：对应地一维基函数编号；4：标记是否为内部单元，
        若为1则为左侧单元，若为2则为右侧单元。
        %}
        LagBaseIndex1D_x=zeros(2*Nx,4);
        LagBaseIndex1D_y=zeros(2*Ny,4);
        for i=1:Nx
            LagBaseIndex1D_x(2*i-1,:)=[T_x(i,:),1,0];
            LagBaseIndex1D_x(2*i,:)=[T_x(i,:),2,0];
        end
        for i=1:Ny
            LagBaseIndex1D_y(2*i-1,:)=[T_y(i,:),1,0];
            LagBaseIndex1D_y(2*i,:)=[T_y(i,:),2,0];
        end
        % 分别合并左右两端点的基函数
        LagBaseIndex1D_x(1,:)=[];
        LagBaseIndex1D_x(end,:)=[];
        LagBaseIndex1D_x(1,4)=1;
        LagBaseIndex1D_x(end,4)=2;

        LagBaseIndex1D_y(1,:)=[];
        LagBaseIndex1D_y(end,:)=[];
        LagBaseIndex1D_y(1,4)=1;
        LagBaseIndex1D_y(end,4)=2;
        % 使用张量积生成两维信息表
        % 基函数信息
        Nx=size(LagBaseIndex1D_x,1);
        Ny=size(LagBaseIndex1D_y,1);
        index1=(1:Nx)';
        index2=(1:Ny)';
        [index1,index2]=meshgrid(index1,index2);
        index1=index1(:);
        index2=index2(:);
        LagBaseIndex=zeros(size(index1,1),8);
        LagBaseIndex(:,[1,2,5,7])=LagBaseIndex1D_x(index1,:);
        LagBaseIndex(:,[3,4,6,8])=LagBaseIndex1D_y(index2,:);
        % 网格信息
        Nx=size(T_x,1);
        Ny=size(T_y,1);
        index1=(1:Nx)';
        index2=(1:Ny)';
        [index1,index2]=meshgrid(index1,index2);
        index1=index1(:);
        index2=index2(:);
        LagMeshIndex=zeros(size(index1,1),8);
        LagMeshIndex(:,1:2)=T_x(index1,:);
        LagMeshIndex(:,3:4)=T_y(index2,:);
        for i=1:size(index1,1)
            indexBase=find(all(LagBaseIndex(:,1:4)==LagMeshIndex(i,1:4),2));
            for j=1:length(indexBase)
                LagMeshIndex(i,4+j)=indexBase(j);
            end
        end
    case "quadratic"
        T_x=[x_step(1:end-1),x_step(2:end)];
        T_y=[y_step(1:end-1),y_step(2:end)];
        Nx=size(T_x,1);
        Ny=size(T_y,1);
        %{
        LagBaseIndex1D_x：一维基函数信息表，每一行表示一个基函数，每一列存储的信息为：
        1~2：所属单元的左右端点；3：对应地一维基函数编号；4：标记是否为内部节点，
        若为1则为左端点，若为2则为右端点。
        %}
        LagBaseIndex1D_x=zeros(3*Nx,4);
        LagBaseIndex1D_y=zeros(3*Ny,4);
        for i=1:Nx
            LagBaseIndex1D_x(3*i-2,:)=[T_x(i,:),1,0]; % 左
            LagBaseIndex1D_x(3*i-1,:)=[T_x(i,:),3,0]; % 中
            LagBaseIndex1D_x(3*i,:)=[T_x(i,:),2,0]; % 右
        end
        for i=1:Ny
            LagBaseIndex1D_y(3*i-2,:)=[T_y(i,:),1,0]; % 左
            LagBaseIndex1D_y(3*i-1,:)=[T_y(i,:),3,0]; % 中
            LagBaseIndex1D_y(3*i,:)=[T_y(i,:),2,0]; % 右
        end
        % 分别合并左右两端点的基函数
        LagBaseIndex1D_x(1,:)=[];
        LagBaseIndex1D_x(end,:)=[];
        LagBaseIndex1D_x(1,4)=1;
        LagBaseIndex1D_x(end,4)=2;

        LagBaseIndex1D_y(1,:)=[];
        LagBaseIndex1D_y(end,:)=[];
        LagBaseIndex1D_y(1,4)=1;
        LagBaseIndex1D_y(end,4)=2;
        % 使用张量积生成两维信息表
        % 基函数信息
        Nx=size(LagBaseIndex1D_x,1);
        Ny=size(LagBaseIndex1D_y,1);
        index1=(1:Nx)';
        index2=(1:Ny)';
        [index1,index2]=meshgrid(index1,index2);
        index1=index1(:);
        index2=index2(:);
        LagBaseIndex=zeros(size(index1,1),8);
        LagBaseIndex(:,[1,2,5,7])=LagBaseIndex1D_x(index1,:);
        LagBaseIndex(:,[3,4,6,8])=LagBaseIndex1D_y(index2,:);
        % 网格信息
        Nx=size(T_x,1);
        Ny=size(T_y,1);
        index1=(1:Nx)';
        index2=(1:Ny)';
        [index1,index2]=meshgrid(index1,index2);
        index1=index1(:);
        index2=index2(:);
        LagMeshIndex=zeros(size(index1,1),13);
        LagMeshIndex(:,1:2)=T_x(index1,:);
        LagMeshIndex(:,3:4)=T_y(index2,:);
        for i=1:size(index1,1)
            indexBase=find(all(LagBaseIndex(:,1:4)==LagMeshIndex(i,1:4),2));
            for j=1:length(indexBase)
                LagMeshIndex(i,4+j)=indexBase(j);
            end
        end
end
end

