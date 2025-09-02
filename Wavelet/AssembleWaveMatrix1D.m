function K = AssembleWaveMatrix1D(Coeff,interval,Table,supp,IntegralBlock_ref,...
    diff_test,diff_trail,type_test,type_trail,order)
%AssembleWaveMatrix1D 组装刚度矩阵
%   interval：求解问题的区域
%   Table:非0元素索引表
%   supp:基函数编号与j,k的对应关系
%   IntegralBlock：一维小波基函数的积分分段，为矩阵，每一行表示一个基函数，非NaN数据表示节点（从左到右依次排列）
%   diff_test: 测试函数的导数阶
%   diff_trail:试探函数的导数阶
%   order:Gauss积分的阶
%   Coeff： 系数函数
K = zeros(2*size(Table,1)-size(supp,1),3);
pos=1;
err=[];
for i = 1:size(Table,1)
    IntegralBlock1_x=IntegralBlock_ref(Table(i,1),:);
    IntegralBlock2_x=IntegralBlock_ref(Table(i,2),:);
    IntegralBlock1_x(isnan(IntegralBlock1_x))=[];
    IntegralBlock2_x(isnan(IntegralBlock2_x))=[];
    if Table(i,1)==Table(i,2)
        testFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,1),2),...
            supp(Table(i,1),3),diff_test,type_test,supp(Table(i,1),1));
        trailFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,2),2),...
            supp(Table(i,2),3),diff_trail,type_trail,supp(Table(i,2),1));
        f = @(x)testFun(x).*trailFun(x).*Coeff(x);
        value=GaussIntForWavelet2(f,order,1,IntegralBlock1_x,IntegralBlock2_x);
        value2=integral(f,Table(i,3),Table(i,4));
        err=[err,[value;value2]];
        K(pos,:) = [Table(i,1),Table(i,2),value]; % blockGauss
        % K(pos,:) = [Table(i,1),Table(i,2),GaussIntForWavelet(f,interval,[Table(i,3),Table(i,4)],order,1,supp,[Table(i,1);Table(i,2)])]; % blockGauss
        pos=pos+1;
    else
        testFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,1),2),...
            supp(Table(i,1),3),diff_test,type_test,supp(Table(i,1),1));
        trailFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,2),2),...
            supp(Table(i,2),3),diff_trail,type_trail,supp(Table(i,2),1));
        f=@(x)testFun(x).*trailFun(x).*Coeff(x);
        value=GaussIntForWavelet2(f,order,1,IntegralBlock1_x,IntegralBlock2_x);
        K(pos,:) = [Table(i,1),Table(i,2),value]; % blockGauss
        % K(pos,:) = [Table(i,1),Table(i,2),GaussIntForWavelet(f,interval,[Table(i,3),Table(i,4)],order,1,supp,[Table(i,1);Table(i,2)])]; % blockGauss
        pos=pos+1;
        if diff_test~=diff_trail
            testFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,2),2),...
                supp(Table(i,2),3),diff_test,type_test,supp(Table(i,2),1));
            trailFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,1),2),...
                supp(Table(i,1),3),diff_trail,type_trail,supp(Table(i,1),1));
            f=@(x)testFun(x).*trailFun(x).*Coeff(x);
            value=GaussIntForWavelet2(f,order,1,IntegralBlock1_x,IntegralBlock2_x);
            K(pos,:) = [Table(i,2),Table(i,1),value]; % blockGauss
            % K(pos,:) = [Table(i,2),Table(i,1),GaussIntForWavelet(f,interval,[Table(i,3),Table(i,4)],order,1,supp,[Table(i,2);Table(i,1)])]; % blockGauss
        else
            K(pos,:) = [Table(i,2),Table(i,1),K(pos-1,3)];
        end
        pos=pos+1;
    end
end
K = sparse(K(:,1),K(:,2),K(:,3),size(supp,1),size(supp,1));
end

