function K = AssembleWaveMatrix1D_2(Coeff,interval,Table,supp,diff_test,diff_trail,type_test,type_trail,order)
%ASSEMBLEMATRIX1D ��װ�նȾ���
%   interval��������������
%   Table:��0Ԫ��������
%   supp:�����������j,k�Ķ�Ӧ��ϵ
%   diff_test: ���Ժ����ĵ�����
%   diff_trail:��̽�����ĵ�����
%   order:Gauss���ֵĽ�
%   Coeff�� ϵ������
K = zeros(2*size(Table,1)-size(supp,1),3);
pos=1;
for i = 1:size(Table,1)
    if Table(i,1)==Table(i,2)
        testFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,1),2),...
            supp(Table(i,1),3),diff_test,type_test,supp(Table(i,1),1));
        trailFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,2),2),...
            supp(Table(i,2),3),diff_trail,type_trail,supp(Table(i,2),1));
        f = @(x)testFun(x).*trailFun(x).*Coeff(x);
        K(pos,:) = [Table(i,1),Table(i,2),GaussIntForWavelet(f,interval,[Table(i,3),Table(i,4)],order,1,supp,[Table(i,1);Table(i,2)])]; % blockGauss
        pos=pos+1;
    else
        testFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,1),2),...
            supp(Table(i,1),3),diff_test,type_test,supp(Table(i,1),1));
        trailFun = @(x) @(x) WaveletBaseFun1D(x,interval,supp(Table(i,2),2),...
            supp(Table(i,2),3),diff_trail,type_trail,supp(Table(i,2),1));
        f=@(x)testFun(x).*trailFun(x).*Coeff(x);
        K(pos,:) = [Table(i,1),Table(i,2),GaussIntForWavelet(f,interval,[Table(i,3),Table(i,4)],order,1,supp,[Table(i,1);Table(i,2)])]; % blockGauss
        pos=pos+1;
        if diff_test~=diff_trail
            testFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,2),2),...
                supp(Table(i,2),3),diff_test,type_test,supp(Table(i,2),1));
            trailFun = @(x) WaveletBaseFun1D(x,interval,supp(Table(i,1),2),...
                supp(Table(i,1),3),diff_trail,type_trail,supp(Table(i,1),1));
            f=@(x)testFun(x).*trailFun(x).*Coeff(x);
            K(pos,:) = [Table(i,2),Table(i,1),GaussIntForWavelet(f,interval,[Table(i,3),Table(i,4)],order,1,supp,[Table(i,2);Table(i,1)])]; % blockGauss
        else
            K(pos,:) = [Table(i,2),Table(i,1),K(pos-1,3)];
        end
        pos=pos+1;
    end
end
K = sparse(K(:,1),K(:,2),K(:,3),size(supp,1),size(supp,1));
end

