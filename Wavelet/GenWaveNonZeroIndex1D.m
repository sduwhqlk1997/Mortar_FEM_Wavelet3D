function Table = GenWaveNonZeroIndex1D(supp)
%GenWaveNonZeroIndex1D ���ɽ�֧�������ǿյ�һά�������������
%   Table: ����Ԫ�������� 1~4��Ԫ������Ϊ����һ����������š��ڶ�����������š���֧��������˵㡢��֧�������Ҷ˵�
%   supp: ÿ��7��Ԫ�أ�����Ϊtype2(�߶Ȼ�(0)orС����(1))��j��k��֧����˵㡢֧���Ҷ˵㡢�����������������ͻ�������� 
[index1,index2]=meshgrid(supp(:,end));
index = [index1(:),index2(:)];
index(index(:,1)>index(:,2),:)=[];
index((supp(index(:,1),5)<=supp(index(:,2),4))|(supp(index(:,1),4)>=supp(index(:,2),5)),:)=[];
Table = [index,max(supp(index(:,1),4),supp(index(:,2),4)),min(supp(index(:,1),5),supp(index(:,2),5))];
end

