function eijk = E(coef,i,j,k)
%E ��ѹ���洢���ѹ����������ȡ����
%   coef��ѹ���洢���ѹ������
%   ijk:ѹ���洢ǰ������Ӧ���±�
if j==k
    m=j;
elseif j+k==5
    m=4;
elseif j+k==4
    m=5;
elseif j+k==3
    m=6;
end
eijk = coef(i,m);
end

