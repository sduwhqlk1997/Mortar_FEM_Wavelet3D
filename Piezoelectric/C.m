function Cijkl = C(coef,i,j,k,l)
%C ��ѹ���洢�ĵ��������а����±���ȡ����
%   coef:ѹ���洢�ĵ�������
%   i,j,k,lԪ��ѹ��ǰ���±�
if i==j
    m=i;
elseif i+j==5
    m=4;
elseif i+j==4
    m=5;
elseif i+j==3
    m=6;
end
if k==l
    n=k;
elseif k+l==5
    n=4;
elseif k+l==4
    n=5;
elseif k+l==3
    n=6;
end
Cijkl=coef(m,n);
end

