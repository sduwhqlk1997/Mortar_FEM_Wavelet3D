function [X,Y,Z] = GenWaveGaussBlock(interval_Sol,interval,dim,BaseIndex_x,DofX_num,BaseIndex_y,DofY_num,BaseIndex_z,DofZ_num)
%GenWaveGaussBlock ���ɷ�ƬGauss���ֵķ�Ƭ��Ϣ
%   interval_Sol:����������
%   interval:�������䣬Ϊdim*2�ľ���ÿһ�д���һ��ά��
%   dim:����ά��
%   BaseIndex_x,BaseIndex_y,BaseIndex_z:�������귽��Ļ�����������ÿһ�д���һ�������������еڶ��б�ʾ�û�������level
%   DofX_num,DofY_num,DofZ_num:����������������ɶȱ��Ϊ������
%   X,Y,Z:��������Ļ��ַ�Ƭ
if dim==1
    jmax = max(BaseIndex_x(DofX_num,2)); % ��ȡ��ά�ȵ����level
    d = 2^(-jmax-1)*(interval_Sol(1,2)-interval_Sol(1,1));
    X = (interval(1,1):d:interval(1,2))';
    if length(X)<=1
        X = interval;
    end
elseif dim==2
    jmax_x = max(BaseIndex_x(DofX_num,2));
    jmax_y = max(BaseIndex_y(DofY_num,2));
    
    d = 2^(-jmax_x-1)*(interval_Sol(1,2)-interval_Sol(1,1));
    X = (interval(1,1):d:interval(1,2))';
    if length(X)<=1
        X = interval(1,:);
    end
    
    d = 2^(-jmax_y-1)*(interval_Sol(2,2)-interval_Sol(2,1));
    Y = (interval(2,1):d:interval(2,2))';
    if length(Y)<=1
        Y = interval(2,:);
    end
elseif dim==3 
    jmax_x = max(BaseIndex_x(DofX_num,2));
    jmax_y = max(BaseIndex_y(DofY_num,2));
    jmax_z = max(BaseIndex_z(DofZ_num,2));
    
    d = 2^(-jmax_x-1)*(interval_Sol(1,2)-interval_Sol(1,1));
    X = (interval(1,1):d:interval(1,2))';
    if length(X)<=1
        % 1
        d = 2^(-jmax_x-2)*(interval_Sol(1,2)-interval_Sol(1,1));
        X = (interval(1,1):d:interval(1,2))';
    end
    
    d = 2^(-jmax_y-1)*(interval_Sol(2,2)-interval_Sol(2,1));
    Y = (interval(2,1):d:interval(2,2))';
    if length(Y)<=1
        % 1
        d = 2^(-jmax_y-2)*(interval_Sol(2,2)-interval_Sol(2,1));
        Y = (interval(2,1):d:interval(2,2))';
    end
    
    d = 2^(-jmax_z-1)*(interval_Sol(3,2)-interval_Sol(3,1));
    Z = (interval(3,1):d:interval(3,2))';
    if length(Z)<=1
        % 1
        d = 2^(-jmax_z-2)*(interval_Sol(3,2)-interval_Sol(3,1));
        Z = (interval(3,1):d:interval(3,2))';
    end
else
    error("ά�Ȳ��ܴ�����ά������С��һά")
end
end

