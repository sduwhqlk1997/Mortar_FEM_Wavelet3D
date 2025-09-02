function interval = GenIntegralInterval(interval1,interval2)
%GENINTEGRALINTERVAL 构造两个被积函数相乘的积分区域
%   此处显示详细说明
%Scale=interval1(:,2)-interval1(:,1);
if any(interval1(:,1)>=interval2(:,2))||any(interval2(:,1)>=interval1(:,2))
    interval=[];
else
    interval=[max(interval1(:,1),interval2(:,1)),min(interval1(:,2),interval2(:,2))];
end
end

