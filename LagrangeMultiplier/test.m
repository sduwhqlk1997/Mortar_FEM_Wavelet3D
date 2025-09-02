%% 测试LagBaseIndex = GenLagBaseIndex(interval,baseType,x_step,y_step,Nunknows)
baseType="linear";
interval=[0,1;0,3];
x_step=linspace(interval(1,1),interval(1,2),9)';
y_step=linspace(interval(2,1),interval(2,2),9)';
[LagBaseIndex,LagMeshIndex] = GenLagBaseIndex(baseType,x_step,y_step);

% yy=(interval(2,1):0.01:interval(2,2))';
% xx=ones(length(yy),1);

figure(1)
ylim([0,2])
index=LagBaseIndex(:,1)==x_step(2)&LagBaseIndex(:,5)==1;
LagBaseIndex=LagBaseIndex(index,:);
for i=1:size(LagBaseIndex,1)
    yy=(LagBaseIndex(i,3):0.001:LagBaseIndex(i,4))';
    xx=x_step(2)*ones(length(yy),1);
    ff=LagBaseFun(xx,yy,LagBaseIndex(i,7:8),baseType,...
        [LagBaseIndex(i,1:2);LagBaseIndex(i,3:4)],...
        LagBaseIndex(i,5:6),[0;0]);
    plot(yy,ff)
    hold on
end
ylim([0,2])
