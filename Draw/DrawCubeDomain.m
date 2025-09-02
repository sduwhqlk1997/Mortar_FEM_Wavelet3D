function DrawCubeDomain(interval)
%DrawCubeDomain 绘制区域图
%   此处显示详细说明
[xx,yy,zz] = meshgrid(interval(1,:),interval(2,:),interval(3,:));
P = [xx(:),yy(:),zz(:)];
T = [1,2,3,4,5,6,7,8];
viewMesh(P,T)
end

