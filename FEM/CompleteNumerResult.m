function [Dof_index,Sol] = CompleteNumerResult(Dof_index,Sol,Bound_Info)
%COMPLETENUMERRESULT 补全边界上的数值解
%   Bound_Info边界上数值解的信息
%       1：自由度类型
%       2：坐标矩阵
%       3：2中各点上的数值解
%       4: 自由度编号

for i=1:size(Bound_Info,1)
    Dof_index = [Dof_index;cell2mat(Bound_Info(i,1)),cell2mat(Bound_Info(i,2)),...
        cell2mat(Bound_Info(i,4))];
    Sol = [Sol;cell2mat(Bound_Info(i,3))];
end
end

