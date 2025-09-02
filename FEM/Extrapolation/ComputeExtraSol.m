function u_extra = ComputeExtraSol(u_CoarseGrid,s,m,Extra_type)
%COMPUTEEXTRASOL 计算外推解
%   u_CoarseGrid：粗网格解，每一列表示一个网格，从左到右依次加密即：(0,0,0)->(1,0,0)->...->(1,1,0)->...
%   m：外推参数
%   Extra_type：外推类型
switch s
    case 3
        switch m
            case 2 % (0,0,0) (1,0,0) (1,1,0) (2,0,0)
                switch Extra_type
                    case 1
                        w=[17/5;-4;16/9;64/45];
                    case 2
                        w=[83/24;-208/45;16/9;81/40];
                end

            case 3

            case 4

        end
    case 4
        switch m
            case 2 % (0,0,0,0) (1,0,0,0) (1,1,0,0) (2,0,0,0)  
                switch Extra_type
                    case 1
                        w=[349/45;-52/9;16/9;64/45];
                    case 2
                        w=[47/6;-32/5;16/9;81/40];
                end
        end
end
u_extra=u_CoarseGrid*w;
end

