function y = WaveletBaseRef(BaseType,SpaceType,x,j,k,diff)
%WAVELETBASEREF : 一维参考小波基，定义在[0,1]上
%   BaseType：基函数类型
%   SpaceType：空间类型，0：尺度基函数，1：小波基函数
%   x：自变量
%   j：level
%   k：基函数编号
%   diff：导数阶

switch BaseType
    case "linear"
        disp("该基函数还未设置")
    case "quadratic"
        switch diff
            case 0
                phi = @(x) (x.^2./2).*(x>=0 & x<1) ...
                    + (-x.^2+3.*x-3/2).*(x>=1 & x<2) ...
                    + (x.^2/2-3.*x+9/2).*(x>=2 & x<=3);
                phi_b1 = @(x) (x.^2-2.*x+1).*(x>=0&x<=1);
                phi_b2 = @(x) (-3.*x.^2./2 + 2.*x).*(x>0&x<1)+(x.^2/2-2.*x+2).*(x>=1&x<=2);
                switch SpaceType
                    case 0 % 尺度基
                        if k==1
                            y = 2^(j/2).*phi_b1(2^(j).*x);
                        elseif k==2
                            y = 2^(j/2).*phi_b2(2^(j).*x);
                        elseif k>=3 && k<=2^j
                            y = 2^(j/2).*phi(2^j.*x-k+3);
                        elseif k==2^j+1
                            y = 2^(j/2).*phi_b2(2^(j).*(1-x));
                        elseif k==2^j+2
                            y = 2^(j/2).*phi_b1(2^(j).*(1-x));
                        else
                            error("k超出范围,该基函数不存在")
                        end
                    case 1 % 小波基
                        psi = @(x) -1/4.*phi(2.*x)+3/4.*phi(2.*x-1)-3/4.*phi(2.*x-2)+1/4.*phi(2.*x-3);
                        psi_bD = @(x) -phi_b2(2.*x)/4 + 47.*phi(2.*x)/120-13.*phi(2.*x-1)/40+phi(2.*x-2)/10;
                        if k==1
                            y = 2^(j/2).*psi_bD(2^j.*x);
                        elseif k>=2&&k<=2^j-1
                            y = 2^(j/2).*psi(2^j.*x-k+2);
                        elseif k==2^j
                            y = 2^(j/2).*psi_bD(2^j.*(1-x));
                        else
                            error("k超出范围,该基函数不存在")
                        end
                    otherwise
                        error('SpaceType的值只能为0或者1')
                end
            case 1
                dphi = @(x) x.*(x>=0 & x<1) ...
                    + (-2.*x+3).*(x>=1 & x<2) ...
                    + (x-3).*(x>=2 & x<=3);
                dphi_b1 = @(x) (2.*x-2).*(x>=0 & x<=1);
                dphi_b2 = @(x) (-3.*x+2).*(x>=0 & x<1)+(x-2).*(x>=1 & x<=2);
                switch SpaceType
                    case 0 % 尺度基
                        if k==1
                            y = 2^(3.*j/2).*dphi_b1(2^(j).*x);
                        elseif k==2
                            y = 2^(3.*j/2).*dphi_b2(2^(j).*x);
                        elseif k>=3 && k<=2^j
                            y = 2^(3.*j/2).*dphi(2^j.*x-k+3);
                        elseif k==2^j+1
                            y = -2^(3.*j/2).*dphi_b2(2^(j).*(1-x));
                        elseif k==2^j+2
                            y = -2^(3.*j/2).*dphi_b1(2^(j).*(1-x));
                        else
                            error("k超出范围,该基函数不存在")
                        end
                    case 1 % 小波基
                        dpsi = @(x) -1/2.*dphi(2.*x)+3/2.*dphi(2.*x-1)-3/2.*dphi(2.*x-2)+1/2.*dphi(2.*x-3);
                        dpsi_bD = @(x) -dphi_b2(2.*x)/2 + 47.*dphi(2.*x)/60-13.*dphi(2.*x-1)/20+dphi(2.*x-2)/5;
                        if k==1
                            y = 2^(3.*j/2).*dpsi_bD(2^j.*x);
                        elseif k>=2&&k<=2^j-1
                            y = 2^(3.*j/2).*dpsi(2^j.*x-k+2);
                        elseif k==2^j
                            y = -2^(3.*j/2).*dpsi_bD(2^j.*(1-x));
                        else
                            error("k超出范围,该基函数不存在")
                        end
                    otherwise
                        error('SpaceType的值只能为0或者1')
                end
        end
    otherwise
        error("基函数不存在")
end
end

