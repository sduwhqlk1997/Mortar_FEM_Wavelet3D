function supp = GenSuppofWavelet1D(interval,type1,type2,j,k)
%GENSUPPOFWAVELET1D 生成一维尺度or小波基函数的支集
%   此处显示详细说明
if type1=="linear"
    if type2==0
        supp_phi = [0,2];
        supp_phi_b = [0,1];
        if k==1
            supp=supp_phi_b.*2^(-j);
        elseif k==2^j+1
            supp=[1-2^(-j),1];
        else
            supp=(supp_phi+k-2).*2^(-j);
        end
    else
        supp_psi = [0,3];
        supp_psi_b = [0,2];
        if k==1
            supp=supp_psi_b.*2^(-j);
        elseif k==2^j
            supp=[1-2^(1-j),1];
        else
            supp=(supp_psi+k-2).*2^(-j);
        end
    end
elseif type1=="quadratic"
    supp_phi = [0,3];
    supp_phi_b1 = [0,1];
    supp_phi_b2 = [0,2];
    supp_psi = [0,3];
    supp_psi_bD = [0,2.5];
    if type2==0
        if k==1
            supp=supp_phi_b1.*2^(-j);
        elseif k==2
            supp=supp_phi_b2.*2^(-j);
        elseif k==2^j+1
            supp=[1-2^(-j+1),1];
        elseif k==2^j+2
            supp=[1-2^(-j),1];
        else
            supp=(supp_phi+k-3).*2^(-j);
        end
    elseif type2==1
        if k==1
            supp=supp_psi_bD.*2^(-j);
        elseif k==2^j
            supp=[1-5/2^(j+1),1];
        else
            supp=(supp_psi+k-2).*2^(-j);
        end
    end
else
    error("该基函数不存在")
end
supp = supp.*(interval(2)-interval(1))+interval(1);
end

