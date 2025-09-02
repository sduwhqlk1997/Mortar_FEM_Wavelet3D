function [f1,f2,f3,g] = GenPiezoeleModel(ModelCoeff,u,v,w,x,y,z,phi)
load(ModelCoeff)
c_LN = cell2mat(materials(2));e_LN = cell2mat(materials(5));epcl0 = cell2mat(materials(6));epcl_LN = cell2mat(materials(7));
epcl_LN = epcl_LN*epcl0; 
% c_LN = c_LN*1e-10;
% epcl_LN = epcl_LN*1e10;

% syms u(x,y,z) v(x,y,z) w(x,y,z) phi(x,y,z) f1(x,y,z) f2(x,y,z) f3(x,y,z) g(x,y,z)
% 
% u(x,y,z) = x.^3+x.^2+1;
% v(x,y,z) = 1;
% w(x,y,z) = z.^3+z.^2+1;
% phi(x,y,z) = (x.^3+y.^3+z.^3)+1;
% var = [x;y;z];
sigma = @(i,j) C(c_LN,i,j,1,1)*diff(u,x,1)+C(c_LN,i,j,1,2)*(diff(v,x,1)+diff(u,y,1))+C(c_LN,i,j,1,3)*(diff(w,x,1)+diff(u,z,1))+...
    C(c_LN,i,j,2,2)*diff(v,y,1)+C(c_LN,i,j,2,3)*(diff(w,y,1)+diff(v,z,1))+...
    C(c_LN,i,j,3,3)*diff(w,z,1)+...
    E(e_LN,1,i,j)*diff(phi,x,1)+E(e_LN,2,i,j)*diff(phi,y,1)+E(e_LN,3,i,j)*diff(phi,z,1);
D = @(i) E(e_LN,i,1,1)*diff(u,x,1)+E(e_LN,i,1,2)*(diff(v,x,1)+diff(u,y,1))+E(e_LN,i,1,3)*(diff(w,x,1)+diff(u,z,1))+...
    E(e_LN,i,2,2)*diff(v,y,1)+E(e_LN,i,2,3)*(diff(w,y,1)+diff(v,z,1))+...
    E(e_LN,i,3,3)*diff(w,z,1)-...
    epcl_LN(i,1)*diff(phi,x,1)-epcl_LN(i,2)*diff(phi,y,1)-epcl_LN(i,3)*diff(phi,z,1);


f1 = matlabFunction(-(diff(sigma(1,1),x,1)+diff(sigma(1,2),y,1)+diff(sigma(1,3),z,1)));
f2 = matlabFunction(-(diff(sigma(2,1),x,1)+diff(sigma(2,2),y,1)+diff(sigma(2,3),z,1)));
f3 = matlabFunction(-(diff(sigma(3,1),x,1)+diff(sigma(3,2),y,1)+diff(sigma(3,3),z,1)));
g = matlabFunction(-(diff(D(1),x,1)+diff(D(2),y,1)+diff(D(3),z,1)));
end
% coef_phi = 1;
% u_exact = @(x,y,z) x.^3+x.^2+1;
% v_exact = @(x,y,z) zeros(size(x,1),1)+1;
% w_exact = @(x,y,z) z.^3+z.^2+1;
% phi_exact = @(x,y,z) coef_phi.*(x.^3+y.^3+z.^3)+1;
% ff1 =  -(6*(C(c_LN,1,1,1,1)+coef_phi.*E(e_LN,1,1,1)).*x+6*coef_phi.*E(e_LN,2,2,1).*y + 6*(C(c_LN,3,1,3,3)+coef_phi.*E(e_LN,3,3,1)).*z + 2*C(c_LN,1,1,1,1) + 2*C(c_LN,3,1,3,3));
% ff2 =  -(6*(C(c_LN,1,2,1,1)+coef_phi.*E(e_LN,1,1,2)).*x+6*coef_phi.*E(e_LN,2,2,2).*y + 6*(C(c_LN,3,2,3,3)+coef_phi.*E(e_LN,3,3,2)).*z + 2*C(c_LN,1,2,1,1) + 2*C(c_LN,3,2,3,3));
% ff3 =  -(6*(C(c_LN,1,3,1,1)+coef_phi.*E(e_LN,1,1,3)).*x+6*coef_phi.*E(e_LN,2,2,3).*y + 6*(C(c_LN,3,3,3,3)+coef_phi.*E(e_LN,3,3,3)).*z + 2*C(c_LN,1,3,1,1) + 2*C(c_LN,3,3,3,3));
% gg=  -(6*(E(e_LN,1,1,1)-coef_phi.*epcl_LN(1,1)).*x-6*coef_phi.*epcl_LN(2,2).*y+6*(E(e_LN,3,3,3)-coef_phi.*epcl_LN(3,3)).*z+2*E(e_LN,3,3,3)+2*E(e_LN,1,1,1));
% % f1-ff1
% % f2-ff2
% % f3-ff3
% g-gg