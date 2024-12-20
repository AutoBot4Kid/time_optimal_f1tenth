function [casadi_State] = solGen(optParam,svec,N,Si,So,Sid,Sod,Sidd,Sodd)
del_s = svec(2)-svec(1);
L = 0.34;

z1 =  optParam(1:5:end);
z2 =  optParam(2:5:end);
z3 =   optParam(3:5:end);
u1 =   optParam(4:5:end);
u2 =  optParam(5:5:end);
cost = 0;
for ind = 1:N

    xi = Si(ind,1);    yi = Si(ind,2);
    xo = So(ind,1);    yo = So(ind,2);

    xi_p = Sid(ind,1);    yi_p = Sid(ind,2);
    xo_p = Sod(ind,1);    yo_p = Sod(ind,2);

    xi_dp = Sidd(ind,1);    yi_dp = Sidd(ind,2);
    xo_dp = Sodd(ind,1);    yo_dp = Sodd(ind,2);

    lam = z1(ind);
    lmbda_p = z2(ind);
    lmbda_dp= u1(ind);
    b = z3(ind)^2;
    bp = 2*u2(ind);

    x_ds = lam.*xi_p + xi.*lmbda_p - xo.*lmbda_p - (lam -1).*xo_p;
    x_dds = lam.*xi_dp + xi.*lmbda_dp - xo.*lmbda_dp - (lam - 1).*xo_dp + 2.*lmbda_p.*xi_p - 2.*lmbda_p.*xo_p;

    y_ds= lam.*yi_p+ yi.*lmbda_p - yo.*lmbda_p - (lam - 1).*yo_p;
    y_dds =lam.*yi_dp + yi.*lmbda_dp - yo.*lmbda_dp - (lam - 1).*yo_dp + 2.*lmbda_p.*yi_p- 2.*lmbda_p.*yo_p;

    x_dT = lam*xi_p*sqrt(b) + xi*lmbda_p*sqrt(b) - xo*lmbda_p*sqrt(b) - xo_p*(lam - 1)*sqrt(b);
    x_ddT = lam*xi_p*bp/2 -xo_dp*(lam - 1)*sqrt(b)^2 + xi*lmbda_p*bp/2 - xo*lmbda_p*bp/2 + lam*xi_dp*sqrt(b)^2 + xi*lmbda_dp*sqrt(b)^2 + 2*lmbda_p*xi_p*sqrt(b)^2 - xo*lmbda_dp*sqrt(b)^2 - 2*lmbda_p*xo_p*sqrt(b)^2 - xo_p*(lam - 1)*bp/2;

    y_dT = lam*yi_p*sqrt(b) + yi*lmbda_p*sqrt(b) - yo*lmbda_p*sqrt(b) - yo_p*(lam - 1)*sqrt(b);
    y_ddT = lam*yi_p*bp/2 - yo_dp*(lam - 1)*sqrt(b)^2 + yi*lmbda_p*bp/2 - yo*lmbda_p*bp/2 + lam*yi_dp*sqrt(b)^2 + yi*lmbda_dp*sqrt(b)^2 + 2*lmbda_p*yi_p*sqrt(b)^2 - yo*lmbda_dp*sqrt(b)^2 - 2*lmbda_p*yo_p*sqrt(b)^2 - yo_p*(lam - 1)*bp/2;
    

    x =  lam*xi + (1-lam)*xo;
    y =  lam*yi + (1-lam)*yo;
    th = atan2(y_dT,x_dT);
    vT = sqrt(x_dT.^2 + y_dT.^2);
    wT = (x_dT.*y_ddT-y_dT.*x_ddT)./(vT).^2;
    phi = atan2(L*wT,vT);
    Along = (x_dT*x_ddT+ y_dT*y_ddT)/(vT);
    Alat = (x_dT*y_ddT- y_dT*x_ddT)/(vT);
    cost = cost + (1/z3(ind))*del_s;
    
    casadi_State(ind,:) = [svec(ind),cost,x,y,th,vT,phi,Along,Alat];
end


