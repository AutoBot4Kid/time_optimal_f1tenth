function optsol= casadishooting(svec,N,Si,So,Sid,Sod,Sidd,Sodd)
import casadi.*

del_s = svec(2)-svec(1);
L = 0.34;
vmin = 0.2;
vmax =  20;
amin = -3.0;
amax = 3.0;
mu = 0.3;
grav = 9.81;
phimin = -0.3;
phimax = 0.3;

%%%%%%%%%%%%%%%%% Casadi
Z =  SX.sym ('z',3);% z(1) = lambda, z(2) = lambda_p, z(3) = sdot
U =  SX.sym('u',2); % u(1) = lambda_dp, u(2) = a
f = [Z(2);U(1);U(2)];

dae = struct('x',Z,'p',U ,'ode',f,'alg',[],'quad',[]);
optionODE = struct('tf',svec(end)/N,'abstol',1e-8);
F = integrator('F','idas',dae,optionODE);

w ={};
w0 = [];
lbw = [];
ubw = [];
g = {};
lbg  =[];
ubg =[];
cost = 0;


Zk = MX.sym(['z' num2str(1)],3);
w= [w,{Zk}];
w0 = [w0;0.5;0.5;0.5];
lbw = [lbw;0;-inf;0];
ubw = [ubw;1;inf;inf];


for ind = 1:length(Si)-1
    Uk = MX.sym(['u' num2str(ind)],2);
    w = {w{:},Uk};
    w0 = [w0;100;100];
    lbw = [lbw;-inf;-inf];
    ubw = [ubw;inf;inf];

    Fk = F('x0',Zk,'p',Uk);
    Zk_end = Fk.xf;

    Zk = MX.sym(['z' num2str(ind+1)],3);
    w= [w, {Zk}];
    w0 = [w0;0.5;0;1];
    lbw = [lbw;0;-inf;0];
    ubw = [ubw;1;inf;inf];

    xi = Si(ind+1,1);    yi = Si(ind+1,2);
    xo = So(ind+1,1);    yo = So(ind+1,2);

    xi_p = Sid(ind+1,1);    yi_p = Sid(ind+1,2);
    xo_p = Sod(ind+1,1);    yo_p = Sod(ind+1,2);

    xi_dp = Sidd(ind+1,1);    yi_dp = Sidd(ind+1,2);
    xo_dp = Sodd(ind+1,1);    yo_dp = Sodd(ind+1,2);

    lam = Zk(1);
    lmbda_p = Zk(2);
    lmbda_dp= Uk(1);
    b = Zk(3)^2;
    bp = 2*Uk(2);

    x_ds = lam.*xi_p + xi.*lmbda_p - xo.*lmbda_p - (lam -1).*xo_p;
    x_dds = lam.*xi_dp + xi.*lmbda_dp - xo.*lmbda_dp - (lam - 1).*xo_dp + 2.*lmbda_p.*xi_p - 2.*lmbda_p.*xo_p;

    y_ds= lam.*yi_p+ yi.*lmbda_p - yo.*lmbda_p - (lam - 1).*yo_p;
    y_dds =lam.*yi_dp + yi.*lmbda_dp - yo.*lmbda_dp - (lam - 1).*yo_dp + 2.*lmbda_p.*yi_p- 2.*lmbda_p.*yo_p;

    x_dT = lam*xi_p*sqrt(b) + xi*lmbda_p*sqrt(b) - xo*lmbda_p*sqrt(b) - xo_p*(lam - 1)*sqrt(b);
    x_ddT = lam*xi_p*bp/2 -xo_dp*(lam - 1)*sqrt(b)^2 + xi*lmbda_p*bp/2 - xo*lmbda_p*bp/2 + lam*xi_dp*sqrt(b)^2 + xi*lmbda_dp*sqrt(b)^2 + 2*lmbda_p*xi_p*sqrt(b)^2 - xo*lmbda_dp*sqrt(b)^2 - 2*lmbda_p*xo_p*sqrt(b)^2 - xo_p*(lam - 1)*bp/2;

    y_dT = lam*yi_p*sqrt(b) + yi*lmbda_p*sqrt(b) - yo*lmbda_p*sqrt(b) - yo_p*(lam - 1)*sqrt(b);
    y_ddT = lam*yi_p*bp/2 - yo_dp*(lam - 1)*sqrt(b)^2 + yi*lmbda_p*bp/2 - yo*lmbda_p*bp/2 + lam*yi_dp*sqrt(b)^2 + yi*lmbda_dp*sqrt(b)^2 + 2*lmbda_p*yi_p*sqrt(b)^2 - yo*lmbda_dp*sqrt(b)^2 - 2*lmbda_p*yo_p*sqrt(b)^2 - yo_p*(lam - 1)*bp/2;

    th = atan2(y_dT,x_dT);
    vT = sqrt(x_dT.^2 + y_dT.^2);
    wT = (x_dT.*y_ddT-y_dT.*x_ddT)./(vT).^2;
    phi = atan2(L*wT,vT);
    Along = (x_dT*x_ddT+ y_dT*y_ddT)/(vT);
    Alat = (x_dT*y_ddT- y_dT*x_ddT)/(vT);
    cost = cost + (1/Zk(3))*del_s;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    shooting = Zk_end-Zk;
    g= [g,{shooting}];
    lbg = [lbg;0;0;0];
    ubg = [ubg;0;0;0];

    g =[g,{vT}];
    lbg =[lbg;vmin];
    ubg =[ubg;vmax];

    g =[g,{phi}];
    lbg =[lbg;phimin];
    ubg =[ubg;phimax];

    g = [g,{Along}];
    lbg = [lbg;amin];
    ubg = [ubg;amax];

    g = [g,{Along^2+Alat^2}];
    lbg = [lbg;0];
    ubg = [ubg;(mu*grav)^2];
end

disp ('SIZE CHECKING ')
disp(strcat('w : ',num2str(size(w0'))))
disp(strcat('lbw : ',num2str(size(lbw))))
disp(strcat('ubw : ',num2str(size(ubw))))

disp(strcat('g : ',num2str(size(g'))))
disp(strcat('lbg : ',num2str(size(lbg))))
disp(strcat('ubg : ',num2str(size(ubg))))
disp(strcat('cost : ',num2str(size(cost))))

nlp= struct('x',vertcat(w{:}),'f',cost,'g',vertcat(g{:}));
options = struct;
options.ipopt.tol = 10^-10;
options.ipopt.max_iter = 50000;
options.ipopt.hessian_approximation = 'exact';
options.ipopt.jacobian_approximation= 'exact';
options.ipopt.gradient_approximation= 'exact';

S = nlpsol('S','ipopt',nlp,options);
sol = S('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
optsol = full(sol.x);
end