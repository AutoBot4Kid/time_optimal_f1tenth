function optsol= casadicollocation(svec,N,Si,So,Sid,Sod,Sidd,Sodd,paramGuess)
import casadi.*

del_s = svec(2)-svec(1);
L = 0.34;
vmin = 0.2;
vmax =  20;
amin = -2.0;
amax = 2.0;
mu = 0.3;
grav = 9.81;
phimin = -0.3;
phimax = 0.3;

%%%%%%%%%%%%%%%%% Casadi

z1 =  SX.sym ('z1',N,1);
z2 =  SX.sym ('z2',N,1);
z3 =  SX.sym ('z3',N,1);
u1 =  SX.sym ('u1',N,1);
u2 = SX.sym ('u2',N,1);

w ={};
w0 = [];
lbw = [];
ubw = [];
g = {};
lbg  =[];
ubg =[];
cost = 0;
%%%%%%%%%%%%

X = [z1,z2,z3];
U = [u1,u2];

Xdot = rhs(svec,X,U);

Xll= X(1:end-1,:);
Xrr= X(2:end,:);

Ull = U(1:end-1,:);
Urr = U(2:end,:);

Xdotll = Xdot(1:end-1,:);
Xdotrr = Xdot(2:end,:);

Xc = 0.5*(Xll+Xrr)+ del_s/8 * (Xdotll-Xdotrr);
Uc = (Ull+Urr)*0.5;
XdotC = rhs(svec,Xc,Uc);
Ceq = (Xll-Xrr)+del_s/6*(Xdotll +4*XdotC +Xdotrr);
Ceq = [Ceq;[0,0,0]];

for ind = 1:N
    w = {w{:},z1(ind),z2(ind),z3(ind),u1(ind),u2(ind)};
    w0 = [w0;paramGuess(ind,1);paramGuess(ind,2);paramGuess(ind,3);paramGuess(ind,4);paramGuess(ind,5)];
    lbw =[lbw;0;-inf;0;-inf;-inf];
    ubw =[ubw;1;inf;inf;inf;inf];

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

    th = atan2(y_dT,x_dT);
    vT = sqrt(x_dT.^2 + y_dT.^2);
    wT = (x_dT.*y_ddT-y_dT.*x_ddT)./(vT).^2;
    phi = atan2(L*wT,vT);
    Along = (x_dT*x_ddT+ y_dT*y_ddT)/(vT);
    Alat = (x_dT*y_ddT- y_dT*x_ddT)/(vT);
    cost = cost + (1/z3(ind))*del_s;

    g = [g,{Ceq(ind,1)}];
    lbg = [lbg;0];
    ubg = [ubg;0];

    g = [g,{Ceq(ind,2)}];
    lbg = [lbg;0];
    ubg = [ubg;0];

    g = [g,{Ceq(ind,3)}];
    lbg = [lbg;0];
    ubg = [ubg;0];

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
disp(strcat('w0 : ',num2str(size(w0))))
disp(strcat('lbw : ',num2str(size(lbw))))
disp(strcat('ubw : ',num2str(size(ubw))))

disp(strcat('g : ',num2str(size(g'))))
disp(strcat('lbg : ',num2str(size(lbg))))
disp(strcat('ubg : ',num2str(size(ubg))))
disp(strcat('cost : ',num2str(size(cost))))

nlp= struct('x',vertcat(w{:}),'f',cost,'g',vertcat(g{:}));
options = struct;
options.ipopt.tol = 10^-10;
options.ipopt.max_iter = 10000;
options.ipopt.hessian_approximation = 'exact';
options.ipopt.jacobian_approximation= 'exact';
options.ipopt.gradient_approximation= 'exact';


S = nlpsol('S','ipopt',nlp,options);
sol = S('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
optsol = full(sol.x);

function [mdot] = rhs(s,m,u);
lambda = m(:,1);
lambdaD = m(:,2);
sdot = m(:,3);

gamma = u(:,1);
alpha = u(:,2);

lambda_Dot = lambdaD;
lambdaD_Dot = gamma ;
sdot_Dot =alpha;

mdot = [lambda_Dot, lambdaD_Dot, sdot_Dot ];

end


end