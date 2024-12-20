clear all ; clc ; close all

%%%%%% Track
track_name = 'YasMarina';
[parentDir, ~, ~] = fileparts(pwd);
track_dir = strcat(parentDir,'\map\',track_name,'\',track_name,'_raceline.csv');
track_load = readtable(track_dir);
sdata = track_load.x_S_m;
xdata = track_load.x_m(:);
ydata = track_load.y_m(:);


width = 2 ;N = 501; nf = 30; 
%flag =  1; % Cubic Spline Fit
%flag =  2 ; % Fourier Fit
flag = 3;
[svec,inner,outer,center] = getracetrack(xdata,ydata,width,N);
[Si,So,Sid,Sod,Sidd,Sodd] = getPathParam(flag,inner,outer,nf,svec);
   
%%%%%%% Optimal Path
paramGuess = [0.5*ones(N,1),1*ones(N,1),0.3*ones(N,1),100*ones(N,1),100*ones(N,1)];

for i  =  1 : 2 
optsol_collocation = casadicollocation(svec,N,Si,So,Sid,Sod,Sidd,Sodd,paramGuess);
z1 =  optsol_collocation(1:5:end);
z2 =  optsol_collocation(2:5:end);
z3 =   optsol_collocation(3:5:end);
u1 =   optsol_collocation(4:5:end);
u2 =  optsol_collocation(5:5:end);
paramGuess = [z1,z2,z3,u1,u2];
end
[state_casadi] = solGen(optsol_collocation,svec,N,Si,So,Sid,Sod,Sidd,Sodd);

%%%%%%% Check
odeOptions = odeset('RelTol',10^-8,'AbsTol',10^-8);
paramTab = [state_casadi(:,1),optsol_collocation(4:5:end),optsol_collocation(5:5:end)];
[ss,paramSol] = ode45(@paramODE,linspace(state_casadi(1,1),state_casadi(end,1),N*2),[optsol_collocation(1),optsol_collocation(2),optsol_collocation(3),0],odeOptions,paramTab);

uTab = [state_casadi(:,2),state_casadi(:,6),state_casadi(:,7)];
[tt,odeSol] = ode45(@biODE,linspace(state_casadi(1,2),state_casadi(end,2),N*2),[state_casadi(1,3),state_casadi(1,4),state_casadi(1,5)],odeOptions,uTab);

filename = strcat(track_name,'_N',num2str(N),'_F',num2str(flag));
save(filename)

%%%%% Limits 
L = 0.34;
vmin = 0.2;
vmax =  20;
amin = -2.0;
amax = 2.0;
mu = 0.3;
grav = 9.81;
phimin = -0.3;
phimax = 0.3;

%%%%%%% plot
data = load(filename);

Si = data.Si;
So = data.So;
state_casadi = data.state_casadi;
odeSol = data.odeSol;
paramSol = data.paramSol;
ss = data.ss;
tt = data.tt;
optsol_collocation = data.optsol_collocation;

figure(1)
subplot(3,2,1)
plot(Si(:,1),Si(:,2),'r','LineWidth',2);
hold on
plot(So(:,1),So(:,2),'b','LineWidth',2);
grid on 
hold on 
plot(state_casadi(:,3),state_casadi(:,4),'--','LineWidth',2);
hold on 
plot(odeSol(:,1),odeSol(:,2),':','LineWidth',2);
axis equal
title('XY')

subplot(3,2,2)
plot(state_casadi(:,2),state_casadi(:,5),'LineWidth',2);
hold on
plot(tt,wrapToPi(odeSol(:,3)),'--','LineWidth',2);
title('theta');
grid on

subplot(3,2,3)
plot(state_casadi(:,2),state_casadi(:,6),'LineWidth',2);
hold on 
plot(state_casadi(:,2),ones(1,length(state_casadi(:,2)))*vmax,'k','LineWidth',2);
hold on
plot(state_casadi(:,2),ones(1,length(state_casadi(:,2)))*vmin,'k','LineWidth',2);
grid on
title('v')

subplot(3,2,4)
plot(state_casadi(:,2),state_casadi(:,7),'LineWidth',2);
hold on
plot(state_casadi(:,2),ones(1,length(state_casadi(:,2)))*phimax,'k','LineWidth',2);
hold on
plot(state_casadi(:,2),ones(1,length(state_casadi(:,2)))*phimin,'k','LineWidth',2);
grid on
title('phi')

subplot(3,2,5)
plot(state_casadi(:,2),state_casadi(:,8),'LineWidth',2);
hold on
plot(state_casadi(:,2),ones(1,length(state_casadi(:,2)))*amax,'k','LineWidth',2);
hold on
plot(state_casadi(:,2),ones(1,length(state_casadi(:,2)))*amin,'k','LineWidth',2);
grid on
title('Acc')

subplot(3,2,6)
plot(mu*grav*cos(linspace(0,2*pi)),mu*grav*sin(linspace(0,2*pi)),'k','LineWidth',2);
hold on 
plot(state_casadi(:,8),state_casadi(:,9),'o','LineWidth',2);
grid on
title('GG')

figure(2)
subplot(3,1,1)
plot(state_casadi(:,1),optsol_collocation(1:5:end),'LineWidth',2);
hold on 
plot(ss,paramSol(:,1),'--','LineWidth',2);
grid on 
title('Lammbda')

subplot(3,1,2)
plot(state_casadi(:,1),optsol_collocation(2:5:end),'LineWidth',2);
hold on
plot(ss,paramSol(:,2),'--','LineWidth',2);
grid on 
title('Labmbda_d')

subplot(3,1,3)
plot(state_casadi(:,1),optsol_collocation(3:5:end),'LineWidth',2);
hold on
plot(ss,paramSol(:,3),'--','LineWidth',2);
grid on 
title('sdot')

figure(3)
plot(state_casadi(:,1),state_casadi(:,2),'LineWidth',2);
hold on
plot(ss,paramSol(:,4),'--','LineWidth',2);
hold on 
plot(ss,(paramSol(end,4)-paramSol(1,4)).*ss,':','LineWidth',2);
grid on 
title('Time Scaling')
