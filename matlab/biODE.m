
function [mdot] = biODE(t,m,tab);
x = m(1);
y = m(2);
th = m(3);

v = interp1(tab(:,1),tab(:,2),t);
phi = interp1(tab(:,1),tab(:,3),t);
L = 0.34;

xdot = v*cos(th);
ydot = v*sin(th);
thdot = v/L*tan(phi);
mdot = [xdot;ydot;thdot];
end
