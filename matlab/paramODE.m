
function [mdot] = paramODE(t,m,tab);
z1 = m(1);
z2 = m(2);
z3 = m(3);
j = m(4);

u1 = interp1(tab(:,1),tab(:,2),t);
u2 = interp1(tab(:,1),tab(:,3),t);

z1dot = z2;
z2dot = u1;
z3dot = u2;
jdot = 1/(z3);

mdot = [z1dot;z2dot;z3dot;jdot];
end

