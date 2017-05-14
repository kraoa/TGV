function [f,g] = mpnrage_2fa(TR,T1,fa,b,TI,TD,num_ti,k)
% clear
% N = 375;
% M = 375;
% TR = 4.9;
% T1 = 1200;
% TI = [0 0];
% TD = [250 250];
% b = [180 180];
% fa = [4 8];
% k = 1;

N = num_ti(1);
M = num_ti(2);

n = [0:N-1]';
m = [0:M-1]';

fa1 = fa(1);
fa2 = fa(2);

TD1 = TD(1);
TD2 = TD(2);

b1 = b(1);
b2 = b(2);


E1 = exp(-TR/T1);


EI = exp(-TI/T1);
EI1 = EI(1);
EI2 = EI(2);

ED = exp(-TD/T1);
ED1 = ED(1);
ED2 = ED(2);

c1 = cosd(k.*fa2).*cosd(b1).*ED2.*EI1;
c2 = cosd(k.*fa1).*cosd(b2).*ED1.*EI2;

a1 = cosd(k.*fa1).*E1;
a2 = cosd(k.*fa2).*E1;

B1 = (1-ED2).*EI1.*cosd(b1) + (1-EI1);
B2 = (1-ED1).*EI2.*cosd(b2) + (1-EI2);


A1 = (1-E1)./(1-a1);%sind(k.*fa1).*
A2 = (1-E1)./(1-a2);%sind(k.*fa2).*

f = A1.*(1-(a1.^n)) + A2.*c1.*(1-a2.^(M-1)).*(a1.^n) + B2.*c1.*(a2.^(M-1)).*(a1.^n) + B1.*(a1.^n);
f_N1 = f(end)./(1-c1.*c2.*(a2.^(M-1)).*(a1.^(N-1)));
f = f + f_N1.*c1.*c2.*(a1.^n).*(a2.^(M-1));


g = A2.*(1-a2.^m) + A1.*c2.*(1-a1.^(N-1)).*(a2.^m) + B1.*c2.*(a1.^(M-1)).*(a2.^m) + B2.*(a2.^m);
g_M1 = g(end)./(1-c2.*c1.*(a1.^(N-1)).*(a2.^(M-1)));
g = g + g_M1.*c2.*c1.*(a2.^m).*(a1.^(N-1));


[s1 s2] = mpnrage_dir_signal(TR,T1,fa.*k,b,TI,TD,[N M],[1 1]);

f = f.*sind(k.*fa1);
g = g.*sind(k.*fa2);

s1 = s1.*sind(k.*fa1);
s2 = s2.*sind(k.*fa2);

% hold on
% plot(f)
% plot(g)
% 
% 
% plot(s1,'o')
% plot(s2,'o')