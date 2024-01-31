clear all
close all

% The error is set to zero at 2N distinct frequency points:
% wn = pi * [1 - 0.8*exp((-1.5*(k-1)*pi)/(2*N - 1))]

syms b0 b1 b2 b3 p1 p2

for i=1:6
    top1 = b0 + b1*cos(wn(i)) + b2*cos(2*wn(i)) + b3*cos(3*wn(i));
    top2 = b1*sin(wn(i)) + b2*sin(2*wn(i)) + b3*sin(3*wn(i));
%     top1 = 1 - (b1 + b2 + b3)*cos(wn(i)) + (b1*b2 + b2*b3 + b1*b3)*cos(2*wn(i))...
%            - b1*b2*b3*cos(3*wn(i));
% 
%     top2 = (b1 + b2 + b3)*sin(wn(i)) - (b1*b2 + b2*b3 + b1*b3)*sin(2*wn(i)) ...
%         + b1*b2*b3*sin(3*wn(i));

    bottom1 = 1 - (p1 + p2 + 1)*cos(wn(i)) + (p1 + p2 + p1*p2)*cos(2*wn(i)) ...
        - p1*p2*cos(3*wn(i));

    bottom2 = (p1 + p2 + 1)*sin(wn(i)) - (p1 + p2 + p1*p2)*sin(2*wn(i)) ...
        + p1*p2*sin(3*wn(i));

    top = top1^2 + top2^2;
    bottom = bottom1^2 + bottom2^2;

    uslov(i) = sqrt(top/bottom) - 1/wn(i) == 0;
%     uslov(i) = b0*sqrt(top/bottom) - 1/wn(i) == 0;
end

S = vpasolve(uslov, [b0, b1, b2, b3, p1, p2], [-100 100; -100 100; -100 100; -100 100; -1 1; -1 1])
 
b0d = double(S.b0)
b1d = double(S.b1)
b2d = double(S.b2)
b3d = double(S.b3)
p1d = double(S.p1)
p2d = double(S.p2)

% nule = [b1d b2d b3d]'
% b = b0d * poly(nule)

b = [b0d b1d b2d b3d]
nule=roots(b)

polovi=[1 p1d p2d]'
a = poly(polovi)
figure
zplane(nule,polovi)

[h,www]=freqz(b,a,1000);
figure
plot(www,www.*(abs(h)-1./www),'b','LineWidth',3)
figure
plot(www,(abs(h)-1./www),'r','LineWidth',3)
grid  

faza=unwrap(angle(h));
faza_id=-pi/2-0.5*www;
figure
plot(www,faza,'r','LineWidth',3)
title('Faza na pocetku')
grid

figure
plot(www(2:1000),180*(faza(2:1000)-faza_id(2:1000))/pi,'r','LineWidth',3)
title('Greska faze')
grid
