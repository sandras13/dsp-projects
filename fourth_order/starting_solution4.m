clear all
close all

% teta0 = 0.6585
% teta1 = 1.4144
% teta2 = 2.0584
% teta3 = 2.4890
% teta4 = 2.8370
% teta5 = 3.036
% teta6 = 3.1078
% teta7 = 3.1356
% 
% teta0 = 0.6
% teta1 = 1.41
% teta2 = 2.04
% teta3 = 2.4
% teta4 = 2.8
% teta5 = 3.0
% teta6 = 3.1
% teta7 = 3.13

% 180 greska
% teta0 = 0.6
% teta1 = 1.41
% teta2 = 2.041
% teta3 = 2.4
% teta4 = 2.8
% teta5 = 3.0
% teta6 = 3.0998
% teta7 = 3.13

% teta0 = 0.6
% teta1 = 1.41
% teta2 = 2.041
% teta3 = 2.42
% teta4 = 2.8
% teta5 = 3.0
% teta6 = 3.0589
% teta7 = 3.13

% % 90 stepeni greska
% teta0 = 0.6
% teta1 = 1.41
% teta2 = 2.041
% teta3 = 2.42
% teta4 = 2.8
% teta5 = 2.95
% teta6 = 3.0589
% teta7 = 3.13
% 
teta0 = 0.6
teta1 = 1.2
teta2 = 1.75
teta3 = 2.408
teta4 = 2.85
teta5 = 3.048
teta6 = 3.10
teta7 = 3.13

% teta0 = 0.6
% teta1 = 1.2
% teta2 = 1.75
% teta3 = 2.408
% teta4 = 2.845
% teta5 = 3.048
% teta6 = 3.10
% teta7 = 3.13
% 
% 
% teta0 = 0.6
% teta1 = 1.2
% teta2 = 1.757
% teta3 = 2.40
% teta4 = 2.8
% teta5 = 2.95
% teta6 = 3.05
% teta7 = 3.13


wn=[teta0; teta1; teta2; teta3; teta4; teta5; teta6; teta7]
syms b0 b1 b2 b3 b4 p1 p2 p3

for i=1:8
    top1 = b0 + b1*cos(wn(i)) + b2*cos(2*wn(i)) + b3*cos(3*wn(i)) + b4*cos(4*wn(i));
    top2 = b1*sin(wn(i)) + b2*sin(2*wn(i)) + b3*sin(3*wn(i)) + b4*sin(4*wn(i));

%     top1 = 1 - (b1 + b2 + b3 + b4)*cos(wn(i)) + (b1*b2 + b2*b3 + b1*b3 + b1*b4 + b2*b4 + b3*b4)*cos(2*wn(i))...
%         - (b1*b2*b3 + b1*b2*b4 + b1*b3*b4 + b2*b3*b4)*cos(3*wn(i)) + b1*b2*b3*b4*cos(4*wn(i));
% 
%     top2 = (b1 + b2 + b3 + b4)*sin(wn(i)) - (b1*b2 + b2*b3 + b1*b3 + b1*b4 + b2*b4 + b3*b4)*sin(2*wn(i)) ...
%         + (b1*b2*b3 + b1*b2*b4 + b1*b3*b4 + b2*b3*b4)*sin(3*wn(i)) - b1*b2*b3*b4*sin(4*wn(i));

    bottom1 = 1 - (p1 + p2 + p3 + 1)*cos(wn(i)) + (p1*p2 + p2*p3 + p3*p1 + p1 + p2 + p3)*cos(2*wn(i)) ...
        - (p1*p2 + p1*p2*p3 + p1*p3 + p2*p3)*cos(3*wn(i)) + p1*p2*p3*cos(4*wn(i));

    bottom2 = (p1 + p2 + p3 + 1)*sin(wn(i)) - (p1*p2 + p2*p3 + p3*p1 +  p1 + p2 + p3)*sin(2*wn(i)) ...
        + (p1*p2 + p1*p2*p3 + p1*p3 + p2*p3)*sin(3*wn(i)) - p1*p2*p3*sin(4*wn(i));

    top = top1^2 + top2^2;
    bottom = bottom1^2 + bottom2^2;

%     uslov(i) = b0*sqrt(top/bottom) - 1/wn(i) == 0;
    uslov(i) = sqrt(top/bottom) - 1/wn(i) == 0;
end

S = vpasolve(uslov, [b0, b1, b2, b3, b4, p1, p2, p3], ...
    [-100 100; -100 100; -100 100; -100 100; -100 100; -1 1; -1 1; -1 1])

 
b0d = double(S.b0)
b1d = double(S.b1)
b2d = double(S.b2)
b3d = double(S.b3)
b4d = double(S.b4)
p1d = double(S.p1)
p2d = double(S.p2)
p3d = double(S.p3)

b = [b0d b1d b2d b3d b4d]
nule=roots(b)

% nule = [b1d b2d b3d b4d]'
% b = b0d * poly(nule)

polovi=[1 p1d p2d p3d]'
a = poly(polovi)
figure
zplane(nule,polovi)

[h,www]=freqz(b,a,1000);
figure
plot(www,www.*(abs(h)-1./www),'b','LineWidth',3)
%figure
%plot(www,(abs(h)-1./www),'r','LineWidth',3)
grid  

faza=unwrap(angle(h));
faza_id=-pi/2+0.5*www;
figure
plot(www/pi,faza,'r','LineWidth',3)
title('Faza na pocetku')
grid

figure
plot(www(2:1000)/pi,180*(faza(2:1000)-faza_id(2:1000))/pi,'r','LineWidth',3)
title('Greska faze')
grid