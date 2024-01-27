clear all
close all


% teta0=0.005
% teta1=0.347
% teta2=0.706
% teta3=0.965

teta0 = 0.059
teta1 = 1.799
teta2 = 2.623
teta3 = 3.097

teta0 = 0.6605
teta1 = 2.2605
teta2 = 2.7305
teta3 = 3.105

% teta0 = 0.4
% teta1 = 2.25
% teta2 = 2.7305
% teta3 = 3.105

% teta0 = 0.702
% teta1 = 2.25
% teta2 = 2.75
% teta3 = 3.097


wn=[teta0; teta1; teta2; teta3]
syms b0 b1 b2 p
j=sqrt(-1)

% S=solve(abs(b0+b1*exp(-j*wn(1))+b2*exp(-j*2*wn(1)))/abs(1-(p+1)*exp(-j*wn(1))+p*exp(-j*2*wn(1)))- 1/wn(1),...
%     abs(b0+b1*exp(-j*wn(2))+b2*exp(-j*2*wn(2)))/abs(1-(p+1)*exp(-j*wn(2))+p*exp(-j*2*wn(2)))- 1/wn(2),...
%     abs(b0+b1*exp(-j*wn(3))+b2*exp(-j*2*wn(3)))/abs(1-(p+1)*exp(-j*wn(3))+p*exp(-j*2*wn(3)))- 1/wn(3),...
%     abs(b0+b1*exp(-j*wn(4))+b2*exp(-j*2*wn(4)))/abs(1-(p+1)*exp(-j*wn(4))+p*exp(-j*2*wn(4)))- 1/wn(4))

% S=solve(((b0+b1*cos(wn(1))+b2*cos(2*wn(1)))^2+(b1*sin(wn(1))+b2*sin(2*wn(1)))^2)...
% /(((1-(p+1)*cos(wn(1))+p*cos(2*wn(1)))^2+((p+1)*sin(wn(1))-p*sin(2*wn(1)))^2))- 1/wn(1)^2,...
% ((b0+b1*cos(wn(2))+b2*cos(2*wn(2)))^2+(b1*sin(wn(2))+b2*sin(2*wn(2)))^2)...
% /(((1-(p+1)*cos(wn(2))+p*cos(2*wn(2)))^2+((p+1)*sin(wn(2))-p*sin(2*wn(2)))^2))- 1/wn(2)^2,...   
% ((b0+b1*cos(wn(3))+b2*cos(2*wn(3)))^2+(b1*sin(wn(3))+b2*sin(2*wn(3)))^2)...
% /(((1-(p+1)*cos(wn(3))+p*cos(2*wn(3)))^2+((p+1)*sin(wn(3))-p*sin(2*wn(3)))^2))- 1/wn(3)^2,... 
% ((b0+b1*cos(wn(4))+b2*cos(2*wn(4)))^2+(b1*sin(wn(4))+b2*sin(2*w n(4)))^2)...
% /(((1-(p+1)*cos(wn(4))+p*cos(2*wn(4)))^2+((p+1)*sin(wn(4))-p*sin(2*wn(4)))^2))- 1/wn(4)^2)




% S=solve(sqrt(((b0+b1*cos(wn(1))+b2*cos(2*wn(1)))^2+(b1*sin(wn(1))+b2*sin(2*wn(1)))^2)...
% /(((1-(p+1)*cos(wn(1))+p*cos(2*wn(1)))^2+((p+1)*sin(wn(1))-p*sin(2*wn(1)))^2)))- 1/wn(1),...
% sqrt(((b0+b1*cos(wn(2))+b2*cos(2*wn(2)))^2+(b1*sin(wn(2))+b2*sin(2*wn(2)))^2)...
% /(((1-(p+1)*cos(wn(2))+p*cos(2*wn(2)))^2+((p+1)*sin(wn(2))-p*sin(2*wn(2)))^2)))- 1/wn(2),...   
% sqrt(((b0+b1*cos(wn(3))+b2*cos(2*wn(3)))^2+(b1*sin(wn(3))+b2*sin(2*wn(3)))^2)...
% /(((1-(p+1)*cos(wn(3))+p*cos(2*wn(3)))^2+((p+1)*sin(wn(3))-p*sin(2*wn(3)))^2)))- 1/wn(3),... 
% sqrt(((b0+b1*cos(wn(4))+b2*cos(2*wn(4)))^2+(b1*sin(wn(4))+b2*sin(2*wn(4)))^2)...
% /(((1-(p+1)*cos(wn(4))+p*cos(2*wn(4)))^2+((p+1)*sin(wn(4))-p*sin(2*wn(4)))^2)))- 1/wn(4))


S=solve((sqrt(((b0+b1*cos(wn(1))+b2*cos(2*wn(1)))^2+(b1*sin(wn(1))+b2*sin(2*wn(1)))^2)...
/(((1-(p+1)*cos(wn(1))+p*cos(2*wn(1)))^2+((p+1)*sin(wn(1))-p*sin(2*wn(1)))^2)))- 1/wn(1))/(1/wn(1)),...
(sqrt(((b0+b1*cos(wn(2))+b2*cos(2*wn(2)))^2+(b1*sin(wn(2))+b2*sin(2*wn(2)))^2)...
/(((1-(p+1)*cos(wn(2))+p*cos(2*wn(2)))^2+((p+1)*sin(wn(2))-p*sin(2*wn(2)))^2)))- 1/wn(2))/(1/wn(2)),...   
(sqrt(((b0+b1*cos(wn(3))+b2*cos(2*wn(3)))^2+(b1*sin(wn(3))+b2*sin(2*wn(3)))^2)...
/(((1-(p+1)*cos(wn(3))+p*cos(2*wn(3)))^2+((p+1)*sin(wn(3))-p*sin(2*wn(3)))^2)))- 1/wn(3))/(1/wn(3)),... 
(sqrt(((b0+b1*cos(wn(4))+b2*cos(2*wn(4)))^2+(b1*sin(wn(4))+b2*sin(2*wn(4)))^2)...
/(((1-(p+1)*cos(wn(4))+p*cos(2*wn(4)))^2+((p+1)*sin(wn(4))-p*sin(2*wn(4)))^2)))- 1/wn(4))/(1/wn(4)))


 
 b0d=double(S.b0)
 b1d=double(S.b1)
 b2d=double(S.b2)
 pd=double(S.p)
 
 
 %   samo prvih 8 od 16 resenja imaju p unutar jedinicnog kruga sto je pol
 %   prenosne funkcije integratora. Treba za svih 8 funkcija nacrtati gde
 %   su nule i polovi i pogledati gresku amplitude da bi se odabralo
 %   "najbolje" resenje
 
 for i=1:1
     nule=roots([b0d(i) b1d(i) b2d(i)])
     polovi=[1 ;pd(i)]  
     figure
     zplane(nule,polovi)
     [h,www]=freqz([b0d(i) b1d(i) b2d(i)],[1 -(pd(i)+1) pd(i)],1000);
     figure
     plot(www,www.*(abs(h)-1./www),'b','LineWidth',3)
     figure
     plot(www,(abs(h)-1./www),'r','LineWidth',3)
     grid  
     
     tse = 0;
     k = 0;
     b = [b0d(i) b1d(i) b2d(i)]
     p = pd(i)

    for w=pi/1000:pi/1000:pi
        k = k + 1;
        Hid(k) = 1 / w;

        %top = b(1)^2 + b(2)^2 + b(3)^2 + 2*( b(2)*(b(1)+b(3))*cos(w) + b(1)*b(3)*cos(2*w))
        %bottom = 2 * (p^2 + p + 1 - cos(w)*(1+p)^2 + p*cos(2*w))

        top = (b(1) + b(2)*cos(w) + b(3)*cos(2*w))^2 + (b(2)*sin(w) + b(3)*sin(2*w))^2;
        bottom = (1 - (p + 1)*cos(w) + p*cos(2*w))^2 + ((p + 1)*sin(w) - p*sin(2 * w))^2;
        
        H(k) = sqrt(top/bottom);
        err(k) = H(k) - Hid(k);
        tse = tse + err(k)^2;
    end 

    tse
    %pause
 end