clear all
close all

% The error is set to zero at 2N distinct frequency points:
% wn = pi * [1 - 0.8*exp((-1.5*(k-1)*pi)/(2*N - 1))]

syms b0 b1 b2 p
j=sqrt(-1)

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
