clear all
close all

% load b, p, wn from starting_solution

polovi = [1 p]'
a = poly(polovi)

[h,www]=freqz(b, a, 1000);
figure
plot(www,www.*(abs(h)-1./www),'b','LineWidth',3)
title('Pocetno resenje')

nule = roots(b)
p = polovi(2)
figure
zplane(nule,polovi)
title('Nule i polovi na pocetku')

epsilon = 10e-9
broj_iteracija = 0
current_max = 1
eps_val = [0.1 0.5 0.5 0.5 1]

while current_max > epsilon
    for i=1:5
        top = (b(1) + b(2)*cos(wn(i)) + b(3)*cos(2*wn(i)))^2 + (b(2)*sin(wn(i)) + b(3)*sin(2 * wn(i)))^2;
        bottom = (1 - (p + 1)*cos(wn(i)) + p*cos(2*wn(i)))^2 + ((p + 1)*sin(wn(i)) - p*sin(2 * wn(i)))^2;

        df_db0(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom * 2 * (b(1) + b(2)*cos(wn(i)) + b(3)*cos(2*wn(i)));
    
        df_db1(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * (b(1) + b(2)*cos(wn(i)) + b(3)*cos(2*wn(i))) * (cos(wn(i))) ...
            + 2 * (b(2)*sin(wn(i)) + b(3)*sin(2 * wn(i))) * sin(wn(i)));

        df_db2(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * (b(1) + b(2)*cos(wn(i)) + b(3)*cos(2*wn(i))) * (cos(2*wn(i))) ...
            + 2 * (b(2)*sin(wn(i)) + b(3)*sin(2 * wn(i))) * sin(2*wn(i)));

        df_dp(i) =  - 1/(2*sqrt(top/bottom)) * top / bottom^2 ...
            * (2 * (1 - (p + 1)*cos(wn(i)) + p*cos(2*wn(i))) * (-cos(wn(i)) + cos(2*wn(i))) ...
            + 2 * ((p + 1)*sin(wn(i)) - p*sin(2 * wn(i))) * (-sin(wn(i)) + sin(2*wn(i))));

        ff(i) = -sqrt(top/bottom) + 1/wn(i);
        fe(i) = -eps_val(i)/wn(i) * (-1)^(1+i);
    end

    B = [df_db0' df_db1' df_db2' df_dp' fe'];
    C = ff';

    A = B \ C

    b(1) = b(1) + A(1);
    b(2) = b(2) + A(2);
    b(3) = b(3) + A(3);
    p = p + A(4);
    eps = A(5);

    current_max = max(abs(A(1:4)))
    broj_iteracija = broj_iteracija + 1
end

b
p
eps

a = [1 -(p+1) p]
[h,www]=freqz(b, a, 1000);
figure
plot(www,www.*(abs(h)-1./www),'g','LineWidth',3)
title(['Nakon ' ,num2str(broj_iteracija), ' iteracija'])

nule = roots(b)
polovi = roots(a)
figure
zplane(nule,polovi)
title('Nule i polovi na kraju')

faza=unwrap(angle(h));
faza_id=-pi/2+0.5*www;
figure
plot(www/pi,faza,'r','LineWidth',3)
title(['Faza posle' ,num2str(broj_iteracija), ' iteracija'])
grid

figure
plot(www(2:1000)/pi,180*(faza(2:1000)-faza_id(2:1000))/pi,'r','LineWidth',3)
title(['Greska faze posle' ,num2str(broj_iteracija), ' iteracija u stepenima'])
grid
