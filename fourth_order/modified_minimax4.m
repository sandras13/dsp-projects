clear all
close all

% load b, p, wn from starting_solution
load N4_START_SOL

polovi = [1 p]'

a = poly(polovi)

[h,www]=freqz(b, a, 1000);
figure
plot(www,www.*(abs(h)-1./www),'b','LineWidth',3)
title('Pocetno resenje')

nule = roots(b)
figure
zplane(nule,polovi)
title('Nule i polovi na pocetku')

faza=unwrap(angle(h));
faza_id=-pi/2+0.5*www;
figure
plot(www/pi,faza,'r','LineWidth',3)
title('Faza na pocetku')
grid

epsilon = 10e-9
broj_iteracija = 0
current_max = 1
eps_val = [0.1 0.3 0.3 0.3 0.5 0.5 0.5 0.5 1]

while current_max > epsilon
    for i=1:9
        top1 = b(1) + b(2)*cos(wn(i)) + b(3)*cos(2*wn(i)) + b(4)*cos(3*wn(i)) + b(5)*cos(4*wn(i));
        top2 = b(2)*sin(wn(i)) + b(3)*sin(2*wn(i)) + b(4)*sin(3*wn(i)) + b(5)*sin(4*wn(i));
        
        bottom1 = 1 - (p(1) + p(2) + p(3) + 1)*cos(wn(i)) + (p(1)*p(2) + p(2)*p(3)+ p(3)*p(1) + p(1) + p(2) + p(3))*cos(2*wn(i)) ...
        - (p(1)*p(2) + p(1)*p(2)*p(3) + p(1)*p(3) + p(2)*p(3))*cos(3*wn(i)) + p(1)*p(2)*p(3)*cos(4*wn(i));

        bottom2 = (p(1) + p(2) + p(3) + 1)*sin(wn(i)) - (p(1)*p(2) + p(2)*p(3)+ p(3)*p(1) + p(1) + p(2) + p(3))*sin(2*wn(i)) ...
        + (p(1)*p(2) + p(1)*p(2)*p(3) + p(1)*p(3) + p(2)*p(3))*sin(3*wn(i)) - p(1)*p(2)*p(3)*sin(4*wn(i));

        top = top1^2 + top2^2;
        bottom = bottom1^2 + bottom2^2;

        df_db0(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom * 2 * top1;

        df_db1(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(wn(i)) + 2 * top2 * sin(wn(i)));

        df_db2(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(2*wn(i)) + 2 * top2 * sin(2*wn(i)));

        df_db3(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(3*wn(i)) + 2 * top2 * sin(3*wn(i)));

        df_db4(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(4*wn(i)) + 2 * top2 * sin(4*wn(i)));

        df_dp1(i) =  - 1/(2*sqrt(top/bottom)) * top / bottom^2 ...
            * (2 * bottom1 * (-cos(wn(i)) + (1 + p(2) + p(3)) * cos(2*wn(i))...
            - (p(2) + p(3) + p(2)*p(3)) * cos(3*wn(i))...
            + p(2)*p(3)*cos(4*wn(i)))...
            + 2 * bottom2 * (sin(wn(i)) - (1 + p(2) + p(3)) * sin(2*wn(i))...
            + (p(2) + p(3) + p(2)*p(3)) * sin(3*wn(i))...
            - p(2)*p(3)*sin(4*wn(i))));

        df_dp2(i) =  - 1/(2*sqrt(top/bottom)) * top / bottom^2 ...
            * (2 * bottom1 * (-cos(wn(i)) + (1 + p(1) + p(3)) * cos(2*wn(i))...
            - (p(1) + p(3) + p(1)*p(3)) * cos(3*wn(i))...
            + p(1)*p(3)*cos(4*wn(i)))...
            + 2 * bottom2 * (sin(wn(i)) - (1 + p(1) + p(3)) * sin(2*wn(i))...
            + (p(1) + p(3) + p(1)*p(3)) * sin(3*wn(i))...
            - p(1)*p(3)*sin(4*wn(i))));

        df_dp3(i) =  - 1/(2*sqrt(top/bottom)) * top / bottom^2 ...
            * (2 * bottom1 * (-cos(wn(i)) + (1 + p(2) + p(1)) * cos(2*wn(i))...
            - (p(2) + p(1) + p(2)*p(1)) * cos(3*wn(i))...
            + p(2)*p(1)*cos(4*wn(i)))...
            + 2 * bottom2 * (sin(wn(i)) - (1 + p(2) + p(1)) * sin(2*wn(i))...
            + (p(2) + p(1) + p(2)*p(1)) * sin(3*wn(i))...
            - p(2)*p(1)*sin(4*wn(i))));

        ff(i) = -sqrt(top/bottom) + 1/wn(i);
        fe(i) = -eps_val(i)/wn(i) * (-1)^(1+i);
    end

    B = [df_db0' df_db1' df_db2' df_db3' df_db4' df_dp1' df_dp2' df_dp3' fe']
    C = ff'

    A = B \ C

    b(1) = b(1) + A(1);
    b(2) = b(2) + A(2);
    b(3) = b(3) + A(3);
    b(4) = b(4) + A(4);
    b(5) = b(5) + A(5);
    p(1) = p(1) + A(6);
    p(2) = p(2) + A(7);
    p(3) = p(3) + A(8);
    eps = A(9);

    current_max = max(abs(A(1:8)))
    broj_iteracija = broj_iteracija + 1
end

b
p
eps

polovi = [1 p]'

a = poly(polovi)
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
title(['Greska faze posle ' ,num2str(broj_iteracija), ' iteracija u stepenima'])
grid
