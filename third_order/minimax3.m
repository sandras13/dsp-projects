clear all
close all
% 
% b = [0.8655    1.2191    0.4404    0.0283]
% wn = [0.001 1.3504  2.331 2.75  3.008 3.115 3.13845]
% p = [-0.8535 -0.3757]

b = [0.071167882067755   0.960602064542676   1.186928071675837   0.352874433598843]
p = [ -0.859902859737884 -0.381044315703754]
wn = [0.001 1.21265 2.07659 2.62323 2.92482 3.08819 3.13845]

b = [.079283883661398   0.978350937672573   1.271248924983076   0.407305762269729]
p = [-0.902352690780659  -0.436304953393463]
wn = [0.001 1.30062 2.2054 2.73319 3.00336 3.11018 3.13845]

b = [0.083900277592316   0.986912408210583   1.295911111218375   0.425364823694206]
p = [-0.905405539537736  -0.462458150666399]
wn = [0.001 1.35403 2.29022 2.78659 3.01278 3.10704 3.13845]

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

% (H - Hid) / Hid =  eps * (-1)^(i + 1)
% H - Hid = Hid *  eps * (-1)^(i + 1)

% ff + df_db0 * delta_b0 + df_db1 * delta_b1 + df_db2 * delta_b2
% + df_db3 * delta_b3 + df_dp1 * delta_p1 + df_dp2 * delta_p2
% - 1/w = 1/w * eps * (-1)^(i + 1)

% df_db0 * delta_b0 + df_db1 * delta_b1 + df_db2 * delta_b2
% + df_db3 * delta_b3 + df_dp1 * delta_p1 + df_dp2 * delta_p2
% - 1/w * eps * (-1)^(i + 1) = -ff + 1 /w

broj_iteracija = 20

for cnt=1:broj_iteracija

    for i=1:7
        top1 = b(1) + b(2)*cos(wn(i)) + b(3)*cos(2*wn(i)) + b(4)*cos(3*wn(i));
        top2 = b(2)*sin(wn(i)) + b(3)*sin(2*wn(i)) + b(4)*sin(3*wn(i));

        bottom1 = 1 - (1 + p(1) + p(2))*cos(wn(i)) + (p(1) + p(2) + p(1)*p(2))*cos(2*wn(i)) ...
                  - p(1) * p(2) * cos(3*wn(i));

        bottom2 = (1 + p(1) + p(2))*sin(wn(i)) - (p(1) + p(2) + p(1)*p(2))*sin(2*wn(i)) ...
                  + p(1) * p(2) * sin(3*wn(i));

        top = top1^2 + top2^2;
        bottom = bottom1^2 + bottom2^2;

        df_db0(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom * 2 * top1;

        df_db1(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(wn(i)) + 2 * top2 * sin(wn(i)));

        df_db2(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(2*wn(i)) + 2 * top2 * sin(2*wn(i)));

        df_db3(i) = 1/(2*sqrt(top/bottom)) * 1 / bottom ...
            * (2 * top1 * cos(3*wn(i)) + 2 * top2 * sin(3*wn(i)));

        df_dp1(i) =  - 1/(2*sqrt(top/bottom)) * top / bottom^2 ...
            * (2 * bottom1 * (-cos(wn(i)) + (1 + p(2)) * cos(2*wn(i))...
            - p(2) * cos(3*wn(i)))...
            + 2 * bottom2 * (sin(wn(i)) - (1 + p(2)) * sin(2*wn(i))...
            + p(2) * sin(3*wn(i))));

        df_dp2(i) =  - 1/(2*sqrt(top/bottom)) * top / bottom^2 ...
            * (2 * bottom1 * (-cos(wn(i)) + (1 + p(1)) * cos(2*wn(i))...
            - p(1) * cos(3*wn(i)))...
            + 2 * bottom2 * (sin(wn(i)) - (1 + p(1)) * sin(2*wn(i))...
            + p(1) * sin(3*wn(i))));

        ff(i) = -sqrt(top/bottom) + 1/wn(i);
        fe(i) = -1/wn(i) * (-1)^(1+i);
    end

    B = [df_db0' df_db1' df_db2' df_db3' df_dp1' df_dp2' fe'];
    C = ff';

    A = B \ C; % [delta_b0 delta_b1 delta_b2 delta_b3 delta_p1 delta_p2 eps]'

    b(1) = b(1) + A(1);
    b(2) = b(2) + A(2);
    b(3) = b(3) + A(3);
    b(4) = b(4) + A(4);
    p(1) = p(1) + A(5);
    p(2) = p(2) + A(6);
    eps = A(7);

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
faza_id=-pi/2-0.5*www;
figure
plot(www/pi,faza,'r','LineWidth',3)
title(['Faza posle' ,num2str(broj_iteracija), ' iteracija'])
grid

figure
plot(www(2:1000)/pi,180*(faza(2:1000)-faza_id(2:1000))/pi,'r','LineWidth',3)
title(['Greska faze posle ' ,num2str(broj_iteracija), ' iteracija u stepenima'])
grid