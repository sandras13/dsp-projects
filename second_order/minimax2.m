clear all
close all

% Opcija 1
% b = [-0.8631   -0.6350   -0.0580]
% p =  -0.5560
% wn = [0.001 1.2409 2.3185 2.9373 3.1385]

% Opcija 2
b = [0.8634    0.6883    0.0718]
p = -0.6166
wn = [0.001 1.5833 2.5541 2.9907 3.1385]


a = [1 -(p+1) p]
[h,www]=freqz(b, a, 1000);
figure
plot(www/pi,www.*(abs(h)-1./www),'b','LineWidth',3)
title('Pocetno resenje')

nule = roots(b)
polovi = roots(a)
figure
zplane(nule,polovi)
title('Nule i polovi na pocetku')

% (H - Hid) / Hid =  eps * (-1)^(i + 1)
% H - Hid = Hid *  eps * (-1)^(i + 1)

% ff + df_db0 * delta_b0 + df_db1 * delta_b1 + df_db2 * delta_b2
% + df_dp * delta_p - 1/w = 1/w * eps * (-1)^(i + 1)

% df_db0 * delta_b0 + df_db1 * delta_b1 + df_db2 * delta_b2
% + df_dp * delta_p - 1/w * eps * (-1)^(i + 1) = -ff + 1 /w

broj_iteracija = 20

for cnt=1:broj_iteracija

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
        fe(i) = -1/wn(i) * (-1)^(1+i);
    end

    B = [df_db0' df_db1' df_db2' df_dp' fe'];
    C = ff';

    A = B \ C % [delta_b0 delta_b1 delta_b2 delta_p eps]'

    b(1) = b(1) + A(1);
    b(2) = b(2) + A(2);
    b(3) = b(3) + A(3);
    p = p + A(4);
    eps = A(5);

end

b
p
eps

a = [1 -(p+1) p]
[h,www]=freqz(b, a, 1000);
figure
plot(www/pi,www.*(abs(h)-1./www),'g','LineWidth',3)
title(['Nakon ' ,num2str(broj_iteracija), ' iteracija'])

nule = roots(b)
polovi = roots(a)
figure
zplane(nule,polovi)
title('Nule i polovi na kraju')