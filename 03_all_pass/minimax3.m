clear all
close all

b = [1.156069364161850  -1.156069364161850]
a = [1.000000000000000   0.222030960372211]

pts = 50000
[H1,ww] = freqz(b, a, pts);
faza1 = unwrap(angle(H1));

faza1_id = pi/2 - 0.5*ww;

N = 3
broj_iteracija = 30

a1 = [-0.222799753747073   0.049599932168010  -0.010362718537241]

for cnt=1:broj_iteracija
    a = [1 a1];
    b = [a(4) a(3) a(2) a(1)];
    [H2,ww] = freqz(b, a, pts);
    faza2 = unwrap(angle(H2));

    fazaA = (-N*ww + faza2)/2;

    faza_id = pi/2 - (N + 0.5)*ww;
    faza_uk = faza1 + fazaA;
    ffaz = faza_uk - faza_id;

    xn = nadji_ekstremume(N, ffaz, pts)
    wn = xn/pts * pi

    for i=1:4
        top = a1(1) * sin(wn(i)) + a1(2) * sin(2*wn(i)) + a1(3) * sin(3*wn(i));
        bottom = a1(1) * cos(wn(i)) + a1(2) * cos(2*wn(i)) + a1(3) * cos(3*wn(i));

        faza2 = -N*wn(i) + 2*atan(top/(1 + bottom));
        fazaA = (-N*wn(i) + faza2)/2;

        faza_id = pi/2 - (N + 0.5)*wn(i);
        faza_uk = faza1(xn(i)) + fazaA;
        
        f11 = sin(wn(i)) * (1 + a1(1) * cos(wn(i)) + a1(2) * cos(2*wn(i)) + a1(3) * cos(3*wn(i)))...
            - (a1(1)*sin(wn(i)) + a1(2) * sin(2*wn(i)) + a1(3) * sin(3*wn(i)))*cos(wn(i));

        f12 = (1 + a1(1) * cos(wn(i)) + a1(2) * cos(2*wn(i)) + a1(3) * cos(3*wn(i)))^2;

        f21 = sin(2*wn(i)) * (1 + a1(1) * cos(wn(i)) + a1(2) * cos(2*wn(i)) + a1(3) * cos(3*wn(i)))...
            - (a1(1)*sin(wn(i)) + a1(2) * sin(2*wn(i)) + a1(3) * sin(3*wn(i)))*cos(2*wn(i));

        f31 = sin(3*wn(i)) * (1 + a1(1) * cos(wn(i)) + a1(2) * cos(2*wn(i)) + a1(3) * cos(3*wn(i)))...
            - (a1(1)*sin(wn(i)) + a1(2) * sin(2*wn(i)) + a1(3) * sin(3*wn(i)))*cos(3*wn(i));

        df_da1(i) = 2/(1 + (top/(1+bottom))^2) * (f11/f12);
        df_da2(i) = 2/(1 + (top/(1+bottom))^2) * (f21/f12);
        df_da3(i) = 2/(1 + (top/(1+bottom))^2) * (f31/f12);

        ff(i) = -faza_uk + faza_id;
        fe(i) = -(-1)^(i+1);
    end

    B = [df_da1' df_da2' df_da3' fe'];
    C = ff';

    A = B \ C

    a1(1) = a1(1) + A(1);
    a1(2) = a1(2) + A(2);
    a1(3) = a1(3) + A(3);
    eps = A(4);
end

a = [1 a1]
b = [a(4) a(3) a(2) a(1)];

[H2,ww] = freqz(b, a, pts);
faza2 = unwrap(angle(H2));

fazaA = (-N*ww + faza2)/2;

faza_id = pi/2 - (N + 0.5)*ww;
faza_uk = faza1 + fazaA;

figure
plot(ww(2:pts), 180*(faza_uk(2:pts)-faza_id(2:pts))/pi)
title('Greska ukupne faze')
grid