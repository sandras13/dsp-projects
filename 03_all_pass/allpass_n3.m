clear all
close all

b = [1.156069364161850  -1.156069364161850]
a = [1.000000000000000   0.222030960372211]

[H1,ww] = freqz(b, a, 1000);
faza1 = unwrap(angle(H1));
faza1_id = pi/2 - 0.5*ww;

xn = [200 500 700]
wn = xn/1000 * pi

N = 3

syms a1 a2 a3

for i=1:3
    top = a1 * sin(wn(i)) + a2 * sin(2*wn(i)) + a3 * sin(3*wn(i));
    bottom = a1 * cos(wn(i)) + a2 * cos(2*wn(i)) + a3 * cos(3*wn(i));

    faza2 = -N*wn(i) + 2*atan(top/(1 + bottom));
    fazaA = (-N*wn(i) + faza2)/2;

    faza_id = pi/2 - (N + 0.5)*wn(i);
    faza_uk = faza1(xn(i)) + fazaA;

    uslov(i) = faza_uk - faza_id == 0; 
end

S = solve(uslov);

a1d = double(S.a1)
a2d = double(S.a2)
a3d = double(S.a3)

a = [1 a1d a2d a3d]
b = [a(4) a(3) a(2) a(1)]

nule = roots(b)
polovi = roots(a)

figure
zplane(nule, polovi)
title('Nule i polovi - all pass')

[H2,ww] = freqz(b, a, 1000);
faza2 = unwrap(angle(H2));

fazaA = (-N*ww + faza2)/2;
Ha = cos((-N*ww - faza2)/2);

Huk = H1.*Ha;

faza_id = pi/2 - (N + 0.5)*ww;
faza_uk = faza1 + fazaA;

figure
plot(ww(2:1000), 180*(faza_uk(2:1000)-faza_id(2:1000))/pi)
title('Greska ukupne faze')
grid

figure
plot(ww, (abs(Huk)-ww)./ww)
title('Greska amplitude')
