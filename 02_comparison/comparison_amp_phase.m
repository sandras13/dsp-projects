clear all
close all

fs = 200000;
Ts = 1/fs;

td = [0:199] * Ts;
Nt = length(td);

max_err_rel = 0.1;
max_err_ph_stepeni = 20;
max_err_ph = max_err_ph_stepeni / 180 * pi;
np = 4; % broj ekstremuma - 1

diff = 0.05;
num_periods = diff:(np-diff)/100:np;

freq = 1000 * num_periods;
tdd = [0:100] * Ts;

min_val = 0.1;
koef = [min_val:(1-min_val)/100:1];
koef_flip = fliplr(koef);

noise_amp = cos(2*pi*freq.*tdd).*koef;
noise_ph = cos(2*pi*freq.*tdd).*koef_flip;

figure
plot([0:100], noise_amp)
title('Signal greske amplitude')

figure
plot([0:100], noise_ph)
title('Signal greske faze')

td2 = td*1000;
xval = 0.001;

% Kvadratna f-ja
x(1:100) = 2*xval*td2(1:100).^2 - xval*td2(1:100);
x(101:200) = -2*xval*td2(101:200).^2 + 3*xval*td2(101:200) - xval;

X = fft(x);

figure
stem([0:Nt-1], x)
grid
title('Ulazni signal')

figure 
stem([0:Nt-1],abs(X));
grid
title('Spektar ulaznog signala')

ww = pi*[0:100]/100;
Hid = ww;
faza_id = pi/2;

Yid = odredi_odziv(noise_amp, 0, 0, Hid, faza_id, X);
yid = real(ifft(Yid)) * fs;

figure
stem([0:199], yid)
grid
title('Izlaz idealnog diferencijatora')

k = 1;
for ii = 0:1/100:1
    gre_amp = ii * max_err_rel;
    gre_faz = ii * max_err_ph;
    
    Y1 = odredi_odziv_razl(noise_amp, noise_ph, gre_amp, gre_faz, Hid, faza_id, X);
    y1 = real(ifft(Y1)) * fs;

    eps1(k) = 1/Nt * sum(abs(y1 - yid));

    delta(k) = ii;
    k = k + 1;
end

figure
plot(delta, eps1)
title('Prosecno odstupanje od idealnog izlaza')