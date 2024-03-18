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
% koef = fliplr(koef);
noise = cos(2*pi*freq.*tdd).*koef;

figure
plot([0:100], noise)
title('Signal greske')

td2 = td*1000;
xval = 0.001;

%% Pravougaona f-ja
% x(1:100) = -xval;
% x(101:200) = xval;

% x1(1:20) = -xval/5;
% x1(21:40) = xval/5;
% x = repmat(x1, 1, 5);

% x1(1:10) = -xval/10;
% x1(11:20) = xval/10;
% x = repmat(x1, 1, 10);

% x1(1:2) = -xval/50;
% x1(3:4) = xval/50;
% x = repmat(x1, 1, 50);

%% Linearna f-ja
% x(1:100) = -xval/4 + xval*td2(1:100);
% x(101:200) = xval*3/4 - xval*td2(101:200);

% x1(1:20) = -xval/20 + xval*td2(1:20);
% x1(21:40) = xval*3/20 - xval*td2(21:40);
% x = repmat(x1, 1, 5);

% x1(1:10) = -xval/40 + xval*td2(1:10);
% x1(11:20) = xval*3/40 - xval*td2(11:20);
% x = repmat(x1, 1, 10);

% x1(1:2) = -xval/200 + xval*td2(1:2);
% x1(3:4) = xval*3/200 - xval*td2(3:4);
% x = repmat(x1, 1, 50);

%% Kvadratna f-ja
% x(1:100) = 2*xval*td2(1:100).^2 - xval*td2(1:100);
% x(101:200) = -2*xval*td2(101:200).^2 + 3*xval*td2(101:200) - xval;

% x1(1:20) = 10*xval*td2(1:20).^2 - xval*td2(1:20);
% x1(21:40) = -10*xval*td2(21:40).^2 + 3*xval*td2(21:40) - xval/5;
% x = repmat(x1, 1, 5);

% x1(1:10) = 20*xval*td2(1:10).^2 - xval*td2(1:10);
% x1(11:20) = -20*xval*td2(11:20).^2 + 3*xval*td2(11:20) - xval/10;
% x = repmat(x1, 1, 10);

x1(1:2) = 100*xval*td2(1:2).^2 - xval*td2(1:2);
x1(3:4) = -100*xval*td2(3:4).^2 + 3*xval*td2(3:4) - xval/50;
x = repmat(x1, 1, 50);

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

Yid = odredi_odziv(noise, 0, 0, Hid, faza_id, X);
yid = real(ifft(Yid)) * fs;

figure
stem([0:199], yid)
grid
title('Izlaz idealnog diferencijatora')

k = 1;
for ii = 0:1/100:1
    gre_amp = ii * max_err_rel;
    gre_faz = ii * max_err_ph;
    
    Y1 = odredi_odziv(noise, gre_amp, 0, Hid, faza_id, X);
    y1 = real(ifft(Y1)) * fs;

    Y2 = odredi_odziv(noise, 0, gre_faz, Hid, faza_id, X);
    y2 = real(ifft(Y2)) * fs;

    eps1(k) = 1/Nt * sum(abs(y1 - yid));
    eps2(k) =  1/Nt * sum(abs(y2 - yid));

    delta(k) = ii;
    k = k + 1;
end

figure
plot(delta, eps1, delta, eps2)
title('Prosecno odstupanje od idealnog izlaza')
legend('amplituda', 'faza')

amp_greska = delta*max_err_rel*100;
fazna_greska = delta(eps2<max(eps1))*max_err_ph_stepeni;

size_diff = length(amp_greska) - length(fazna_greska);
t = linspace(1, length(fazna_greska), length(fazna_greska) + size_diff);
fazna_greska_interp = interp1(1:length(fazna_greska), fazna_greska, t, 'linear');

figure
plot(amp_greska, fazna_greska_interp);
xlabel('rel. magn. error [%]')
ylabel('phase error [deg]')