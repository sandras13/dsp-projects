function [Y] = odredi_odziv_razl(korekc_amp, korekc_ph, amp_err, ph_err, Hid, faza_id, X)
    H = (1 + amp_err*korekc_amp) .* Hid;
    faza = faza_id + korekc_ph*ph_err;

    Hodz = H .* exp(i*faza);
    Y = [Hodz conj(fliplr(Hodz(2:100)))] .* X;
end