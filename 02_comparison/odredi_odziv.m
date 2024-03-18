function [Y] = odredi_odziv(korekc, amp_err, ph_err, Hid, faza_id, X)
    H = (1 + amp_err*korekc) .* Hid;
    faza = faza_id + korekc*ph_err;

    Hodz = H .* exp(i*faza);
    Y = [Hodz conj(fliplr(Hodz(2:100)))] .* X;
end