function [ekstr] = nadji_ekstremume(N, faza, cnt)
    i = 1;
    ekstr = [];
    for k = 3:cnt-2
        if (faza(k) > faza(k + 1)) && (faza(k) > faza(k - 1))
            ekstr(i) = k;
            i = i + 1;
        end

        if (faza(k) < faza(k + 1)) && (faza(k) < faza(k - 1))
            ekstr(i) = k;
            i = i + 1;
        end
    end
end