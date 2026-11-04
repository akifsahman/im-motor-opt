function esson = esson_constant(p1_cift, Pn)
%ESSON_SABITI  Esson constant (C0) table interpolation (range-checked).
%   Uses the ORIGINAL 33-row / 6-column table (guc_skalasi 100..50000).
%   Instead of returning 0 when out of range, a defined error
%   (imDesign:esson:*) is thrown. The signature (p1_cift, Pn) and the
%   output scaling are unchanged, so no changes are needed in callers.

    kutup_sayisi_skalasi = [0.90 1.00 1.03 1.05 1.08 1.10;
                            0.93 1.03 1.05 1.08 1.10 1.13;
                            0.95 1.05 1.08 1.10 1.13 1.15;
                            0.98 1.08 1.10 1.13 1.15 1.18;
                            1.00 1.10 1.13 1.15 1.18 1.20;
                            1.06 1.14 1.17 1.20 1.23 1.26;
                            1.12 1.18 1.22 1.25 1.29 1.32;
                            1.18 1.22 1.26 1.30 1.34 1.38;
                            1.24 1.26 1.31 1.35 1.40 1.44;
                            1.30 1.30 1.35 1.40 1.45 1.50;
                            1.33 1.40 1.50 1.60 1.70 1.80;
                            1.36 1.50 1.63 1.75 1.88 2.00;
                            1.40 1.70 1.85 2.00 2.15 2.30;
                            1.46 1.90 2.03 2.16 2.28 2.41;
                            1.53 1.98 2.12 2.25 2.39 2.53;
                            1.59 2.06 2.20 2.35 2.49 2.64;
                            1.65 2.14 2.29 2.45 2.60 2.75;
                            1.71 2.22 2.38 2.54 2.70 2.86;
                            1.78 2.30 2.47 2.64 2.81 2.98;
                            1.84 2.55 2.68 2.82 2.95 3.09;
                            1.90 2.80 2.90 3.00 3.10 3.20;
                            1.95 3.00 3.08 3.16 3.24 3.32;
                            2.00 3.10 3.18 3.27 3.35 3.43;
                            2.01 3.15 3.25 3.29 3.38 3.49;
                            2.05 3.18 3.28 3.33 3.56 3.54;
                            2.06 3.21 3.30 3.38 3.59 3.58;
                            2.07 3.25 3.34 3.42 3.59 3.59;
                            2.10 3.30 3.39 3.48 3.59 3.67;
                            2.11 3.35 3.48 3.53 3.62 3.71;
                            2.13 3.40 3.50 3.59 3.68 3.78;
                            2.14 3.48 3.55 3.61 3.70 3.80;
                            2.15 3.49 3.58 3.68 3.79 3.86;
                            2.17 3.50 3.60 3.70 3.80 3.90];

    guc_skalasi = [100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 ...
        6000 7000 8000 9000 10000 15000 20000 25000 30000 32000 34000 36000 37000 ...
        40000 42000 45000 46000 48000 50000];

    if ~isscalar(Pn) || ~isfinite(Pn) || ~isreal(Pn)
        error('imDesign:esson:badInput','Pn must be scalar and finite (given: %g).', Pn);
    end
    nCol = size(kutup_sayisi_skalasi,2);          % = 6
    if ~ismember(p1_cift, 1:nCol)
        error('imDesign:esson:badColumn','p1_cift must be in 1..%d (given: %g).', nCol, p1_cift);
    end
    if Pn < guc_skalasi(1) || Pn > guc_skalasi(end)
        error('imDesign:esson:outOfRange','Pn = %g is outside the table range [%g, %g].', ...
            Pn, guc_skalasi(1), guc_skalasi(end));
    end
    tablo = interp1(guc_skalasi, kutup_sayisi_skalasi(:,p1_cift), Pn, 'linear');
    esson = 10^3*60*tablo;
    if ~isfinite(esson)
        error('imDesign:esson:interpFailed','Interpolation did not produce a valid value (Pn = %g).', Pn);
    end
end
