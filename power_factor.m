function pf = power_factor(p1_cift, Pn)
%GUC_FAKTORU  Power-factor table interpolation (range-checked).
%   If the input is OUTSIDE the table range, a defined error
%   (imDesign:gucFaktoru:*) is thrown instead of returning 0; the try/catch
%   in the optimization wrapper can catch it and assign a large penalty.

    % --- Tables ---
    kutup_sayisi_skalasi = [0.8600 0.7600; 0.8750 0.7850; 0.8820 0.8000;
                            0.8850 0.8080; 0.8870 0.8130; 0.8875 0.8180;
                            0.8880 0.8220; 0.8885 0.8240; 0.8890 0.8260;
                            0.8895 0.8290; 0.8900 0.8310; 0.8953 0.8450;
                            0.8977 0.8520; 0.9000 0.8573; 0.9020 0.8600;
                            0.9040 0.8617; 0.9060 0.8633; 0.9065 0.8647;
                            0.9070 0.8657; 0.9075 0.8663];

    guc_skalasi = [0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 ...
        20000 30000 40000 50000 60000 70000 80000 90000 100000];

    % --- Input validation ---
    if ~isscalar(Pn) || ~isfinite(Pn) || ~isreal(Pn)
        error('imDesign:gucFaktoru:badInput', ...
            'Pn must be scalar and finite (given: %g).', Pn);
    end

    nCol = size(kutup_sayisi_skalasi, 2);          % = 2
    if ~ismember(p1_cift, 1:nCol)
        error('imDesign:gucFaktoru:badColumn', ...
            'p1_cift must be in 1..%d (given: %g).', nCol, p1_cift);
    end

    % --- Range check (preserves the original 0/negative/upper-bound behavior) ---
    if Pn <= guc_skalasi(1) || Pn > guc_skalasi(end)
        error('imDesign:gucFaktoru:outOfRange', ...
            'Pn = %g W is outside the table range (%g, %g].', ...
            Pn, guc_skalasi(1), guc_skalasi(end));
    end

    % --- Linear interpolation (NO extrapolation -> out of range gives NaN) ---
    pf = interp1(guc_skalasi, kutup_sayisi_skalasi(:, p1_cift), Pn, 'linear');

    if ~isfinite(pf)
        error('imDesign:gucFaktoru:interpFailed', ...
            'Interpolation did not produce a valid value (Pn = %g).', Pn);
    end
end
