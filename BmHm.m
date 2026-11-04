function aki_siddeti = BmHm(aki_yogunlugu)
%BMHM  B-H curve (M19-24G) table interpolation (range-checked).
%   If the flux density is OUTSIDE the table range, a defined error
%   (imDesign:BmHm:*) is thrown instead of returning 0. Because the inputs
%   of this function (Bts, Btr, Bcr, Bcs, Bg2) are derived from the
%   optimization variables, this is the function most likely to go
%   out of range at extreme design points.

    % --- M19-24G B-H table ---
    aki_yogunlugu_skalasi = [0.00 0.05 0.10 0.15 0.36 0.54 0.65 0.99 1.20 ...
        1.28 1.33 1.36 1.44 1.52 1.58 1.63 1.67 1.80 1.90 2.00 2.10 2.30 2.50];

    aki_siddeti_skalasi = [0.00 22.28 25.46 31.83 47.74 63.66 79.57 159.15 ...
        318.30 477.46 636.61 795.77 1591.50 3183.00 4774.60 6366.10 7957.70 ...
        15915.00 31830.00 111407.00 190984.00 350138.00 509292.00];

    % --- Input validation ---
    if ~isscalar(aki_yogunlugu) || ~isfinite(aki_yogunlugu) || ~isreal(aki_yogunlugu)
        error('imDesign:BmHm:badInput', ...
            'Flux density must be scalar and finite (given: %g).', aki_yogunlugu);
    end

    % --- Range check ---
    if aki_yogunlugu < aki_yogunlugu_skalasi(1) || ...
       aki_yogunlugu > aki_yogunlugu_skalasi(end)
        error('imDesign:BmHm:outOfRange', ...
            'B = %g T is outside the table range [%g, %g].', ...
            aki_yogunlugu, aki_yogunlugu_skalasi(1), aki_yogunlugu_skalasi(end));
    end

    % --- Linear interpolation (NO extrapolation) ---
    aki_siddeti = interp1(aki_yogunlugu_skalasi, aki_siddeti_skalasi, ...
        aki_yogunlugu, 'linear');

    if ~isfinite(aki_siddeti)
        error('imDesign:BmHm:interpFailed', ...
            'Interpolation did not produce a valid value (B = %g).', aki_yogunlugu);
    end
end
