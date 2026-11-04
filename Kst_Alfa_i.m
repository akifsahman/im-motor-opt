function Alfa_i = Kst_Alfa_i(Kst)
%KST_ALFA_I  Kst -> Alfa_i (pole-arc coefficient) interpolation (range-checked).
%   Kst is an optimization variable; at extreme design points it may fall
%   out of range. In the original code, when out of range 'Alfa_i_ekseni'
%   was zeroed and the function returned with Alfa_i never assigned ("Output
%   argument not assigned" error). Here a defined error (imDesign:KstAlfaI:*)
%   is thrown instead.

    % Original convention: the table is defined on the Kst+1 axis.
    Kst = Kst + 1;

    Kst_ekseni = [1.00 1.03 1.07 1.11 1.15 1.20 1.25 1.30 1.35 1.42 1.49 ...
        1.56 1.66 1.77 1.90 2.05 2.23 2.42 2.63 2.85 3.00];

    Alfa_i_ekseni = [0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 ...
        0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84];

    % --- Input validation ---
    if ~isscalar(Kst) || ~isfinite(Kst) || ~isreal(Kst)
        error('imDesign:KstAlfaI:badInput', ...
            'Kst must be scalar and finite (given: %g).', Kst - 1);
    end

    % --- Range check (axis is on Kst+1) ---
    if Kst < Kst_ekseni(1) || Kst > Kst_ekseni(end)
        error('imDesign:KstAlfaI:outOfRange', ...
            'Kst = %g is outside the table range [%g, %g].', ...
            Kst - 1, Kst_ekseni(1) - 1, Kst_ekseni(end) - 1);
    end

    % --- Linear interpolation (NO extrapolation) ---
    Alfa_i = interp1(Kst_ekseni, Alfa_i_ekseni, Kst, 'linear');

    if ~isfinite(Alfa_i)
        error('imDesign:KstAlfaI:interpFailed', ...
            'Interpolation did not produce a valid value (Kst = %g).', Kst - 1);
    end
end
