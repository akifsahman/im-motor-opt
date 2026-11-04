function Kf = Kst_Kf(Kst)
%KST_KF  Kst -> Kf (form factor) table interpolation (range-checked).
%   Kst is an optimization variable; at extreme design points it may fall
%   out of range. In the original code, when out of range 'Kf_ekseni' was
%   zeroed and the function returned with Kf never assigned ("Output
%   argument not assigned" error). Here a defined error (imDesign:KstKf:*)
%   is thrown instead.
%
%   NOTE: along the curve Kf DECREASES as Kst increases (negative slope);
%   interp1 handles the monotonically decreasing y-axis correctly.

    % Original convention: the table is defined on the Kst+1 axis.
    Kst = Kst + 1;

    Kst_ekseni = [1.00 1.13 1.30 1.52 1.76 2.04 2.34 2.65 2.97];
    Kf_ekseni  = [1.11 1.10 1.09 1.08 1.07 1.06 1.05 1.04 1.03];

    % --- Input validation ---
    if ~isscalar(Kst) || ~isfinite(Kst) || ~isreal(Kst)
        error('imDesign:KstKf:badInput', ...
            'Kst must be scalar and finite (given: %g).', Kst - 1);
    end

    % --- Range check (axis is on Kst+1) ---
    if Kst < Kst_ekseni(1) || Kst > Kst_ekseni(end)
        error('imDesign:KstKf:outOfRange', ...
            'Kst = %g is outside the table range [%g, %g].', ...
            Kst - 1, Kst_ekseni(1) - 1, Kst_ekseni(end) - 1);
    end

    % --- Linear interpolation (NO extrapolation) ---
    Kf = interp1(Kst_ekseni, Kf_ekseni, Kst, 'linear');

    if ~isfinite(Kf)
        error('imDesign:KstKf:interpFailed', ...
            'Interpolation did not produce a valid value (Kst = %g).', Kst - 1);
    end
end
