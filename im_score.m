function Sc = im_score(Ver, Viol, Obj, isDeb)
%IM_SKOR  In-population ranking score (lower = better).
%   isDeb -> Deb surrogate (feasible first); otherwise -> penalty score Obj.
    if isDeb, Sc = (Viol>1e-9)*1e6 + Viol*100 - Ver; else, Sc = Obj; end
end
