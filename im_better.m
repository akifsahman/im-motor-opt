function w = im_better(cv, cf, co, lv, lf, lo, isDeb)
%IM_BETTER  Is the candidate (viol cv, eff cf, score co) better than the
%   current best (lv, lf, lo)?
%   isDeb=true  -> Deb's rule (feasible first; among feasible, higher
%                  efficiency; among infeasible, smaller violation).
%   isDeb=false -> scalar penalty: co < lo.
    if isDeb
        fa=(cv<=0); fb=(lv<=0);
        if     fa && fb, w = (cf >= lf);
        elseif fa,       w = true;
        elseif fb,       w = false;
        else             w = (cv <= lv);
        end
    else
        w = (co < lo);
    end
end
