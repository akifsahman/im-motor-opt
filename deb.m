function [viol, verim, x] = deb(viol_a, verim_a, x_a, viol_b, verim_b, x_b)
%DEB  Deb's constraint-handling rule (feasibility-based selection).
%   Selects the "better" of two candidates. Objective: MAXIMIZE efficiency,
%   MINIMIZE violation.
%   viol_* : total constraint violation (>=0, feasible if 0)
%   verim_*: efficiency (%), x_* : design vector
%
%   Rule (Deb 2000):
%     1) Both feasible    -> the one with higher efficiency
%     2) One feasible     -> the feasible one
%     3) Both infeasible  -> the one with smaller violation
    fa = (viol_a <= 0);
    fb = (viol_b <= 0);
    if fa && fb
        sec_a = (verim_a >= verim_b);      % higher efficiency wins
    elseif fa && ~fb
        sec_a = true;
    elseif ~fa && fb
        sec_a = false;
    else
        sec_a = (viol_a <= viol_b);        % smaller violation wins
    end
    if sec_a
        viol = viol_a; verim = verim_a; x = x_a;
    else
        viol = viol_b; verim = verim_b; x = x_b;
    end
end
