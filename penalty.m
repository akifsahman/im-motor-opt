function [toplam, tviol] = penalty(verim, gviol, method, iter, maxiter)
%CEZALAR  Single hub for 8 different penalty (constraint-handling) approaches,
%   including Deb's rule.
%
%   Single objective: MAXIMIZE eff_h (efficiency). For scalar methods the
%   minimized value is: toplam = -verim + penalty.
%
%   INPUT:
%     verim   : im_core output (%), to be maximized
%     gviol   : 1xNc constraint-violation vector (each element >=0, 0=satisfied)
%     method  : 1..8 (see below)
%     iter    : current outer iteration (for dynamic methods), 1..maxiter
%     maxiter : total number of outer iterations (= max_fes/n)
%
%   OUTPUT:
%     toplam  : scalar to be minimized (NaN for method 1/Deb; driver uses deb())
%     tviol   : raw total violation (for feasibility and Deb) = sum(gviol)
%
%   METHODS:
%     1  Deb's rule               (selection via deb.m; toplam unused)
%     2  Static constant penalty  -verim + K*(#violated constraints)
%     3  Violation sum (linear)   -verim + sum(gviol)
%     4  Sum of squared viol.     -verim + sum(gviol.^2)
%     5  Normalized viol. count   -verim + K*(#violated/Nc)
%     6  Joines-Houck (dynamic)   -verim + (C*iter)^alpha * sum(gviol.^beta)
%     7  Liu et al. (dyn, linear) -verim + mult(iter) * sum(gviol)
%     8  Liu et al. (dyn, square) -verim + mult(iter) * sum(gviol.^2)

    if nargin < 4 || isempty(iter),    iter = 1;    end
    if nargin < 5 || isempty(maxiter), maxiter = 1; end

    tol = 1e-9;
    tviol = sum(gviol);
    nc    = numel(gviol);
    nihl  = sum(gviol > tol);          % number of violated constraints

    switch method
        case 1   % --- Deb: no scalar penalty ---
            toplam = NaN;

        case 2   % --- Static constant penalty (Penalty Method-1) ---
            K = 1000;
            toplam = -verim + K*nihl;

        case 3   % --- Violation sum, linear (Penalty Method-2) ---
            toplam = -verim + sum(gviol);

        case 4   % --- Sum of squared violations (Penalty Method-3) ---
            toplam = -verim + sum(gviol.^2);

        case 5   % --- Normalized violation-count penalty ---
            K = 1000;
            toplam = -verim + K*(nihl/nc);

        case 6   % --- Joines-Houck dynamic (Penalty Method-4) ---
            C = 1e-3; alpha = 1; beta = 1;      % same ramp as ceza_6 in the reference code
            toplam = -verim + ((C*iter)^alpha) * sum(gviol.^beta);

        case 7   % --- Liu et al. dynamic, linear (Penalty Method-5) ---
            mult = liu_mult(iter, maxiter);
            toplam = -verim + mult * sum(gviol);

        case 8   % --- Liu et al. dynamic, square ---
            mult = liu_mult(iter, maxiter);
            toplam = -verim + mult * sum(gviol.^2);

        otherwise
            error('penalty:badMethod', 'method must be 1..8 (given: %g).', method);
    end
end

function mult = liu_mult(gen, G)
%LIU_MULT  Liu et al. dynamic penalty multiplier (ceza_7/8 formula):
%   10^( (theta2-theta1)/(1+exp(20*(-gen+G/4)/G)) + theta1 ),  theta1=2, theta2=6
    theta1 = 2; theta2 = 6;
    mult = 10.^( (theta2-theta1)./(1+exp(20*(-gen+G/4)/G)) + theta1 );
end
