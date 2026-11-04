% RUN_ALL  IM design optimization: 10 algorithms x 8 penalty approaches.
%   Single objective: MAXIMIZE eff_h (efficiency).
%   Requirement: all .m files in this folder + the hardened lookup functions
%   (power_factor, BmHm, Kst_Kf, Kst_Alfa_i, esson_constant, conductor_diameter,
%   efficiency_lookup) must be on the MATLAB path.
clear; clc;

pop     = 40;          % population size (n)
D       = 16;          % number of variables (16: Lamda_k, KD_k added)
ST      = 0.2;         % TSA seeking tendency
max_fes = 100000;      % max number of objective-function evaluations
maxrun  = 30;          % number of independent runs (for statistics)

% 10 algorithms. Common signature: motorXXX(maxrun,pop,D,ST,method,max_fes)
algolar = {'TSA','PSO','DSA','WOA','GWO','JAYA','FOA','LFD','ABC','AAA'};
motorfn = {@motorCTSA,@motorPso,@motorDsa,@motorWoa,@motorGwo, ...
           @motorJaya,@motorFoa,@motorLfd,@motorAbc,@motorAaa};
% Penalty approaches:
%  1 Deb | 2 Static constant | 3 Violation sum | 4 Squared violation
%  5 Normalized count | 6 Joines-Houck | 7 Liu linear | 8 Liu square
yontemler = 1:8;

sonucTablo = {};
for a = 1:numel(algolar)
    for m = yontemler
        [~,~,oz] = motorfn{a}(maxrun,pop,D,ST,m,max_fes);
        sonucTablo(end+1,:) = {algolar{a}, m, oz.eniyiVerim, oz.ortVerim, ...
                               oz.sapmaVerim, oz.ortIhlal};
        fprintf('=== %s  Method %d : best=%.4f  mean=%.4f  meanViol=%.4f ===\n', ...
                algolar{a}, m, oz.eniyiVerim, oz.ortVerim, oz.ortIhlal);
    end
end

if ~exist('sonuclar','dir'), mkdir('sonuclar'); end
save(fullfile('sonuclar','sonucTablo.mat'), 'sonucTablo');
disp('--- ALL RUNS DONE (10 algorithms x 8 methods). Summary: sonuclar/sonucTablo.mat ---');
