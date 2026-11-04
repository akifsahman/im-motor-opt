function ozet = save_results(algo, method, minVerim, verimCeza, bestlerrun, iterVerim, iterCeza, iterParams)
%KAYDET  Saves per-run statistics + iteration-based convergence history.
%   Outputs are written as .mat files into the 'sonuclar/' folder.
%   File name pattern: ALGO_Mx_*.mat
%
%   WHAT IS SAVED:
%     *_verimCeza.mat     : maxrun x 2  [best efficiency, best violation] (per run)
%     *_bestlerrun.mat    : maxrun x D  (best design vector of each run)
%     *_iterVerim.mat     : maxrun x iter  (per-run, iteration-wise best efficiency)
%     *_iterCeza.mat      : maxrun x iter  (per-run, iteration-wise best violation)
%     *_ortIterVerim.mat  : 1 x iter  (mean over runs -> convergence curve, efficiency)
%     *_ortIterCeza.mat   : 1 x iter  (mean over runs -> convergence curve, violation)
%     *_ortalamalar.mat   : 2 x 5  row1=efficiency, row2=violation;
%                           columns = [mean, std, median, best, worst]
%     *_iterParams.mat    : iter x 14  iteration-wise design vector of the BEST run
%                           (parameter evolution over iterations)
%     *_ozet.mat          : struct (summary below)

    if ~exist('sonuclar','dir'), mkdir('sonuclar'); end
    etiket = sprintf('%s_M%d', algo, method);
    P = @(ad) fullfile('sonuclar', [etiket '_' ad '.mat']);

    % --- end-of-run summary statistics ---
    ozet.algo=algo; ozet.method=method;
    ozet.eniyiVerim = max(verimCeza(:,1));
    ozet.enkotuVerim= min(verimCeza(:,1));
    ozet.ortVerim   = mean(verimCeza(:,1));
    ozet.medyanVerim= median(verimCeza(:,1));
    ozet.sapmaVerim = std(verimCeza(:,1));
    ozet.ortIhlal   = mean(verimCeza(:,2));
    [~,ix] = max(verimCeza(:,1));
    ozet.eniyiParams = bestlerrun(ix,:);

    % summary table ('ortalamalar' format from the reference code)
    v = verimCeza(:,1); c = verimCeza(:,2);
    ortalamalar = [ mean(v) std(v) median(v) max(v) min(v);   % efficiency: max=best
                    mean(c) std(c) median(c) min(c) max(c) ];  % violation:  min=best

    save(P('verimCeza'),   'verimCeza');
    save(P('bestlerrun'),  'bestlerrun');
    save(P('ortalamalar'), 'ortalamalar');
    save(P('ozet'),        'ozet');

    % --- iteration-based convergence history (for plotting) ---
    if nargin >= 7 && ~isempty(iterVerim)
        ortIterVerim = mean(iterVerim, 1);   % convergence curve (efficiency)
        ortIterCeza  = mean(iterCeza,  1);   % convergence curve (violation)
        save(P('iterVerim'),    'iterVerim');
        save(P('iterCeza'),     'iterCeza');
        save(P('ortIterVerim'), 'ortIterVerim');
        save(P('ortIterCeza'),  'ortIterCeza');
    end

    % --- parameter evolution of the best run ---
    if nargin >= 8 && ~isempty(iterParams)
        save(P('iterParams'), 'iterParams');
    end
end
