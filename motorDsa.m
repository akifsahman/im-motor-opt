function [fbestVerim, bestParams, ozet] = motorDsa(maxrun, n, D, ST, method, max_fes)
%MOTORDSA  IM design optimization with the Differential Search Algorithm (DSA).
%   Uses Elitist-1 (E1-DSA) direction generation.
%   method 1..8 : penalty approach (1 = Deb's rule).

    [dmin, dmax] = im_bounds();
    up = ones(n,D).*dmax;  low = ones(n,D).*dmin;
    isDeb   = (method == 1);
    maxiter = floor(max_fes/n);

    minVerim = zeros(1,maxrun);
    verimCeza = zeros(maxrun,2);
    bestlerrun = zeros(maxrun,D);
    iterVerim = zeros(maxrun, maxiter);   % iteration-wise best efficiency
    iterCeza  = zeros(maxrun, maxiter);   % iteration-wise best violation
    paramHistBest = zeros(maxiter, D);    % parameter evolution of the BEST run
    bestVerimSoFar = -inf;

    for rx = 1:maxrun
        rng('shuffle');
        fes = n;
        super = low + (up-low).*rand(n,D);

        Ver = zeros(1,n); Viol = zeros(1,n); Obj = zeros(1,n);
        for i=1:n
            [v,gv] = im_core(super(i,:));
            Ver(i)=v; Viol(i)=sum(gv);
            if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end
        end

        if isDeb
            cbest=Viol(1); fbest=Ver(1); bestParams=super(1,:);
            for i=2:n
                [cbest,fbest,bestParams]=deb(Viol(i),Ver(i),super(i,:),cbest,fbest,bestParams);
            end
        else
            [obest,ib]=min(Obj); fbest=Ver(ib); cbest=Viol(ib); bestParams=super(ib,:);
        end

        iter = 1;
        pHist = zeros(maxiter, D);
        while fes < max_fes
            p1 = 0.3*rand;  p2 = 0.3*rand;
            % --- direction: Elitist-1 (toward a random top-best) ---
            [~,jind] = sort(-Ver);                 % efficiency high to low
            ibest = jind(ceil(rand*n));
            direction = repmat(super(ibest,:),[n 1]);
            % --- active-individuals map ---
            map = aktif_harita(n,D,p1,p2);
            R = 1./(-0.5*log(rand));   % Gamma(1,0.5)=Exponential(mean=0.5); no toolbox needed
            stopover = super + (R.*map).*(direction - super);
            % boundary clipping
            stopover = max(stopover, low); stopover = min(stopover, up);

            for i=1:n
                [vv,gv] = im_core(stopover(i,:));
                tv = sum(gv);
                if isDeb
                    [Viol(i),Ver(i),super(i,:)] = deb(tv,vv,stopover(i,:), Viol(i),Ver(i),super(i,:));
                    [cbest,fbest,bestParams]    = deb(Viol(i),Ver(i),super(i,:), cbest,fbest,bestParams);
                else
                    oc = penalty(vv,gv,method,iter,maxiter);
                    if oc < Obj(i)
                        super(i,:)=stopover(i,:); Ver(i)=vv; Viol(i)=tv; Obj(i)=oc;
                    end
                    if Obj(i) < obest
                        obest=Obj(i); fbest=Ver(i); cbest=Viol(i); bestParams=super(i,:);
                    end
                end
            end
            fes = fes + n;
            if iter <= maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter = iter + 1;
            fprintf('DSA m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n', ...
                    method, rx, fes, iter, fbest, cbest);
        end

        if iter <= maxiter
            iterVerim(rx,iter:end)=fbest; iterCeza(rx,iter:end)=cbest;
            pHist(iter:end,:)=repmat(bestParams, maxiter-iter+1, 1);
        end
        if fbest > bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams;
        verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end

    fbestVerim = max(minVerim);
    ozet = save_results('DSA', method, minVerim, verimCeza, bestlerrun, iterVerim, iterCeza, paramHistBest);
end

function map = aktif_harita(n, clan, p1, p2)
    map = zeros(n, clan);
    if rand < rand
        if rand < p1
            for i=1:n, map(i,:) = rand(1,clan) < rand; end          % Random-mut #1
        else
            for i=1:n, map(i,randi(clan)) = 1; end                  % Differential-mut
        end
    else
        for i=1:n, map(i,randi(clan,1,ceil(p2*clan))) = 1; end      % Random-mut #2
    end
end
