function [fbestVerim, bestParams, ozet] = motorCTSA(maxrun, n, D, ST, method, max_fes)
%MOTORCTSA  IM design optimization with the Tree-Seed Algorithm (TSA).
%   ST : search tendency parameter.
%   method 1..8 : penalty approach (1 = Deb's rule).

    [dmin, dmax] = im_bounds();
    isDeb   = (method == 1);
    maxiter = floor(max_fes/n);

    minVerim = zeros(1,maxrun);
    verimCeza = zeros(maxrun,2);
    bestlerrun = zeros(maxrun,D);
    % TSA's outer-iteration count is variable; maxiter is a safe upper bound.
    iterVerim = zeros(maxrun, maxiter);   % iteration-wise best efficiency
    iterCeza  = zeros(maxrun, maxiter);   % iteration-wise best violation
    paramHistBest = zeros(maxiter, D);    % parameter evolution of the BEST run
    bestVerimSoFar = -inf;

    for rx = 1:maxrun
        rng('shuffle');
        low  = round(n*0.1);   high = round(n*0.25);
        fes  = n;
        trees = dmin + (dmax-dmin).*rand(n,D);

        Ver = zeros(1,n); Viol = zeros(1,n); Obj = zeros(1,n);
        for i=1:n
            [v,gv] = im_core(trees(i,:));
            Ver(i)=v; Viol(i)=sum(gv);
            if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end
        end

        if isDeb
            cbest=Viol(1); fbest=Ver(1); bestParams=trees(1,:);
            for i=2:n
                [cbest,fbest,bestParams]=deb(Viol(i),Ver(i),trees(i,:),cbest,fbest,bestParams);
            end
        else
            [obest,ib]=min(Obj); fbest=Ver(ib); cbest=Viol(ib); bestParams=trees(ib,:);
        end

        iter = 1;
        pHist = zeros(maxiter, D);
        while fes < max_fes
            % best (current) tree - for direction in seed generation
            if isDeb
                bestIdx = 1;
                for i=2:n
                    [~,~,~] = deb(Viol(i),Ver(i),trees(i,:),Viol(bestIdx),Ver(bestIdx),trees(bestIdx,:));
                    if debWins(Viol(i),Ver(i),Viol(bestIdx),Ver(bestIdx)), bestIdx=i; end
                end
            else
                [~,bestIdx] = min(Obj);
            end
            bestPar = trees(bestIdx,:);

            for i=1:n
                ns = min(high, fix(low+(high-low)*rand)+1);
                for j=1:ns
                    r = fix(rand*n)+1; while r==i, r=fix(rand*n)+1; end
                    seed = trees(i,:);
                    for d=1:D
                        if rand < ST
                            seed(d) = trees(i,d) + (bestPar(d)-trees(r,d))*(2*rand-1);
                        else
                            seed(d) = trees(i,d) + (trees(i,d)-trees(r,d))*(2*rand-1);
                        end
                        if seed(d) > dmax(d) || seed(d) < dmin(d)
                            seed(d) = dmin(d) + (dmax(d)-dmin(d))*rand;   % resample on boundary violation
                        end
                    end
                    [vv,gv] = im_core(seed); tv = sum(gv);
                    if isDeb
                        [Viol(i),Ver(i),trees(i,:)] = deb(tv,vv,seed, Viol(i),Ver(i),trees(i,:));
                    else
                        oc = penalty(vv,gv,method,iter,maxiter);
                        if oc < Obj(i)
                            trees(i,:)=seed; Ver(i)=vv; Viol(i)=tv; Obj(i)=oc;
                        end
                    end
                end
                fes = fes + ns;
            end

            % update global best
            if isDeb
                for i=1:n
                    [cbest,fbest,bestParams]=deb(Viol(i),Ver(i),trees(i,:),cbest,fbest,bestParams);
                end
            else
                for i=1:n
                    if Obj(i) < obest
                        obest=Obj(i); fbest=Ver(i); cbest=Viol(i); bestParams=trees(i,:);
                    end
                end
            end

            if iter <= maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter = iter + 1;
            fprintf('TSA m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n', ...
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
    ozet = save_results('TSA', method, minVerim, verimCeza, bestlerrun, iterVerim, iterCeza, paramHistBest);
end

function tf = debWins(va, fa, vb, fb)
% does a beat b by Deb's rule?
    fea_a = (va<=0); fea_b = (vb<=0);
    if fea_a && fea_b, tf = (fa >= fb);
    elseif fea_a && ~fea_b, tf = true;
    elseif ~fea_a && fea_b, tf = false;
    else tf = (va <= vb); end
end
