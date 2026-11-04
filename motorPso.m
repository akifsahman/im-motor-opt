function [fbestVerim, bestParams, ozet] = motorPso(maxrun, n, D, ST, method, max_fes)
%MOTORPSO  IM design optimization with PSO (single objective: maximize eff_h).
%   method 1..8 : penalty approach (1 = Deb's rule).
%   Output files are saved into the 'sonuclar/' folder.

    [dmin, dmax] = im_bounds();
    UB = ones(n,D).*dmax;  LB = ones(n,D).*dmin;
    c1 = 1.5; c2 = 2.5;
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
        trees = LB + (UB-LB).*rand(n,D);

        % --- initial evaluation ---
        Ver = zeros(1,n); Viol = zeros(1,n); Obj = zeros(1,n);
        for i=1:n
            [v,gv] = im_core(trees(i,:));
            Ver(i)=v; Viol(i)=sum(gv);
            if ~isDeb, Obj(i) = penalty(v,gv,method,1,maxiter); end
        end

        % --- initial gbest/pbest ---
        pbest = trees; pVer = Ver; pViol = Viol; pObj = Obj;
        if isDeb
            cbest=Viol(1); fbest=Ver(1); bestParams=trees(1,:);
            for i=2:n
                [cbest,fbest,bestParams]=deb(Viol(i),Ver(i),trees(i,:),cbest,fbest,bestParams);
            end
        else
            [obest,ib]=min(Obj); fbest=Ver(ib); cbest=Viol(ib); bestParams=trees(ib,:);
        end

        v = 0.1*trees; iter = 1;
        pHist = zeros(maxiter, D);
        while fes < max_fes
            w = 0.1 + rand*0.4;
            % --- velocity and position update ---
            for i=1:n
                v(i,:) = w*v(i,:) + c1*rand(1,D).*(pbest(i,:)-trees(i,:)) ...
                                  + c2*rand(1,D).*(bestParams-trees(i,:));
            end
            cand = trees + v;
            % boundary clipping
            cand = max(cand, LB); cand = min(cand, UB);

            for i=1:n
                [vv,gv] = im_core(cand(i,:));
                tv = sum(gv);
                if isDeb
                    % update individual (Deb)
                    [Viol(i),Ver(i),trees(i,:)] = deb(tv,vv,cand(i,:), Viol(i),Ver(i),trees(i,:));
                    % update pbest
                    [pViol(i),pVer(i),pbest(i,:)] = deb(Viol(i),Ver(i),trees(i,:), pViol(i),pVer(i),pbest(i,:));
                    % update gbest
                    [cbest,fbest,bestParams] = deb(pViol(i),pVer(i),pbest(i,:), cbest,fbest,bestParams);
                else
                    oc = penalty(vv,gv,method,iter,maxiter);
                    if oc < Obj(i)
                        trees(i,:)=cand(i,:); Ver(i)=vv; Viol(i)=tv; Obj(i)=oc;
                    end
                    if Obj(i) < pObj(i)
                        pbest(i,:)=trees(i,:); pVer(i)=Ver(i); pViol(i)=Viol(i); pObj(i)=Obj(i);
                    end
                    if pObj(i) < obest
                        obest=pObj(i); fbest=pVer(i); cbest=pViol(i); bestParams=pbest(i,:);
                    end
                end
            end
            fes = fes + n;
            if iter <= maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter = iter + 1;
            fprintf('PSO m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n', ...
                    method, rx, fes, iter, fbest, cbest);
        end

        % pad remaining columns with the last (converged) value
        if iter <= maxiter
            iterVerim(rx,iter:end)=fbest; iterCeza(rx,iter:end)=cbest;
            pHist(iter:end,:)=repmat(bestParams, maxiter-iter+1, 1);
        end
        if fbest > bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams;
        verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end

    fbestVerim = max(minVerim);
    ozet = save_results('PSO', method, minVerim, verimCeza, bestlerrun, iterVerim, iterCeza, paramHistBest);
end
