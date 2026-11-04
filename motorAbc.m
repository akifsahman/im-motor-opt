function [fbestVerim, bestParams, ozet] = motorAbc(maxrun, n, D, ST, method, max_fes)
%MOTORABC  Artificial Bee Colony (ABC) for the IM problem: employed +
%   onlooker (roulette) + scout phases, with Deb/scalar acceptance.
%   method 1 = Deb, 2..8 = penalty methods.
    [dmin,dmax]=im_bounds(); isDeb=(method==1); maxiter=floor(max_fes/n);
    workers=floor(n/2); onlooker=workers; max_limit=workers*D;
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;

    for rx=1:maxrun
        rng('shuffle'); limit=zeros(1,workers);
        bees=dmin+(dmax-dmin).*rand(workers,D);
        Ver=zeros(1,workers); Viol=zeros(1,workers); Obj=zeros(1,workers);
        for i=1:workers, [v,gv]=im_core(bees(i,:)); Ver(i)=v; Viol(i)=sum(gv); if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end, end
        [fbest,cbest,obest,bestParams]=im_firstbest(Ver,Viol,Obj,bees,isDeb);
        fes=workers; iter=1; pHist=zeros(maxiter,D);
        while fes<max_fes
            % --- employed bees ---
            for i=1:workers
                [ns,vv,tv,oc]=abc_komsu(bees,i,workers,D,dmin,dmax,method,iter,maxiter,isDeb); fes=fes+1;
                if im_better(tv,vv,oc, Viol(i),Ver(i),Obj(i), isDeb)
                    bees(i,:)=ns; Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end; limit(i)=0;
                else, limit(i)=limit(i)+1; end
                [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,ns,fbest,cbest,obest,bestParams,isDeb);
            end
            % --- onlooker bees (roulette) ---
            Sc = secim_skoru(Ver,Viol,Obj,isDeb);
            fit = 1./(1 + Sc - min(Sc) + eps); ih = fit/sum(fit);
            t=0; i=1;
            while t<workers && fes<max_fes
                if rand<ih(i)
                    [ns,vv,tv,oc]=abc_komsu(bees,i,workers,D,dmin,dmax,method,iter,maxiter,isDeb); fes=fes+1;
                    if im_better(tv,vv,oc, Viol(i),Ver(i),Obj(i), isDeb)
                        bees(i,:)=ns; Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end; limit(i)=0;
                    else, limit(i)=limit(i)+1; end
                    [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,ns,fbest,cbest,obest,bestParams,isDeb);
                    t=t+1;
                end
                i=i+1; if i>onlooker, i=1; end
            end
            % --- scout bee ---
            [~,wi]=max(limit);
            if limit(wi)>max_limit
                bees(wi,:)=dmin+(dmax-dmin).*rand(1,D);
                [v,gv]=im_core(bees(wi,:)); Ver(wi)=v; Viol(wi)=sum(gv); if ~isDeb, Obj(wi)=penalty(v,gv,method,iter,maxiter); end
                limit(wi)=0; fes=fes+1;
                [fbest,cbest,obest,bestParams]=im_bestupd(Ver(wi),Viol(wi),Obj(wi),bees(wi,:),fbest,cbest,obest,bestParams,isDeb);
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('ABC m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('ABC',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
function [ns,vv,tv,oc]=abc_komsu(bees,i,workers,D,dmin,dmax,method,iter,maxiter,isDeb)
% neighbor (new candidate) generation for ABC
    rp=ceil(rand*D); rk=ceil(rand*workers); while rk==i, rk=ceil(rand*workers); end
    ns=bees(i,:); ns(rp)=bees(i,rp)+(bees(i,rp)-bees(rk,rp))*(2*rand-1);
    ns(rp)=min(max(ns(rp),dmin(rp)),dmax(rp));
    [vv,gv]=im_core(ns); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end
end
function Sc=secim_skoru(Ver,Viol,Obj,isDeb)
% selection score for roulette (lower = better)
    if isDeb, Sc=(Viol>1e-9)*1e6 + Viol*100 - Ver; else, Sc=Obj; end
end
