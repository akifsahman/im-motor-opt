function [fbestVerim, bestParams, ozet] = motorFoa(maxrun, n, D, ST, method, max_fes)
%MOTORFOA  Red Fox Optimization (FOA) for the IM problem. method 1 = Deb, 2..8 = penalty.
%   Positions are updated unconditionally (gbest is kept separately).
    [dmin,dmax]=im_bounds(); isDeb=(method==1); maxiter=floor(max_fes/n);
    c1=0.18; c2=0.82;
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;
    for rx=1:maxrun
        rng('shuffle');
        X=dmin+(dmax-dmin).*rand(n,D);
        Ver=zeros(1,n); Viol=zeros(1,n); Obj=zeros(1,n);
        for i=1:n, [v,gv]=im_core(X(i,:)); Ver(i)=v; Viol(i)=sum(gv); if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end, end
        [fbest,cbest,obest,bestParams]=im_firstbest(Ver,Viol,Obj,X,isDeb);
        fes=n; iter=1; pHist=zeros(maxiter,D); MinT=inf;
        while fes<max_fes
            a=2*(1-iter/maxiter);
            for i=1:n
                r=rand; p=rand;
                if r>=0.5
                    Time=rand(1,D); DFR=0.5*(bestParams./Time).*Time;
                    tt=sum(Time)/D; Jump=0.5*9.81*(tt/2)^2;
                    if p>0.18, X(i,:)=DFR*Jump*c1; else, X(i,:)=DFR*Jump*c2; end
                    if MinT>tt, MinT=tt; end
                else
                    mt=MinT; if ~isfinite(mt), mt=1; end
                    X(i,:)=bestParams+randn(1,D)*(mt*a);
                end
                X(i,:)=min(max(X(i,:),dmin),dmax);
                [vv,gv]=im_core(X(i,:)); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end
                Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end; fes=fes+1;
                [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,X(i,:),fbest,cbest,obest,bestParams,isDeb);
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('FOA m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('FOA',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
