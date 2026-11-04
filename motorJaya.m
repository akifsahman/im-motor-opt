function [fbestVerim, bestParams, ozet] = motorJaya(maxrun, n, D, ST, method, max_fes)
%MOTORJAYA  Jaya algorithm (greedy) for the IM problem. method 1 = Deb, 2..8 = penalty.
%   xnew = x + r1*(Best-|x|) - r2*(worst-|x|); accept if the candidate is better.
    [dmin,dmax]=im_bounds(); isDeb=(method==1); maxiter=floor(max_fes/n);
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;
    for rx=1:maxrun
        rng('shuffle');
        x=dmin+(dmax-dmin).*rand(n,D);
        Ver=zeros(1,n); Viol=zeros(1,n); Obj=zeros(1,n);
        for i=1:n, [v,gv]=im_core(x(i,:)); Ver(i)=v; Viol(i)=sum(gv); if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end, end
        [fbest,cbest,obest,bestParams]=im_firstbest(Ver,Viol,Obj,x,isDeb);
        fes=n; iter=1; pHist=zeros(maxiter,D);
        while fes<max_fes
            Sc=im_score(Ver,Viol,Obj,isDeb); [~,bi]=min(Sc); [~,wi]=max(Sc);
            Best=x(bi,:); worst=x(wi,:);
            for i=1:n
                yeni=x(i,:)+rand(1,D).*(Best-abs(x(i,:))) - rand(1,D).*(worst-abs(x(i,:)));
                yeni=min(max(yeni,dmin),dmax);
                [vv,gv]=im_core(yeni); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end
                fes=fes+1;
                if im_better(tv,vv,oc, Viol(i),Ver(i),Obj(i), isDeb)
                    x(i,:)=yeni; Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end
                end
                [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,yeni,fbest,cbest,obest,bestParams,isDeb);
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('JAYA m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('JAYA',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
