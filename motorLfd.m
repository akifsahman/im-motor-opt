function [fbestVerim, bestParams, ozet] = motorLfd(maxrun, n, D, ST, method, max_fes)
%MOTORLFD  Levy Flight Distribution (LFD) for the IM problem. method 1 = Deb, 2..8 = penalty.
%   Neighbor attraction (normalized distance threshold) + Levy flight toward gbest. Unconditional update.
    [dmin,dmax]=im_bounds(); isDeb=(method==1); maxiter=floor(max_fes/n);
    thr=0.25; rng_=dmax-dmin;
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;
    for rx=1:maxrun
        rng('shuffle');
        P=dmin+(dmax-dmin).*rand(n,D);
        Ver=zeros(1,n); Viol=zeros(1,n); Obj=zeros(1,n);
        for i=1:n, [v,gv]=im_core(P(i,:)); Ver(i)=v; Viol(i)=sum(gv); if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end, end
        [fbest,cbest,obest,bestParams]=im_firstbest(Ver,Viol,Obj,P,isDeb);
        fes=n; iter=1; pHist=zeros(maxiter,D);
        while fes<max_fes
            Sc=im_score(Ver,Viol,Obj,isDeb); smin=min(Sc); smax=max(Sc);
            Pn=(P-dmin)./rng_;
            for i=1:n
                S_i=zeros(1,D); nn=0;
                for j=1:n
                    if i~=j
                        dis=sqrt(mean((Pn(i,:)-Pn(j,:)).^2));
                        if dis<thr
                            nn=nn+1; vf=(-1)^(rand>0.5);
                            w=(Sc(j)-smin)/(smax-smin+eps)*0.9+0.1;
                            S_i=S_i+vf*w.*P(j,:);
                        end
                    end
                end
                if nn>0, S_i=S_i/nn; end
                Xr=P(floor(n*rand+1),:);
                yeni=bestParams+0.01*S_i+rand*5e-5*((bestParams+0.005*Xr)/2-P(i,:));
                yeni=lfd_levy(yeni,bestParams);
                yeni=min(max(yeni,dmin),dmax);
                [vv,gv]=im_core(yeni); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end
                P(i,:)=yeni; Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end; fes=fes+1;
                [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,yeni,fbest,cbest,obest,bestParams,isDeb);
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('LFD m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('LFD',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
function s=lfd_levy(s,DP)
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(size(s))*sigma; v=randn(size(s)); step=u./abs(v).^(1/beta);
    s=s+0.01*step.*(s-DP).*randn(size(s));
end
