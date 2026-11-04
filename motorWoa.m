function [fbestVerim, bestParams, ozet] = motorWoa(maxrun, n, D, ST, method, max_fes)
%MOTORWOA  Whale Optimization Algorithm (WOA) + DE-mutant improvement, for the IM problem.
%   method 1 = Deb, 2..8 = penalty methods. Positions are updated
%   unconditionally by WOA; when stuck (flag>=2 and the candidate is not
%   better) a DE trial vector is applied (CR,F1,F2).
    [dmin,dmax]=im_bounds(); low=ones(n,D).*dmin; up=ones(n,D).*dmax;
    isDeb=(method==1); maxiter=floor(max_fes/n);
    CR=0.9; F1=0.1; F2=0.9;
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;

    for rx=1:maxrun
        rng('shuffle');
        P=low+(up-low).*rand(n,D);
        Ver=zeros(1,n); Viol=zeros(1,n); Obj=zeros(1,n);
        for i=1:n, [v,gv]=im_core(P(i,:)); Ver(i)=v; Viol(i)=sum(gv); if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end, end
        [fbest,cbest,obest,bestParams]=im_firstbest(Ver,Viol,Obj,P,isDeb);
        fes=n; iter=1; pHist=zeros(maxiter,D); flag=0;
        while fes<max_fes
            a=2-iter*(2/maxiter); a2=-1+iter*(-1/maxiter);
            for i=1:n
                A=2*a*rand-a; C=2*rand; b=1; l=(a2-1)*rand+1; p=rand;
                yeni=P(i,:);
                for j=1:D
                    if p<0.5
                        if abs(A)>=1
                            Xr=P(floor(n*rand+1),:); yeni(j)=Xr(j)-A*abs(C*Xr(j)-P(i,j));
                        else
                            yeni(j)=bestParams(j)-A*abs(C*bestParams(j)-P(i,j));
                        end
                    else
                        yeni(j)=abs(bestParams(j)-P(i,j))*exp(b*l)*cos(l*2*pi)+bestParams(j);
                    end
                end
                yeni=min(max(yeni,dmin),dmax);
                [vv,gv]=im_core(yeni); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end
                P(i,:)=yeni; Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end
                fes=fes+1;
                % update gbest
                [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,yeni,fbest,cbest,obest,bestParams,isDeb);
                % --- DE-mutant improvement (when stuck) ---
                iyi = im_better(tv,vv,oc, cbest,fbest,obest, isDeb);
                if flag>=2 && ~iyi
                    r1=fix(rand*n)+1; while r1==i, r1=fix(rand*n)+1; end
                    r2=fix(rand*n)+1; while r2==i||r2==r1, r2=fix(rand*n)+1; end
                    mut=P(i,:)+(bestParams-P(i,:))*F2+(P(r1,:)-P(r2,:))*F1;
                    tr=P(i,:); kr=fix(rand*D)+1;
                    for k=1:D, if rand<CR||k==kr, tr(k)=mut(k); end, end
                    tr=min(max(tr,dmin),dmax);
                    [vt,gt]=im_core(tr); tvt=sum(gt); oct=0; if ~isDeb, oct=penalty(vt,gt,method,iter,maxiter); end
                    fes=fes+1;
                    if im_better(tvt,vt,oct, Viol(i),Ver(i),Obj(i), isDeb)
                        P(i,:)=tr; Ver(i)=vt; Viol(i)=tvt; if ~isDeb, Obj(i)=oct; end
                        [fbest,cbest,obest,bestParams]=im_bestupd(vt,tvt,oct,tr,fbest,cbest,obest,bestParams,isDeb);
                    end
                    flag=0;
                end
                flag=flag+1;
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('WOA m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('WOA',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
