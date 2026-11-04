function [fbestVerim, bestParams, ozet] = motorGwo(maxrun, n, D, ST, method, max_fes)
%MOTORGWO  Grey Wolf Optimizer (GWO) for the IM problem. alpha/beta/delta are
%   updated PERSISTENTLY (cascaded, at every evaluation).
%   method 1 = Deb, 2..8 = penalty methods.
    [dmin,dmax]=im_bounds(); isDeb=(method==1); maxiter=floor(max_fes/n);
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;

    for rx=1:maxrun
        rng('shuffle');
        P=dmin+(dmax-dmin).*rand(n,D);
        % alpha/beta/delta leaders: [viol eff obj pos]
        Lv=[inf inf inf]; Lf=[-inf -inf -inf]; Lo=[inf inf inf]; Lp=zeros(3,D);
        fes=0; iter=1; pHist=zeros(maxiter,D);
        while fes<max_fes
            % --- evaluate + cascade-update the leaders ---
            for i=1:n
                P(i,:)=min(max(P(i,:),dmin),dmax);
                [vv,gv]=im_core(P(i,:)); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end
                fes=fes+1;
                if im_better(tv,vv,oc, Lv(1),Lf(1),Lo(1), isDeb)
                    Lv(3)=Lv(2);Lf(3)=Lf(2);Lo(3)=Lo(2);Lp(3,:)=Lp(2,:);
                    Lv(2)=Lv(1);Lf(2)=Lf(1);Lo(2)=Lo(1);Lp(2,:)=Lp(1,:);
                    Lv(1)=tv;Lf(1)=vv;Lo(1)=oc;Lp(1,:)=P(i,:);
                elseif im_better(tv,vv,oc, Lv(2),Lf(2),Lo(2), isDeb)
                    Lv(3)=Lv(2);Lf(3)=Lf(2);Lo(3)=Lo(2);Lp(3,:)=Lp(2,:);
                    Lv(2)=tv;Lf(2)=vv;Lo(2)=oc;Lp(2,:)=P(i,:);
                elseif im_better(tv,vv,oc, Lv(3),Lf(3),Lo(3), isDeb)
                    Lv(3)=tv;Lf(3)=vv;Lo(3)=oc;Lp(3,:)=P(i,:);
                end
            end
            fbest=Lf(1); cbest=Lv(1); bestParams=Lp(1,:);
            % --- position update (eqs. 3.3-3.7) ---
            a=2-iter*(2/maxiter);
            for i=1:n
                for j=1:D
                    A1=2*a*rand-a; X1=Lp(1,j)-A1*abs(2*rand*Lp(1,j)-P(i,j));
                    A2=2*a*rand-a; X2=Lp(2,j)-A2*abs(2*rand*Lp(2,j)-P(i,j));
                    A3=2*a*rand-a; X3=Lp(3,j)-A3*abs(2*rand*Lp(3,j)-P(i,j));
                    P(i,j)=(X1+X2+X3)/3;
                end
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('GWO m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('GWO',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
