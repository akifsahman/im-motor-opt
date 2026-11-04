function [fbestVerim, bestParams, ozet] = motorAaa(maxrun, n, D, ST, method, max_fes)
%MOTORAAA  Artificial Algae Algorithm (AAA) for the IM problem. method 1 = Deb, 2..8 = penalty.
%   Energy (chlorophyll) + friction surface + size; neighbor tournament selection; evolution + adaptation.
    [dmin,dmax]=im_bounds(); isDeb=(method==1); maxiter=floor(max_fes/n);
    K=2; le=0.3; Ap=0.5;
    minVerim=zeros(1,maxrun); verimCeza=zeros(maxrun,2); bestlerrun=zeros(maxrun,D);
    iterVerim=zeros(maxrun,maxiter); iterCeza=zeros(maxrun,maxiter);
    paramHistBest=zeros(maxiter,D); bestVerimSoFar=-inf;
    for rx=1:maxrun
        rng('shuffle');
        ALG=dmin+(dmax-dmin).*rand(n,D); STARVE=zeros(1,n); BIG=ones(1,n);
        Ver=zeros(1,n); Viol=zeros(1,n); Obj=zeros(1,n);
        for i=1:n, [v,gv]=im_core(ALG(i,:)); Ver(i)=v; Viol(i)=sum(gv); if ~isDeb, Obj(i)=penalty(v,gv,method,1,maxiter); end, end
        [fbest,cbest,obest,bestParams]=im_firstbest(Ver,Viol,Obj,ALG,isDeb);
        fes=n; iter=1; pHist=zeros(maxiter,D);
        while fes<max_fes
            Sc=im_score(Ver,Viol,Obj,isDeb);
            KLORO=aaa_great(BIG); SURF=aaa_fric(BIG);
            for i=1:n
                istarve=0;
                while KLORO(i)>=0 && fes<max_fes
                    Nb=aaa_turnuva(Sc); while Nb==i, Nb=aaa_turnuva(Sc); end
                    pr=randperm(D,3); NEW=ALG(i,:);
                    NEW(pr(1))=NEW(pr(1))+(ALG(Nb,pr(1))-NEW(pr(1)))*(K-SURF(i))*(rand-0.5)*2;
                    NEW(pr(2))=NEW(pr(2))+(ALG(Nb,pr(2))-NEW(pr(2)))*(K-SURF(i))*cos(rand*2*pi);
                    NEW(pr(3))=NEW(pr(3))+(ALG(Nb,pr(3))-NEW(pr(3)))*(K-SURF(i))*sin(rand*2*pi);
                    NEW=min(max(NEW,dmin),dmax);
                    [vv,gv]=im_core(NEW); tv=sum(gv); oc=0; if ~isDeb, oc=penalty(vv,gv,method,iter,maxiter); end; fes=fes+1;
                    KLORO(i)=KLORO(i)-le/2;
                    if im_better(tv,vv,oc, Viol(i),Ver(i),Obj(i), isDeb)
                        ALG(i,:)=NEW; Ver(i)=vv; Viol(i)=tv; if ~isDeb, Obj(i)=oc; end; Sc(i)=im_score(vv,tv,oc,isDeb); istarve=1;
                    else, KLORO(i)=KLORO(i)-le/2; end
                    [fbest,cbest,obest,bestParams]=im_bestupd(vv,tv,oc,NEW,fbest,cbest,obest,bestParams,isDeb);
                end
                if istarve==0, STARVE(i)=STARVE(i)+1; end
            end
            BIG=aaa_buyukluk(BIG,Sc);
            m1=fix(rand*D)+1; [~,imax]=max(BIG); [~,imin]=min(BIG); ALG(imin,m1)=ALG(imax,m1);
            [~,istrv]=max(STARVE);
            if rand<Ap
                ALG(istrv,:)=min(max(ALG(istrv,:)+(ALG(imax,:)-ALG(istrv,:)).*rand(1,D),dmin),dmax);
            end
            if iter<=maxiter, iterVerim(rx,iter)=fbest; iterCeza(rx,iter)=cbest; pHist(iter,:)=bestParams; end
            iter=iter+1;
            fprintf('AAA m%d Run=%d Fes=%d Iter=%d eff=%.4f viol=%.3f\n',method,rx,fes,iter,fbest,cbest);
        end
        [iterVerim,iterCeza,pHist]=im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams);
        if fbest>bestVerimSoFar, bestVerimSoFar=fbest; paramHistBest=pHist; end
        minVerim(rx)=fbest; bestlerrun(rx,:)=bestParams; verimCeza(rx,1)=fbest; verimCeza(rx,2)=cbest;
    end
    fbestVerim=max(minVerim);
    ozet=save_results('AAA',method,minVerim,verimCeza,bestlerrun,iterVerim,iterCeza,paramHistBest);
end
function y=aaa_great(B)   % energy (chlorophyll) from rank of size
    s2=numel(B); [~,sira]=sort(B); y=zeros(1,s2); for i=1:s2, y(sira(i))=i^2; end; y=aaa_norm(y);
end
function y=aaa_fric(B)    % friction surface (hemispherical) from size
    s2=numel(B); y=zeros(1,s2); for i=1:s2, r=((B(i)*3)/(4*pi))^(1/3); y(i)=2*pi*r^2; end; y=aaa_norm(y);
end
function B=aaa_buyukluk(B,Sc)   % size (bigness) growth update
    O=1-aaa_norm(Sc);
    for i=1:numel(B), fKs=abs(B(i)/2); M=O(i)/(fKs+O(i)); B(i)=B(i)+M*B(i); end
end
function c=aaa_turnuva(src)   % tournament selection (2-way, lower score wins)
    nb=randperm(numel(src),2); if src(nb(1))<src(nb(2)), c=nb(1); else, c=nb(2); end
end
function y=aaa_norm(x)   % min-max normalization to [0,1]
    mn=min(x); mx=max(x); if mx-mn<eps, y=zeros(size(x)); else, y=(x-mn)/(mx-mn); end
end
