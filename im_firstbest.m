function [fbest,cbest,obest,bestParams] = im_firstbest(Ver,Viol,Obj,P,isDeb)
%IM_ILKBEST  Selects the initial best individual (Deb or scalar).
    n=size(P,1); obest=0;
    if isDeb
        cbest=Viol(1); fbest=Ver(1); bestParams=P(1,:);
        for i=2:n
            if im_better(Viol(i),Ver(i),0, cbest,fbest,0, true)
                cbest=Viol(i); fbest=Ver(i); bestParams=P(i,:);
            end
        end
    else
        [obest,ib]=min(Obj); fbest=Ver(ib); cbest=Viol(ib); bestParams=P(ib,:);
    end
end
