function [iterVerim,iterCeza,pHist] = im_queue(iterVerim,iterCeza,pHist,rx,iter,maxiter,fbest,cbest,bestParams)
%IM_KUYRUK  Pads the remainder of the iteration arrays with the last value
%   when a run finishes early (so all runs have equal-length histories).
    if iter<=maxiter
        iterVerim(rx,iter:end)=fbest; iterCeza(rx,iter:end)=cbest;
        pHist(iter:end,:)=repmat(bestParams, maxiter-iter+1, 1);
    end
end
