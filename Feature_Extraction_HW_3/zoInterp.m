function [copies] = zoInterp(x,numInterp)
%Function that uses repmat and reshape functions to copy each value of x
%numInterp times

% mat=(repmat(x,numInterp));
% copies=reshape(mat(1,:),[1,numInterp]);
copies=reshape(repmat(x,numInterp,1),1,length(x)*numInterp);

end

