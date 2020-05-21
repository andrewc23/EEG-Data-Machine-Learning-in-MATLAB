function [features] = MovingWinFeats(x,fs,winLen,winDisp,featFn)
%Function that returns a vector of the values of the feature on the signal x in all possible windows
%Calculates line length in every window

NumWindows=floor((((length(x)/fs)-winLen)/winDisp)+1);
features=zeros(1,NumWindows);

%calculate line length in each window
starter=1;
for i=1:NumWindows
    %set up x
    in_between=x(starter:starter+winLen*fs-1);
    features(i)=featFn(in_between);
    starter=starter+winDisp*fs;
end

end

