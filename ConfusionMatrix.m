function ConfMat = ConfusionMatrix(Y,Z)
%=======================================================================
% ConfMat = ConfusionMatrix(Y,Z)
%-----------------------------------------------------------------------
%
% ConfusionMatrix() returns an [2x2] confusion matrix
%
%  Default threshold of Z is 0.
%  "Y"     : Actual (binary) target {-1,1}, [nx1] vector, ... 
%  "Z"     : Predicted value of Y, real valued, [nx1] vector,... 
%
% (C) Copyright 2004.-, HyunJung (Helen) Shin (2004-07-04).
%
%=======================================================================

    ConfMat = size(2,2);
    ConfMat(1,1) = length(find(Y>0 & Z>0));
    ConfMat(1,2) = length(find(Y>0 & Z<0));

    ConfMat(2,1) = length(find(Y<0 & Z>0));
    ConfMat(2,2) = length(find(Y<0 & Z<0));

