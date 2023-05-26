function [caldata xgaz ygaz]= calib(xyR,centerCorrect,dotinfo,calibType,rect,xraw,yraw)
if strcmp(calibType,'HV9') 
   xyRc = xyR(:,5);
   xyR = xyR-repmat(xyRc,1,size(xyR,2));
   ixC = [2,4,5,6,8];
elseif strcmp(calibType,'HV5')
% calibration dots top,left,center,right,bottom are the ones used for the basic calibration equation
    xyRc = xyR(:,3);
    xyR = xyR-repmat(xyRc,1,size(xyR,2)); 
    ixC = [1,2,3,4,5];
end
A   = [ones(5,1),xyR(1,ixC)',xyR(2,ixC)',xyR(1,ixC).^2',xyR(2,ixC).^2'];
% bx  = dotinfo.calibpos(ixC,1);
% by  = dotinfo.calibpos(ixC,2);
bx  = dotinfo.calibpos(ixC,1)-rect(3)/2;
by  = dotinfo.calibpos(ixC,2)-rect(4)/2;

ux  = A\bx;
uy  = A\by;

% xgazaux = ux'*[ones(1,size(xraw',2));xraw';yraw';xraw'.^2;yraw'.^2];
% ygazaux = uy'*[ones(1,size(yraw',2));xraw';yraw';xraw'.^2;yraw'.^2];
x_centerCorrect = centerCorrect(1);
y_centerCorrect = centerCorrect(2);

xgazaux = ux'*[ones(1,size(xraw'-xyRc(1),2));xraw'-xyRc(1);yraw'-xyRc(2);(xraw'-xyRc(1)).^2;(yraw'-xyRc(2)).^2];
ygazaux = uy'*[ones(1,size(yraw'-xyRc(2),2));xraw'-xyRc(1);yraw'-xyRc(2);(xraw'-xyRc(1)).^2;(yraw'-xyRc(2)).^2];

xyPaux  = ux'*[ones(1,size(xyR,2));xyR(1,:);xyR(2,:);xyR(1,:).^2;xyR(2,:).^2];
xyPaux  = [xyPaux;uy'*[ones(1,size(xyR,2));xyR(1,:);xyR(2,:);xyR(1,:).^2;xyR(2,:).^2]];

xyDriftaux  = ux'*[ones(1,size(x_centerCorrect'-xyRc(1),2));x_centerCorrect'-xyRc(1);y_centerCorrect'-xyRc(2);(x_centerCorrect'-xyRc(1)).^2;(y_centerCorrect'-xyRc(2)).^2];
xyDriftaux  = [xyDriftaux;uy'*[ones(1,size(y_centerCorrect'-xyRc(2),2));(x_centerCorrect'-xyRc(1));y_centerCorrect'-xyRc(2)';(x_centerCorrect'-xyRc(1)).^2;(y_centerCorrect'-xyRc(2)).^2]];
% xgaz    = xgazaux;
% ygaz    = ygazaux;

xgaz    = xgazaux+rect(3)/2;
ygaz    = (ygazaux+rect(4)/2);
xyP     = xyPaux+repmat(rect(3:4)'/2,1,size(xyPaux,2));
xyDrift = xyDriftaux + rect(3:4)'/2;

caldata.ux                  = ux;
caldata.uy                  = uy;
caldata.correctedDotPos     = xyP;
caldata.uncorrectedDotPoscentered   = xyR+repmat(xyRc,1,size(xyR,2));
caldata.uncorrectedDotPos   = xyR;
caldata.xyDrift             = xyDrift;
caldata.rawCenter           = xyRc;
caldata.rect                = rect;
caldata.calibType           = calibType;
caldata.valid_calibration   = 1;
% caldata.samples             = samples;
if strcmp(calibType,'HV9') 
    caldata.m                   = m;
    caldata.n                   = n;
elseif strcmp(calibType,'HV5') 
    caldata.m                   = NaN;
    caldata.n                   = NaN;
end