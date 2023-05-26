function [sac, radius, threshold, sacstruct] = saccade(cfg, data,times)%(x,vel,VFAC,MINDUR,THRESHOLD)
%-------------------------------------------------------------------
%
%  FUNCTION microsacc.m
%  Detection of monocular candidates for microsaccades;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 2.1, 03 OCT 05)
%   29/07/08 JPO modified for saccade detection with a fixed treshold
%               - THREShOLD :    0 , it uses the relative velocity threshold
%                               >0 , it uses THRESHOLD as a fixed one
%   19/07/19 JPO modified for identifying interupted saccades (blinks end
%   start and also provide resutls as a structure
%-------------------------------------------------------------------
%
%  INPUT:
%
%  x(:,1:2)         position vector
%  vel(:,1:2)       velocity vector
%  VFAC             relative velocity threshold
%  MINDUR           minimal saccade duration
%
%  OUTPUT:
%
%  sac(1:num,1)   onset of saccade
%  sac(1:num,2)   end of saccade
%  sac(1:num,3)   peak velocity of saccade (vpeak)
%  sac(1:num,4)   horizontal component     (dx)
%  sac(1:num,5)   vertical component       (dy)
%  sac(1:num,6)   horizontal amplitude     (dX)
%  sac(1:num,7)   vertical amplitude       (dY)
%  sac(1:num,8)   start after or ends before missing data or data collection
%---------------------------------------------------------------------

if nargout>3 && nargin<3
    error('cannot provide sac struct wihtout time info')
end

x           = data;
normx       = data./cfg.resolution;
[vel acel]  = vecvel(normx, cfg.fs, 3);  
acel(find(vel(:,1)>=1000000 | vel(:,1)<=-1000000),:) = 0;
vel(find(vel(:,1)>=1000000 | vel(:,1)<=-1000000),:) = 0;
VFAC        = cfg.vfac;

MINDUR      = cfg.mindur;
THRESHOLD   = cfg.threshold;
if isfield(cfg,'vfaca') & isfield(cfg,'thresholda')
    VFACA        = cfg.vfaca;
    THRESHOLDA   = cfg.thresholda;
end

% compute threshold
msdx = sqrt( nanmedian(vel(:,1).^2) - (nanmedian(vel(:,1)))^2 );
msdy = sqrt( nanmedian(vel(:,2).^2) - (nanmedian(vel(:,2)))^2 );
msdxa = sqrt( nanmedian(acel(:,1).^2) - (nanmedian(acel(:,1)))^2 );
msdya = sqrt( nanmedian(acel(:,2).^2) - (nanmedian(acel(:,2)))^2 );
if msdx<realmin
    msdx = sqrt( nanmean(vel(:,1).^2) - (nanmean(vel(:,1)))^2 );
%     if msdx<realmin
%         error('msdx<realmin in microsacc.m');
%     end
end
if msdy<realmin
    msdy = sqrt( nanmean(vel(:,2).^2) - (nanmean(vel(:,2)))^2 );
%     if msdy<realmin
%         error('msdy<realmin in microsacc.m');
%     end
end
if msdx>realmin & msdy>realmin
    if THRESHOLD == 0
        radiusx = VFAC*msdx;
        radiusy = VFAC*msdy;
        radius = [radiusx radiusy];
    else
        radiusx = THRESHOLD;
        radiusy = THRESHOLD;
        radius = [THRESHOLD THRESHOLD];
    end

    if isfield(cfg,'vfaca') & isfield(cfg,'thresholda')
        if THRESHOLDA == 0
            radiusxa = VFACA*msdxa;
            radiusya = VFACA*msdya;
            radiusa = [radiusxa radiusya];
        else
            radiusxa = THRESHOLDA;
            radiusya = THRESHOLDA;
            radiusa = [THRESHOLDA THRESHOLDA];
        end
    thresholda = min(radiusa);    
    end
    threshold = min(radius);


    % compute test criterion: ellipse equation
    test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2;
    indx = find(test>1);
    if isfield(cfg,'vfaca') & isfield(cfg,'thresholda')
        testa = (acel(:,1)/radiusxa).^2 + (acel(:,2)/radiusya).^2;
        indxa = find(testa>1);
    end

    % determine saccades
    N = length(indx); 
    sac = [];
    nsac = 0;
    dur = 1;
    a = 1;
    k = 1;
    while k<N
        if indx(k+1)-indx(k)==1
            dur = dur + 1;
        else
            if dur>=MINDUR
                b = k; 
                if isfield(cfg,'vfaca') & isfield(cfg,'thresholda')
                    if ~isempty(find(indxa>indx(a) & indxa<indx(b))) 
                        nsac = nsac + 1;
                        sac(nsac,:) = [indx(a) indx(b)];
                    end
                else
                    nsac = nsac + 1;
                   sac(nsac,:) = [indx(a) indx(b)];
                end
            end
            a = k+1;
            dur = 1;
        end
        k = k + 1;
    end

    % check for minimum duration
    if dur>=MINDUR
        nsac = nsac + 1;
        b = k;
        sac(nsac,:) = [indx(a) indx(b)];
    end

    % compute peak velocity, horizonal and vertical components
    for s=1:nsac
        % onset and offset
        a = sac(s,1); 
        b = sac(s,2); 
        % saccade peak velocity (vpeak)
        vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
        sac(s,3) = vpeak;
        % saccade vector (dx,dy)
        dx = x(b,1)-x(a,1); 
        dy = x(b,2)-x(a,2); 
        sac(s,4) = dx;
        sac(s,5) = dy;
        % saccade amplitude (dX,dY)
        i = sac(s,1):sac(s,2);
        [minx, ix1] = min(x(i,1));
        [maxx, ix2] = max(x(i,1));
        [miny, iy1] = min(x(i,2));
        [maxy, iy2] = max(x(i,2));
        dX = sign(ix2-ix1)*(maxx-minx);
        dY = sign(iy2-iy1)*(maxy-miny);
        sac(s,6:7) = [dX dY];
        if sac(s,1)>5 & sac(s,2)<size(x,1)-5
            if any(any(isnan(x(sac(s,1)-5:sac(s,1)-1,:)))) || any(any(isnan(x(sac(s,2)+1:sac(s,2)+5,:))))
                sac(s,8) = 1;
            else
                sac(s,8) = 0;
            end
        else
            sac(s,8) = 1;
        end
    end
    if isfield(cfg,'vfaca') & isfield(cfg,'thresholda')
        threshold = [threshold thresholda];
    end
else
    sac = [];
    if isfield(cfg,'vfaca') & isfield(cfg,'thresholda')
        threshold = [NaN NaN];
    else
        threshold = NaN;
    end
    radius = NaN;
end
if nargout>3
    if isempty(sac)
        sacstruct = struct('start',[],'end',[],'sx',[],'sy',[],'ex',[],'ey',[],'speed',[],'dur',[],'amp',[],'vec',[],'interrupt',[]);
    else
        sacstruct.start  = times(sac(:,1));
        sacstruct.end    = times(sac(:,2));
        sacstruct.sx     = x(sac(:,1),1)';
        sacstruct.sy     = x(sac(:,1),2)';
        sacstruct.ex     = x(sac(:,2),1)';
        sacstruct.ey     = x(sac(:,2),2)';
        sacstruct.speed  = sac(:,3)';
        sacstruct.dur    = sacstruct.end-sacstruct.start;
        sacstruct.amp    = sqrt(sac(:,6).^2+sac(:,7).^2)';
        sacstruct.vec    = sqrt(sac(:,4).^2+sac(:,5).^2)';
        sacstruct.interrupt = sac(:,8)';
    end
end
