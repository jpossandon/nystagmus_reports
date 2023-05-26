function [fix,fixstruct] = fixations(x,y,notfix,minsampleln,times)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fix = fixations(x,y,notfix,minsampleln)
%
%  OUTPUT:
%
%  fix(1:num,1)   onset of fixation
%  fix(1:num,2)   end of fixation
%  fix(1:num,3)   peak velocity of fixation (vpeak)
%  fix(1:num,4)   median horizontal pos     
%  fix(1:num,5)   median vertical pos

%  fix(1:num,6)   horizontal component     (dx)
%  fix(1:num,7)   vertical component       (dy)
%  fix(1:num,8)   horizontal amplitude     (dX)
%  fix(1:num,9)   vertical amplitude       (dY)
%  fix(1:num,10)   start after or ends before missing data or data collection
% 19.07.19 JPO
% Hamburg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ispfix = ~(isnan(x) | isnan(y)) & notfix;
fix    = findsegments(ispfix,minsampleln);
x      = [x',y']; % to reuse the code below from miccrosaccade Engbert's
for s=1:size(fix,1)
    % onset and offset
    a = fix(s,1); 
    b = fix(s,2); 
    i = fix(s,1):fix(s,2);
    % fixcade peak velocity (vpeak)
    %vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
    vpeak = nan; % not yet calculated
    fix(s,3) = vpeak;
    % fixcade vector (dx,dy)
    dx = x(b,1)-x(a,1); 
    dy = x(b,2)-x(a,2); 
    fix(s,4) = median(x(i,1));
    fix(s,5) = median(x(i,2));
    fix(s,6) = dx;
    fix(s,7) = dy;
    % fixcade amplitude (dX,dY)
    
    [minx, ix1] = min(x(i,1));
    [maxx, ix2] = max(x(i,1));
    [miny, iy1] = min(x(i,2));
    [maxy, iy2] = max(x(i,2));
    dX = sign(ix2-ix1)*(maxx-minx);
    dY = sign(iy2-iy1)*(maxy-miny);
    fix(s,8:9) = [dX dY];
    if fix(s,1)>5 && fix(s,2)<size(x,1)-5
        if any(any(isnan(x(fix(s,1)-5:fix(s,1)-1,:)))) || any(any(isnan(x(fix(s,2)+1:fix(s,2)+5,:))))
            fix(s,10) = 1;
        else
            fix(s,10) = 0;
        end
    else
        fix(s,10) = 1;
    end
end

if nargout>1
    if isempty(fix)
        fixstruct = struct('start',[],'end',[],'sx',[],'sy',[],'ex',[],'ey',[],'speed',[],'x',[],'y',[],'dur',[],'amp',[],'vec',[],'interrupt',[]);
    else
    fixstruct.start  = times(fix(:,1));
    fixstruct.end    = times(fix(:,2));
    fixstruct.sx     = x(fix(:,1),1)';
    fixstruct.sy     = x(fix(:,1),2)';
    fixstruct.ex     = x(fix(:,2),1)';
    fixstruct.ey     = x(fix(:,2),2)';
    fixstruct.speed  = fix(:,3)';
    fixstruct.x      = fix(:,4)';
    fixstruct.y      = fix(:,5)';
    fixstruct.dur    = fixstruct.end-fixstruct.start;
    fixstruct.amp    = sqrt(fix(:,8).^2+fix(:,9).^2)';
    fixstruct.vec    = sqrt(fix(:,6).^2+fix(:,7).^2)';
    fixstruct.interrupt = fix(:,10)';
    end
end