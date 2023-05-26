function [caldata] = calib_nystag_test(trial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [caldata xgaz ygaz] = calib_nystag_test(trial);
% Specific calibration for nystagmus test used in the ODP
% Horizontal and Vertical test start showing 5 points, 5 s each in a cross
% grid in the following order:
% 960  ,540 center 5s
% 1296 ,540 right 5s
% 624  ,540 left 5s
% 960  ,204 down 5s
% 960  ,876 up 5s
% The first trial is done with the right eye and the second with the left
% eye
%
% José Ossandon (jose.ossandon@uni-hamburg.de)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scr_siz          = [0 0 1920 1080];
eye_order        = {'right','left'};

% which samples used to caibrate after the apeparance of each dot
samples_to_calib = [500 3500]; 
% velocity percentile to use
velperc          = 20;
% lower limit of sample per point to calibrate in percentage 
sample_limit     = .2;%

% Reoder to top,left,center,right,bottom,center (as used in calib.m)
reorder = [5,3,1,2,4,6];

for tr = 1:length(eye_order)
    
    dotinfo.calibpos        = trial(tr).disp_scr.value(reorder,:);
    dotinfo.times           = trial(tr).disp_scr.time(reorder);
    
    % get raw data values used to calibrate
    % first allvelocities
    allvel       = sqrt(trial(tr).(eye_order{tr}).samples.rawxvel.^2+trial(tr).(eye_order{tr}).samples.rawyvel.^2);
    vellim       = prctile(allvel,velperc);
    allvel(allvel>vellim) = NaN;
    xyR = [];
    calib_valid = 1;
    %check
%      figure,plot(trial(tr).(eye_order{tr}).samples.time,trial(tr).(eye_order{tr}).samples.rawx,'k'), hold on
    
    for calpos = 1:length(dotinfo.times)
        idxtimes        = trial(tr).(eye_order{tr}).samples.time>dotinfo.times(calpos)+samples_to_calib(1) &...
                            trial(tr).(eye_order{tr}).samples.time<dotinfo.times(calpos)+samples_to_calib(2);
        xyR(:,calpos)   = [nanmedian(trial(tr).(eye_order{tr}).samples.rawx(idxtimes)) nanmedian(trial(tr).(eye_order{tr}).samples.rawy(idxtimes))]';
        
        if sum(~isnan(trial(tr).(eye_order{tr}).samples.rawx(idxtimes)))<sum(idxtimes)*sample_limit || ...
                sum(~isnan(trial(tr).(eye_order{tr}).samples.rawy(idxtimes)))<sum(idxtimes)*sample_limit
                calib_valid = 0;
        end
        % check
%          line([dotinfo.times(calpos)+samples_to_calib(1) dotinfo.times(calpos)+samples_to_calib(2)],[xyR(1,calpos) xyR(1,calpos)])
    end
    centerCorrect           = xyR(:,end); % last (6th) drift correction
    xraw                    = trial(tr).(eye_order{tr}).samples.rawx';
    yraw                    = trial(tr).(eye_order{tr}).samples.rawy';
    [caldata.(eye_order{tr}) xgaz ygaz] = calib(xyR,centerCorrect,dotinfo,'HV5',scr_siz,xraw,yraw);
    caldata.(eye_order{tr}).valid_calibration = calib_valid;
end
    