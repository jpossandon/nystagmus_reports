function  [trial] = recalculate_eye_single(trial,caldata,pixxdeg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  [trial] = recalculate_eye(trial,caldata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eyes = fields(caldata);
% for ey = eyes'
     this_struct = trial;
     this_cal    = caldata;
    cfg.vfac = 4;
    cfg.vfaca = 6;
    cfg.mindur = 5;
    cfg.threshold = 0;
    cfg.thresholda = 0;
    cfg.resolution = pixxdeg;
    cfg.fs = 500;
    for tt = 1:length(this_struct)
        [this_struct(tt).samples.x,this_struct(tt).samples.y]  = correct_raw(this_struct(tt).samples.rawx',...
            this_struct(tt).samples.rawy',this_cal);
        
        this_struct(tt).samples.xvel = cfg.fs*conv(this_struct(tt).samples.x,[1 1 1 0 0 0 -1 -1 -1],'same')./18./pixxdeg;
        this_struct(tt).samples.yvel = cfg.fs*conv(this_struct(tt).samples.y,[1 1 1 0 0 0 -1 -1 -1],'same')./18./pixxdeg;
        this_struct(tt).samples.xacel = cfg.fs.^2*conv(this_struct(tt).samples.x,[1 1 0 -1 -2 -1 0 +1 +1],'same')./24./pixxdeg;
        this_struct(tt).samples.yacel = cfg.fs.^2*conv(this_struct(tt).samples.y,[1 1 0 -1 -2 -1 0 +1 +1],'same')./24./pixxdeg;
%         this_struct(tt).samples.xacel = conv(this_struct(tt).samples.xvel,[1 1 0 -1 -1],'same');
%         this_struct(tt).samples.yacel = conv(this_struct(tt).samples.yvel,[1 1 0 -1 -1],'same');
            [sac, radius, threshold(tt,:)] = saccade(cfg, [this_struct(tt).samples.x' this_struct(tt).samples.y']);
        
    end
    cfg.threshold = nanmean(threshold(:,1));
    cfg.thresholda = nanmean(threshold(:,2));
    for tt = 1:length(this_struct)
            [sac,~,threshold,sacstruct] = saccade(cfg, [this_struct(tt).samples.x' this_struct(tt).samples.y'],this_struct(tt).samples.time);
        
            sacstruct.threshold = threshold;
        this_struct(tt).resaccade        = sacstruct;
        this_struct(tt).samples.type     = nan(1,length(this_struct(tt).samples.x));
        for st = 1:size(sac,1)
            if sac(st,8) == 0
                this_struct(tt).samples.type(sac(st,1):sac(st,2)) = 2;  %sac
            else
                this_struct(tt).samples.type(sac(st,1):sac(st,2)) = 4;  %interrupted sac
            end
        end
        % potential fixations
         [fix,fixstruct] = fixations(this_struct(tt).samples.x,this_struct(tt).samples.y,isnan(this_struct(tt).samples.type),5,this_struct(tt).samples.time);
        this_struct(tt).refixation        = fixstruct;
         for st = 1:size(fix,1)
            if fix(st,10) == 0
                this_struct(tt).samples.type(fix(st,1):fix(st,2)) = 1;  %fix
            else
                this_struct(tt).samples.type(fix(st,1):fix(st,2)) = 3;  % interrupted fix
            end
         end
        % unclssified samples are 5, missing samples are NaNs
         this_struct(tt).samples.type(~(isnan(this_struct(tt).samples.x) | isnan(this_struct(tt).samples.y)) & isnan(this_struct(tt).samples.type)) = 5;
    end
    clear trial
    for tt = 1:length(this_struct)
        trial(tt) =  this_struct(tt);
    end
% end
