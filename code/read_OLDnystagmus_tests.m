%function read_nystagmus_tests(patient_id,PATH)
clear
patient_id = inputdlg('Patient ID:','Reports');
%%
% needs totrial.m, calib.m, correctraw.m
% relace perctile
% check that horizontal, vertical, and center files are there
% and open the filed
% calibrate each eye from horizontal and vertical test
% if is not possible to calibrate one eye, replace from the other
% test, if no eye abailable abort
%PATH = '/Users/jossando/trabajo/ODPdesktop/results';
%outPATH = '/Users/jossando/trabajo/ODPdesktop/results';
PATH = 'C:\Users\Experimenter\Desktop\nystagmography\nystagmus_tests\old';
outPATH = 'C:\Users\Experimenter\Desktop\nystagmus_reports\reports';

test_types              = {'horizontal','vertical','center'};
trials_eyes.horizontal  = {{'right'},{'left'},{'right','left'}};
trials_times.horizontal = [30,30,60];
trials_eyes.vertical    = {{'right'},{'left'},{'right','left'}};
trials_times.vertical   = [30,30,60];
trials_eyes.center      = {{'right','left'}};
trials_times.center     = [240];
eye_plot_length         = 30; % duration of  each eye position plot in ms
scr_size                = [1920 1080];
scr_distance            = 61;
scr_wdth                = 53;
PPD  = scr_size(1)/(2*atan(scr_wdth/2/scr_distance)*180/pi); % TODO acutal calcualtion

%PPD                     = 40.6;
srate                   = 500; % in hez
psd_winsiz              = 1000; %in samples
rem_segm_within_s_tgt   = .5; % segments that are within 500 ms of a target are not inculded in spectogram
scrsiz_deg         = scr_size/PPD;
%%
% here we open the different test files and calibrate them
for tt = 1:length(test_types)
    filename = fullfile(PATH,test_types{tt},'results',patient_id{1},[patient_id{1} '.edf']);
    if exist(filename)
        [trial.(test_types{tt})]      = totrial(filename,{'raw','gaze'});
        if ~strcmp(test_types{tt},'center') % center test does not have 'calibration' trials
            [caldata.(test_types{tt})]    = calib_nystag_test(trial.(test_types{tt}));
        end
        flag_exist(tt) = 1; % data exist
    else
        fprintf('\n%s TEST DOES NOT SEEM TO EXIST IN PATH\n%s\n',test_types{tt})
        flag_exist(tt) = 0; % data exist
    end
end


%%
% calibrate each eye from horizontal and vertical test
% if is not possible to calibrate one eye, replace from the other
% test, if no eye abailable abort

for tt = 1:length(test_types)
    for tr = 1:length(trials_eyes.(test_types{tt}))
        this_eyes = trials_eyes.(test_types{tt}){tr};
        for te = 1:length(this_eyes)
            thisCal = [];
            if ismember(test_types{tt},{'horizontal','vertical'})
                % if the calibration of this eye in the same test is not
                % valid use the other one
                if caldata.(test_types{tt}).(this_eyes{te}).valid_calibration
                    thisCal = test_types{tt};
                elseif caldata.(test_types{find(~ismember({'horizontal','vertical'},test_types{tt}))}).(this_eyes{te}).valid_calibration
                    thisCal = test_types{find(~ismember({'horizontal','vertical'},test_types{tt}))};
                end
            elseif ismember(test_types{tt},{'center'})
                % for the center test we use the calibrations from the
                % other test, if both available for a given eye, we use the
                % one with the smallest error for the 6th center point
                if caldata.horizontal.(this_eyes{te}).valid_calibration & caldata.vertical.(this_eyes{te}).valid_calibration  % if both valid
                    if sqrt(sum((caldata.horizontal.(this_eyes{te}).xyDrift-caldata.horizontal.(this_eyes{te}).rect(3:4)'/2).^2))<...
                            sqrt(sum((caldata.vertical.(this_eyes{te}).xyDrift-caldata.vertical.(this_eyes{te}).rect(3:4)'/2).^2))
                        thisCal = 'horizontal';
                    else
                        thisCal = 'vertical';
                    end
                elseif caldata.horizontal.(this_eyes{te}).valid_calibration
                    thisCal = 'horizontal';
                elseif caldata.vertical.(this_eyes{te}).valid_calibration
                    thisCal = 'vertical';
                end
            end
            if ~isempty(thisCal)
                [trial.(test_types{tt})(tr).(this_eyes{te})] = recalculate_eye_single(trial.(test_types{tt})(tr).(this_eyes{te}),...
                    caldata.(thisCal).(this_eyes{te}),PPD);
                
                
                %                 [trial.(test_types{tt})(tr).right.samples.x trial.(test_types{tt})(tr).right.samples.y] = ...
                %                     correct_raw(trial.(test_types{tt})(tr).right.samples.rawx',...
                %                                 trial.(test_types{tt})(tr).right.samples.rawy',caldata.(thisCal).(this_eyes{te}));
                fprintf('%s test %d trial %s eye calibration done from %s test values\n',test_types{tt},tr,this_eyes{te},thisCal)
            else
                fprintf('%s test %d trial %s eye calibration not posible\n',test_types{tt},tr,this_eyes{te})
                
            end
        end
    end
end

%%
% TODO
% correct raw with nans, done
% velocity
% spectra
% axes
%%
xplot_pos    = [.001 .035 .86 .8875];
xplot_width  = [.04 .8 .12 .07];
yplot_pos    = [.04 .15 .29 .4 .54 .65 .77 .88];

% yplot_pos_L  = [.11 .37 .62 .87];
eye_planes   = {'x','y'};
t_segments   = 30;
eyecolors    = [0.89412      0.10196       0.1098;
    0.21569      0.49412      0.72157];

eyecolorsVel    = [ 0.98431      0.70588      0.68235
    0.70196      0.80392       0.8902];
eyesLabels   = {'right','left'};
plotsPerPage = 8;
whichpage=1;
axesFontSize = 6;

% NOFF steps
NOFFwinlenngth   = 4;
NOFFwinstep      = .5;
NOFFxcorrvelperc = 20;
NOFFxmax         = 1;
NOFFxvelmax      = 6;
allMaxNOFF = struct('valuepopt',[],'valueNOFF',[],'test_type',[],'eye',[],'time',[],'eyecond',[]);
iaMN = 1;
for tt = 1:length(test_types)
    fh = figure;
    fh.Units     = 'centimeters';
    fh.Position  = [0 0 21 29.7];
    
    vert_progess = 1; %where we are on the page
%     tplots = 0;
    % every test trial
    for tr = 1:length(trial.(test_types{tt}))
        %trials_eye_times.vertical   = {{30},{30},{60}};
        total_tr_length = trials_times.(test_types{tt})(tr);
        nsegm = total_tr_length/eye_plot_length;
        
        % data is divided in 30 s segments
        % plots
        
        t_start = 0;
        for nplots = 1:nsegm
            if  vert_progess==plotsPerPage+1
                fh = figure;
                fh.Units     = 'centimeters';
                fh.Position  = [0 0 21 29.7];
                vert_progess = 1;
            end
            
            %
            % test type, patient id,other info
            if vert_progess==1
                subplot('Position',[.02 .98 .9 .02])
                axis off
                if ~flag_exist(tt)
                    text(0,1,sprintf('%s Test: %s data does not exist',upper(test_types{tt}),patient_id{1}),'VerticalAlignment','top','FontSize',12)
                else
                    text(0,1,sprintf('%s Test: %s',upper(test_types{tt}),patient_id{1}),'VerticalAlignment','top','FontSize',12)
                end
            end
            for eye_p = 1:2
                % Main plot eyes position separated by horizontal and
                % vertical movement
                subplot('Position',[xplot_pos(2) 1-(yplot_pos(vert_progess)+.07) xplot_width(2) .07])
                hold on
                if nplots==1 &eye_p==1
                    text(0,scrsiz_deg(eye_p)/2,sprintf('%s eye open',cell2mat(upper(join(trials_eyes.(test_types{tt}){tr},' & ')))),'VerticalAlignment','bottom','FontSize',8)
                end
                % plot position of tartgets in this period
                thisTgts = trial.(test_types{tt})(tr).disp_scr;
                [C,IA,thisTgts.IC] = unique(thisTgts.value,'rows');
                thisTgts.time  = thisTgts.time/1000;
                thisTgts.origvalue = thisTgts.value;
                thisTgts.value = [thisTgts.value(:,1)/PPD-scrsiz_deg(1)/2 thisTgts.value(:,2)/PPD-scrsiz_deg(2)/2];
                idx_Tgts = find(thisTgts.time>t_start & thisTgts.time<t_start+eye_plot_length);
                for tTgt = 1:length(idx_Tgts)
                    if tTgt < length(idx_Tgts)
                        lasttime = thisTgts.time(idx_Tgts(tTgt+1));
                    else
                        lasttime = 2*thisTgts.time(idx_Tgts(tTgt))-thisTgts.time(idx_Tgts(tTgt-1));
                    end
                    line([thisTgts.time(idx_Tgts(tTgt)) lasttime],...
                        [thisTgts.value(idx_Tgts(tTgt),eye_p) thisTgts.value(idx_Tgts(tTgt),eye_p)],'Color',[.5 .5 .5])
                    text(thisTgts.time(idx_Tgts(tTgt)),thisTgts.value(idx_Tgts(tTgt),eye_p),num2str(thisTgts.IC(idx_Tgts(tTgt))),...
                        'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6)
                end
                % and plot one or two eyes depending of the trial
                for eye_to_use = 1:length(trials_eyes.(test_types{tt}){tr})
                    thisEye     = trials_eyes.(test_types{tt}){tr}{eye_to_use};
                    thisData    = trial.(test_types{tt})(tr).(thisEye).samples;
                    thisData.time = thisData.time/1000;
                    idx_toplot  = thisData.time>t_start & thisData.time<=t_start+eye_plot_length;
                    plot(thisData.time(idx_toplot),thisData.(eye_planes{eye_p})(idx_toplot)/PPD-scrsiz_deg(eye_p)/2,'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:))
                    axis([t_start t_start+eye_plot_length -scrsiz_deg(eye_p)/2  scrsiz_deg(eye_p)/2])
                    % calculate NOFF if horizontal
                    if strcmp(eye_planes{eye_p},'x')
                        xwholeperiod    = thisData.(eye_planes{eye_p})(idx_toplot)/PPD-scrsiz_deg(eye_p)/2;
                        xvelwholeperiod = thisData.([eye_planes{eye_p} 'vel'])(idx_toplot);
                        twholeperiod    = thisData.time(idx_toplot);
                        thiswin_popt    = [];
                        thiswin_segmst  = [];
                        for tNOFF = twholeperiod(1):NOFFwinstep:twholeperiod(end)-NOFFwinlenngth
                            thisNOFFwinIdx = twholeperiod>tNOFF & twholeperiod<tNOFF+NOFFwinlenngth;
                            thisNOFFx      = xwholeperiod(thisNOFFwinIdx);
                            thisNOFFxvel   = xvelwholeperiod(thisNOFFwinIdx);
                            velthrforx     = prctile(abs(thisNOFFxvel),NOFFxcorrvelperc);
                            thisxcorr      = nanmedian(thisNOFFx(abs(thisNOFFxvel)<velthrforx));
                            thiswin_popt   = [thiswin_popt,sum(abs(thisNOFFx-thisxcorr)<NOFFxmax & abs(thisNOFFxvel)<NOFFxvelmax)./length(thisNOFFx)];
                            thiswin_segmst = [thiswin_segmst,tNOFF];
                        end
                        [maxNOFF ixmaxNOFF] = max(thiswin_popt);
                        if ~isempty(maxNOFF)
                            allMaxNOFF(iaMN).valuepopt = maxNOFF;
                            allMaxNOFF(iaMN).valueNOFF  = log(maxNOFF./(1-maxNOFF));
                            allMaxNOFF(iaMN).test_type  = test_types{tt};
                            allMaxNOFF(iaMN).eye        = trials_eyes.(test_types{tt}){tr}{eye_to_use};
                            allMaxNOFF(iaMN).eyecond    = cell2mat(join(trials_eyes.(test_types{tt}){tr}));
                            allMaxNOFF(iaMN).time       = thiswin_segmst(ixmaxNOFF);
                            iaMN = iaMN+1;
                            line([thiswin_segmst(ixmaxNOFF) thiswin_segmst(ixmaxNOFF)],[-scrsiz_deg(eye_p)/2  scrsiz_deg(eye_p)/2],'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:),'LineWidth',2)
                            line([thiswin_segmst(ixmaxNOFF)+NOFFwinlenngth thiswin_segmst(ixmaxNOFF)+NOFFwinlenngth],[-scrsiz_deg(eye_p)/2  scrsiz_deg(eye_p)/2],'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:),'LineWidth',2)
    %                         
                                if strcmp(trials_eyes.(test_types{tt}){tr}{eye_to_use},'right')
                                    text(t_start+eye_plot_length,scrsiz_deg(eye_p)/2,sprintf('RE popt: %2.2f  NOFF: %2.1f',maxNOFF,log(maxNOFF./(1-maxNOFF))),'FontSize',6,'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:));
                                elseif strcmp(trials_eyes.(test_types{tt}){tr}{eye_to_use},'left') 
                                    text(t_start+eye_plot_length,scrsiz_deg(eye_p)/2,sprintf('\n\nLE popt: %2.2f  NOFF: %2.1f',maxNOFF,log(maxNOFF./(1-maxNOFF))),'FontSize',6,'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:));
                                end
                         end
                    end
                    
                end
                %set axis/labels and other symbols
                if eye_p == 1
                    text(t_start+.2 ,floor(scrsiz_deg(eye_p)/2),'right',...
                        'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','top')
                    text(t_start+.2 ,ceil(-scrsiz_deg(eye_p)/2),'left',...
                        'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','bottom')
                    set(gca,'YTick',[ceil(-scrsiz_deg(eye_p)/2) 0 floor(scrsiz_deg(eye_p)/2)],...
                        'XTick',t_start:5:t_start+30,'XTickLabel',{},'Fontsize',axesFontSize)
                    ylh = ylabel('Horiz. Pos. (\circ)','Fontsize',axesFontSize);
                    ylh.Units = 'Normalized';
                    ylh.Position(1) = -.02;
                else
                    text(t_start+.2 ,floor(scrsiz_deg(eye_p)/2),'up',...
                        'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','top')
                    text(t_start+.2 ,ceil(-scrsiz_deg(eye_p)/2),'down',...
                        'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','bottom')
                    set(gca,'YTick',[ceil(-scrsiz_deg(eye_p)/2) 0 floor(scrsiz_deg(eye_p)/2)],...
                        'XTick',t_start:5:t_start+30,'XTickLabel',{},'Fontsize',axesFontSize)
                    ylh = ylabel('Vert. Pos. (\circ)','Fontsize',axesFontSize);
                    ylh.Units = 'Normalized';
                    ylh.Position(1) = -.02;
                end
                
                
                %VELOCITY PLOT
                % Main plot eyes position separated by horizontal and
                % vertical movement
                subplot('Position',[xplot_pos(2) 1-(yplot_pos(vert_progess)+.1) xplot_width(2) .03])
                hold on
                
                
                % and plot one or two eyes depending of the trial
                for eye_to_use = 1:length(trials_eyes.(test_types{tt}){tr})
                    thisEye     = trials_eyes.(test_types{tt}){tr}{eye_to_use};
                    thisData    = trial.(test_types{tt})(tr).(thisEye).samples;
                    thisData.time = thisData.time/1000;
                    idx_toplot  = thisData.time>t_start & thisData.time<=t_start+eye_plot_length;
                    plot(thisData.time(idx_toplot),thisData.([eye_planes{eye_p} 'vel'])(idx_toplot),'Color',eyecolorsVel(find(ismember(eyesLabels,thisEye )),:))
                    
                    axis([t_start t_start+eye_plot_length -150 150])
                end
                %set axis/labels and other symbols
                if eye_p == 1
                    %                     text(t_start+.2 ,floor(scrsiz_deg(eye_p)/2),'right',...
                    %                         'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','top')
                    %                     text(t_start+.2 ,ceil(-scrsiz_deg(eye_p)/2),'left',...
                    %                         'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','bottom')
                    
                    set(gca,'YTick',[-150 0 150],'YTickLabel',{'-','0','+'},...
                        'XTick',t_start:5:t_start+30,'XTickLabel',{},'Fontsize',axesFontSize)
                    ylh = ylabel('Vel.','Fontsize',axesFontSize);
                    ylh.Units = 'Normalized';
                    ylh.Position(1) = -.02;
                else
                    %                     text(t_start+.2 ,floor(scrsiz_deg(eye_p)/2),'up',...
                    %                         'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','top')
                    %                     text(t_start+.2 ,ceil(-scrsiz_deg(eye_p)/2),'down',...
                    %                         'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','bottom')
                    text(t_start+eye_plot_length+.35,  -220, 's','Fontsize',axesFontSize)
                    
                    set(gca,'YTick',[-150 0 150],'YTickLabel',{'-','0','+'},...
                        'XTick',t_start:5:t_start+30,'Fontsize',axesFontSize)
                    ylh = ylabel('Vel. ','Fontsize',axesFontSize);
                    ylh.Units = 'Normalized';
                    ylh.Position(1) = -.02;
                end
                
                % SPECTROGRAM
                
                
                
%                 % Main plot eyes position separated by horizontal and
%                 % vertical movement
%                 if eye_p == 1
%                     subplot('Position',[xplot_pos(4) 1-(yplot_pos(vert_progess)+.055) xplot_width(4) xplot_width(4)])
%                 else
%                     subplot('Position',[xplot_pos(4) 1-(yplot_pos(vert_progess)+.105) xplot_width(4) xplot_width(4)])
%                 end
%                 axis square
%                 hold on
%                 
%                 % and plot one or two eyes depending of the trial
%                 for eye_to_use = 1:length(trials_eyes.(test_types{tt}){tr})
%                     thisEye     = trials_eyes.(test_types{tt}){tr}{eye_to_use};
%                     thisData    = trial.(test_types{tt})(tr).(thisEye).samples;
%                     thisData.time = thisData.time/1000;
%                     idx_toplot  = thisData.time>t_start & thisData.time<=t_start+eye_plot_length;
%                     thisDatatime =  thisData.time(idx_toplot);
%                     thisData     = thisData.(eye_planes{eye_p})(idx_toplot)/PPD;
%                     
%                     % cut data in 1 second segments with 50% overlap
%                     topwelch = [];
%                     for segments = 1:psd_winsiz/4:length(thisData)-psd_winsiz
%                         
%                         if ~any((thisTgts.time>thisDatatime(segments)&thisTgts.time<thisDatatime(segments+psd_winsiz-1))|...
%                                 (thisTgts.time+rem_segm_within_s_tgt>thisDatatime(segments)&thisTgts.time+rem_segm_within_s_tgt<thisDatatime(segments+psd_winsiz-1)))
%                             if ~any(isnan(thisData(segments:segments+psd_winsiz-1)))
%                                 topwelch = [topwelch,thisData(segments:segments+psd_winsiz-1)'-mean(thisData(segments:segments+psd_winsiz-1)')];
%                             end
%                         end
%                     end
%                     if ~isempty(topwelch)
%                         [px,f]          = periodogram(topwelch,hamming(psd_winsiz),psd_winsiz,srate);
%                         plot(f(1:20),mean(px(1:20,:),2),'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:))
%                         if eye_to_use==1
%                             text(11,-.15*max(mean(px(1:20,:),2)),' Hz','Fontsize',axesFontSize)
%                         end
%                     end
%                     %                     axis([t_start t_start+eye_plot_length -scrsiz_deg(eye_p)/2  scrsiz_deg(eye_p)/2])
%                 end
%                 %                 %set axis/labels and other symbols
%                 %                 if eye_p == 1
%                 set(gca,...%'YTick',[],...
%                     'XTick',0:1:10,'XTickLabel',{'0','','2','','4','','6','','8','','10'},'Fontsize',axesFontSize)
%                 grid on
%                 ylh = ylabel('Power','Fontsize',axesFontSize);
                %                     xlh = xlabel('Freq. (Hz)','Fontsize',axesFontSize);
                %                     ylh.Units = 'Normalized';
                %                     ylh.Position(1) = -.02;
                %                  end
                %
                
                
                
                % set 2D plot andleft box plot showing where the patient have to look
                if eye_p == 1
                    
                    % Corresponding 2Dplot
                    subplot('Position',[xplot_pos(3) 1-(yplot_pos(vert_progess)+.105+xplot_width(3)/1.7778/2) xplot_width(3) xplot_width(3)/1.7778])
                    hold on
                    for eye_to_use = 1:length(trials_eyes.(test_types{tt}){tr})
                        thisEye     = trials_eyes.(test_types{tt}){tr}{eye_to_use};
                        thisData    = trial.(test_types{tt})(tr).(thisEye).samples;
                        thisData.time = thisData.time/1000;
                        idx_toplot  = find(thisData.time>t_start & thisData.time<=t_start+eye_plot_length);
                        plot(thisData.x(idx_toplot(1:5:end)),thisData.y(idx_toplot(1:5:end)),'MarkerSize',3,'Color',eyecolors(find(ismember(eyesLabels,thisEye )),:))
                        axis([0 scr_size(1) 0 scr_size(2)])
                    end
                    set(gca,'XTick',[],'YTick',[],'box','on')
                    
                    %left box
                    uniqTgt = unique(thisTgts.origvalue,'rows');
                    %                     subplot('Position',[xplot_pos(1) 1-(yplot_pos(vert_progess)+.12) xplot_width(1) xplot_width(1)/1.7778])
                    %                     hold on
                    %                     axis([-scrsiz_deg(1)/2 scrsiz_deg(1)/2 -scrsiz_deg(2)/2 scrsiz_deg(2)/2])
                    line([[0;scr_size(1)],[scr_size(1)/2;scr_size(1)/2]],[[scr_size(2)/2;scr_size(2)/2],[0;scr_size(2)]],'LineStyle',':','Color',[.5 .5 .5])
                    text(uniqTgt(:,1),uniqTgt(:,2),cellfun(@num2str,num2cell(1:size(uniqTgt,1)),'UniformOutput',false),...
                        'FontSize',8,'HorizontalAlignment','Center','VerticalAlignment','middle')
                    %                     set(gca,'XTick',[],'YTick',[],'box','on')
                end
                
                
                vert_progess = vert_progess+1; %where we are on the page
            end
            t_start       = t_start+30;
            if vert_progess == plotsPerPage+1 || (tr==length(trial.(test_types{tt})) & nplots==nsegm)%(tplots+nplots==plotsPerPage || (tt==3 & nplots==nsegm) || (tt<3 && tr==3 && nplots==nsegm))% || (tt==3 && (nplots==nsegm/2 || nplots==nsegm))
                fh.Position = [0 0 21 29.7];
                set(gcf, 'PaperPositionMode', 'auto')
                if whichpage==1
                    export_fig(fullfile(outPATH,patient_id{1}),'-pdf','-transparent')
                    whichpage = whichpage+1; 
                else
                    export_fig(fullfile(outPATH,patient_id{1}),'-pdf','-transparent','-append')
                end
                close(fh)
            end
%             if nplots==nsegm
%              tplots = tplots+nplots;
%             end
        end
        
    end
end
%         close all

%%

fh = figure;
fh.Units     = 'centimeters';
fh.Position  = [0 0 21 29.7];

vert_progess = 1; %where we are on the page
t_start = 0;
for eye_p = 1:2
    subplot('Position',[xplot_pos(2) 1-(yplot_pos(vert_progess)+.07) xplot_width(2) .07])
    axis([t_start t_start+eye_plot_length -scrsiz_deg(eye_p)/2  scrsiz_deg(eye_p)/2])
    text(t_start+2.5,0,'a','FontWeight','bold','FontSize',20)
    if eye_p == 1
        text(t_start+.2 ,floor(scrsiz_deg(eye_p)/2),'right',...
            'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','top')
        text(t_start+.2 ,ceil(-scrsiz_deg(eye_p)/2),'left',...
            'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','bottom')
        set(gca,'YTick',[ceil(-scrsiz_deg(eye_p)/2) 0 floor(scrsiz_deg(eye_p)/2)],...
            'XTick',t_start:5:t_start+30,'XTickLabel',{},'Fontsize',axesFontSize)
        ylh = ylabel('Horiz. Pos. (\circ)','Fontsize',axesFontSize);
        ylh.Units = 'Normalized';
        ylh.Position(1) = -.02;
    else
        text(t_start+.2 ,floor(scrsiz_deg(eye_p)/2),'up',...
            'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','top')
        text(t_start+.2 ,ceil(-scrsiz_deg(eye_p)/2),'down',...
            'Fontsize',axesFontSize,'Color',[.4 .4 .4],'HorizontalAlignment','left','VerticalAlignment','bottom')
        set(gca,'YTick',[ceil(-scrsiz_deg(eye_p)/2) 0 floor(scrsiz_deg(eye_p)/2)],...
            'XTick',t_start:5:t_start+30,'XTickLabel',{},'Fontsize',axesFontSize)
        ylh = ylabel('Vert. Pos. (\circ)','Fontsize',axesFontSize);
        ylh.Units = 'Normalized';
        ylh.Position(1) = -.02;
    end
    subplot('Position',[xplot_pos(2) 1-(yplot_pos(vert_progess)+.1) xplot_width(2) .03])
    axis([t_start t_start+eye_plot_length -150 150])
     text(t_start+2.5,0,'b','FontWeight','bold','FontSize',20)
    if eye_p == 1
        set(gca,'YTick',[-150 0 150],'YTickLabel',{'-','0','+'},...
            'XTick',t_start:5:t_start+30,'XTickLabel',{},'Fontsize',axesFontSize)
        ylh = ylabel('Vel.','Fontsize',axesFontSize);
        ylh.Units = 'Normalized';
        ylh.Position(1) = -.02;
    else
        text(t_start+eye_plot_length+.3,  -220, 's','Fontsize',axesFontSize)
        
        set(gca,'YTick',[-150 0 150],'YTickLabel',{'-','0','+'},...
            'XTick',t_start:5:t_start+30,'Fontsize',axesFontSize)
        ylh = ylabel('Vel. ','Fontsize',axesFontSize);
        ylh.Units = 'Normalized';
        ylh.Position(1) = -.02;
    end
    if eye_p == 1
        subplot('Position',[xplot_pos(4) 1-(yplot_pos(vert_progess)+.055) xplot_width(4) xplot_width(4)])
    else
        subplot('Position',[xplot_pos(4) 1-(yplot_pos(vert_progess)+.105) xplot_width(4) xplot_width(4)])
    end
    xlim([0 10])
      text(5,0.5,'c','FontWeight','bold','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle')
    axis square
     set(gca,'YTick',[],...
                    'XTick',0:1:10,'XTickLabel',{'0','','2','','4','','6','','8','','10'},'Fontsize',axesFontSize)
                grid on
                ylh = ylabel('Power','Fontsize',axesFontSize);
    


% set 2D plot andleft box plot showing where the patient have to look
    if eye_p == 1

        % Corresponding 2Dplot
        subplot('Position',[xplot_pos(3) 1-(yplot_pos(vert_progess)+.105+xplot_width(3)/1.7778/2) xplot_width(3) xplot_width(3)/1.7778])
        hold on
        
        set(gca,'XTick',[],'YTick',[],'box','on')

        %left box
        uniqTgt = unique(thisTgts.origvalue,'rows');
        %                     subplot('Position',[xplot_pos(1) 1-(yplot_pos(vert_progess)+.12) xplot_width(1) xplot_width(1)/1.7778])
        %                     hold on
        %                     axis([-scrsiz_deg(1)/2 scrsiz_deg(1)/2 -scrsiz_deg(2)/2 scrsiz_deg(2)/2])
        line([[0;scr_size(1)],[scr_size(1)/2;scr_size(1)/2]],[[scr_size(2)/2;scr_size(2)/2],[0;scr_size(2)]],'LineStyle',':','Color',[.5 .5 .5])
        text(uniqTgt(:,1),uniqTgt(:,2),cellfun(@num2str,num2cell(1:size(uniqTgt,1)),'UniformOutput',false),...
            'FontSize',8,'HorizontalAlignment','Center','VerticalAlignment','middle')
          text(50,scr_size(2),'d','FontWeight','bold','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','top')
   
        %                     set(gca,'XTick',[],'YTick',[],'box','on')
    end
        vert_progess = vert_progess+1;
end

fontSizeText = 12;
 subplot('Position',[xplot_pos(2) .1 .9 .6])
 text(0.01,1,sprintf('a - Horizontal/vertical eye position across 30 second, in visual degrees.\n''0'' is the midlle of the screen and horizontal lines mark the target position.'),'FontSize',fontSizeText,'VerticalAlignment','top')
 text(0.01,.9,sprintf('b - Horizontal/vertical eye velocity in visual degrees/second. Positive and negative values \nindicate movement to the right/upward or left/downward, respectively.'),'FontSize',fontSizeText,'VerticalAlignment','top')
 text(0.01,.8,sprintf('c - Spectrogram of eye positions shown in ''a''. It show the power at different frequencies of eye\ndisplacement, in Hertz. The large peak around 1Hz is, in general, due to large position transients \nand does not usually represent the nystagmus frequency'),'FontSize',fontSizeText,'VerticalAlignment','top')
 text(0.01,.65,sprintf('d - Position of the eye in the screen. Numbers indicate target locations.'),'FontSize',fontSizeText,'VerticalAlignment','top')
 text(0.01,.6,sprintf('Calibration of each eye is done with respect to the 5 targets presented in the first \ntwo trials of HORIZONTAL an VERTICAL test. When is not possible to calibrate \n(due to bad data or for the CENTER test), the calibration of the other test is used.\n\nTHIS CALIBRATION PROCEDURE IS AUTOMATIC AND PRONE TO ERROR, USE \nTHE VALUES (POSITION, VELOCITY and NOFF) ONLY AS A ROUGH APROXIMATION.'),'FontSize',fontSizeText,'VerticalAlignment','top')
 
 % best NOFF
 valuepopts  = [allMaxNOFF.valuepopt];
 ixright     = find(ismember({allMaxNOFF.eyecond},'right'));
 [maxval ix] = max(valuepopts(ixright));
  text(0.01,.3,sprintf('    RIGHT EYE MONOCULAR VIEWING (only 30 seconds)\n    Best foveation fraction (popt) and NOFF: %1.3f / %1.3f', maxval,log(maxval/(1-maxval))),'FontSize',fontSizeText,'VerticalAlignment','top')
 ixleft     = find(ismember({allMaxNOFF.eyecond},'left'));
 [maxval ix] = max(valuepopts(ixleft));
   text(0.01,.2,sprintf('    LEFT EYE MONOCULAR VIEWING (only 30 seconds)\n    Best foveation fraction (popt) and NOFF: %1.3f / %1.3f', maxval,log(maxval/(1-maxval))),'FontSize',fontSizeText,'VerticalAlignment','top')
 ixboth      = ismember({allMaxNOFF.eyecond},'right left');
 ixright     = ismember({allMaxNOFF.eye},'right');
 ixleft      = ismember({allMaxNOFF.eye},'left');
 [maxval ix] = max(valuepopts(ixboth&ixright));
 text(0.01,.1,sprintf('    RIGHT EYE BINOCULAR VIEWING\n    Best foveation fraction (popt) and NOFF: %1.3f / %1.3f', maxval,log(maxval/(1-maxval))),'FontSize',fontSizeText,'VerticalAlignment','top')
 [maxval ix] = max(valuepopts(ixboth&ixleft));
  text(0.01,0,sprintf('    LEFT EYE BINOCULAR VIEWING\n    Best foveation fraction (popt) and NOFF: %1.3f / %1.3f', maxval,log(maxval/(1-maxval))),'FontSize',fontSizeText,'VerticalAlignment','top')

 
 axis off
 export_fig(fullfile(outPATH,patient_id{1}),'-pdf','-transparent','-append')
