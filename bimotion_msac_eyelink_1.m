
clear all
close all force
set(0,'DefaultTextInterpreter','None')

% to do:
% - conduct the same experiment with the texturized version (longer local durations)
% - try to detect msaccs without interpolating blinks
% - stats of perceptual durations
% compute blink rate

plotFlag = 0;

root_dir = '/Data/studies/bimotion_msac/data';
path(path, '/Users/Natalia/Dropbox/stereoCode/');

f = filesep;

thex = [-8*1000 : 0.5*1000 : 8*1000]; % in samples
conditionColors = 'rgbmck';

% ------- subject info ------- %
s = 1;
% subj(s).init = 'nz_kza';
% subj(s).eyefiles = [1:4];
% subj(s).perceptKeys = [37; 39];
% s = s+1;

% subj(s).init = 'nz';
% subj(s).eyefiles = [1:6]; % file 3 and 5 cant be read (NaNs)
% subj(s).perceptKeys = [37; 39];
% s = s+1;
% 
% subj(s).init = 'nz_squares';
% subj(s).eyefiles = [1:8];
% subj(s).perceptKeys = [37; 39];
% s = s+1;
% 
% subj(s).init = 'nz2';
% subj(s).eyefiles = [1:3];
% subj(s).perceptKeys = [39; 37];
% s = s+1;
% 
% subj(s).init = 'pg';
% subj(s).eyefiles = [1:4];
% subj(s).perceptKeys = [39; 37];
% s = s+1;
% 
% subj(s).init = 'pg230614';
% subj(s).eyefiles = [1:3];
% subj(s).perceptKeys = [39; 37];
% s = s+1;

% subj(s).init = 'mn230714';
% subj(s).eyefiles = [1:3];
% subj(s).perceptKeys = [39; 37];
% s = s+1;

% subj(s).init = 'as240614';
% subj(s).eyefiles = [1:3];
% subj(s).perceptKeys = [39; 37];
% s = s+1;

% subj(s).init = 'nz_diamond';
% subj(s).eyefiles = [1:6];
% subj(s).perceptKeys = [37; 39];
% s = s+1;

% subj(s).init = 'pg220714diamond';
% subj(s).eyefiles = [1];
% subj(s).perceptKeys = [39; 37];
% s = s+1;
% 
subj(s).init = 'hg230714diamond';
subj(s).eyefiles = [1:4];
subj(s).perceptKeys = [39; 37];
s = s+1;

for s = 1:length(subj) % loop over subjects
    
    subj(s).sacPSTH = cell(1,3);
    subj(s).blkPSTH = cell(1,3);
    subj(s).pupPSTH = cell(1,3);
    subj(s).pupPBTH = []; % periblink pupil width
    subj(s).sacData = [];
    subj(s).sacRate = [];
    subj(s).blkRate = [];
    subj(s).percepts = [];
    
    eyefnames = dir(fullfile(root_dir, subj(s).init, ['*.asc']));
    
    eyefnames = struct2cell(eyefnames); eyefnames = eyefnames(1,:); eyefnames = str2mat(eyefnames);
    
    fprintf('found files: \n');
    disp(eyefnames);
    
    for f = subj(s).eyefiles
        
        data = readEyelinkAsc(eyefnames(f,:), fullfile(root_dir, subj(s).init), 'plotFlag', plotFlag);        
        
        
        % saccade data
        sac = myEngbertMsac(data, plotFlag);        
        
        % --- trigger on key data --- %
        if ~isempty(data.t)
            exp.trialStartTime = data.t(:,1);
            exp.trialEndTime = [data.t(2:end,1); data.d(end,1)];
        else % temporary solution
            exp.trialStartTime = [0 60 75 135 150].*1000+data.d(1,1);
            exp.trialEndTime = [60 75 135 150 165].*1000+data.d(1,1);
        end
            
        key.idDown = data.k(:,2);
        key.idUp = data.k(:,2);
        key.timeDown = data.k(:,1);
        key.timeUp = data.k(:,1);
        
        
        % defaults:
        settings.responseType      = 'press-release';
        settings.rejectLastPercept = 1;
        settings.ignoreSameKeyPressedAgain = 0;
        settings.perceptKeys = subj(s).perceptKeys; % beware of differences on windows and mac!
        
        [~, allPercepts]= analyzeBistableKeys(exp, key, 'plotFlag', plotFlag, 'settings', settings);
        
        % make one tiral out of several
        allPercepts = {vertcat(allPercepts{:})};
        
        c = [];
        c.minimalDuration         = 0;
        c.maximalDuration         = Inf;
        c.previousDifferent       = 1;
        c.previousMinimalDuration = 0;
        c.previousMaximalDuration = Inf;
        c.nextDifferent           = 0;
        c.nextMinimalDuration     = 0;
        c.nextMaximalDuration     = Inf;
        
        cleanedPercepts = cleanBistableKeys(allPercepts, c, 'plotFlag', 0);
        subj(s).percepts = [subj(s).percepts; cleanedPercepts];
        
        c = [];
        c.alignTo = 'midPoint';
        c.limitType = 'scaled';
        c.limits = [-0.5 0.5];
        epochs = defineEpochsFromPercepts(cleanedPercepts, c, 'plotFlag', 0);
        
        
        for k = 1:length(subj(s).perceptKeys)
            keyIndices{k} = cleanedPercepts{1}(cleanedPercepts{1}(:,3)==subj(s).perceptKeys(k),1);
        end
        
        uniqueTriggers = subj(s).perceptKeys;
        triggerIndices = keyIndices;
        
        % saccade rate for each percept type
        sacRate = getEventRate( sac, epochs, 1000);
        subj(s).sacRate = [subj(s).sacRate; sacRate];
        
        % blink rate for each percept type
        blkRate = getEventRate( data.b.l-data.d(1,1)+1, epochs, 1000);
        subj(s).blkRate = [subj(s).blkRate; blkRate];
        

        
%         % --- trigger on trial onsets --- %
%         presentTrials = unique(data.t(:,2));
%         for k = 1:length(presentTrials)
%             trialIndices{k} = data.t(data.t(:,2)==presentTrials(k),1)-data.d(1,1)+1;
%         end
%         uniqueTriggers = presentTrials;
%         triggerIndices = trialIndices;
%         
%         % saccade rate for each percept type
%         sacRate = getEventRate( sac, {[trialIndices{1} trialIndices{1}+12*1000 ones(size(trialIndices{1}))]}, 1000);
         
        
        if ~isempty(sac)        
            
            % create artificial saccade timecourse with ones at sac onsets
            sacTc = zeros(size(data.d,1),1);
            sacTc(sac(:,1)) = 1; % onsets
            
            % create artificial blink timecourse with ones at blink onsets
            blkTc = zeros(size(data.d,1),1);
            blkTc(data.b.l(:,1)-data.d(1,1)+1) = 1; % onsets
            
            % pupil time course
            pl = interpolateEvents(data.d(:,1), data.d(:,4), data.b.l, 100);
            pr = interpolateEvents(data.d(:,1), data.d(:,7), data.b.r, 100);
            pl = zscore(detrend(pl));
            pr = zscore(detrend(pr));
            pupTc = mean([pl pr],2);
            
            % periblink pupil width        
            [~, ~, pupPBTH] = mypsth(pupTc, data.b.l(:,1)-data.d(1,1)+1, thex, 'doWhat', 'average');
            subj(s).pupPBTH = [subj(s).pupPBTH; pupPBTH]; 
            
            for k = 1:length(uniqueTriggers)
                [~, ~, sacPSTH] = mypsth(sacTc, triggerIndices{k}, thex, 'doWhat', 'count');
                [~, ~, blkPSTH] = mypsth(blkTc, triggerIndices{k}, thex, 'doWhat', 'count');
                [~, ~, pupPSTH] = mypsth(pupTc, triggerIndices{k},thex, 'doWhat', 'average');
                subj(s).sacPSTH{k} = [subj(s).sacPSTH{k}; sacPSTH];
                subj(s).blkPSTH{k} = [subj(s).blkPSTH{k}; blkPSTH];
                subj(s).pupPSTH{k} = [subj(s).pupPSTH{k}; pupPSTH];               
            end
            
            subj(s).sacData = [subj(s).sacData; sac];
            
        else
            warning('No msacs found in this file - check your data. Skipping...')
            
        end
        
    end % for f
    
    
    figure;
    for k = 1:length(uniqueTriggers)
        plot(thex, sum(subj(s).sacPSTH{k}), conditionColors(k)); hold on;
    end
    title('saccade count')
    
    figure;
    for k = 1:length(uniqueTriggers)
        plot(thex, sum(subj(s).blkPSTH{k}), conditionColors(k)); hold on;
    end
    title('blink count')
    
    figure;
    for k = 1:length(uniqueTriggers)
        plot(thex, nanmean(subj(s).pupPSTH{k}), conditionColors(k)); hold on;
        shadedErrorBar(thex, nanmean(subj(s).pupPSTH{k}), sem(subj(s).pupPSTH{k}), conditionColors(k));
    end
    title('pupil area')
    
    
    figure;
    plot(thex, nanmean(subj(s).pupPBTH)); hold on;
    shadedErrorBar(thex, nanmean(subj(s).pupPBTH), sem(subj(s).pupPBTH));    
    title('pupil area at blink onset')
    
    
    figure;
    subplot(1,3,1)
    scatter(subj(s).sacData(:,3), subj(s).sacData(:,4))
    xlabel('velocity')
    ylabel('amplitude')
    
    subplot(1,3,2)
    scatter(subj(s).sacData(:,4), subj(s).sacData(:,2)-subj(s).sacData(:,1));
    xlabel('amplitude')
    ylabel('duration')
    
    subplot(1,3,3)
    polar(subj(s).sacData(:,5), subj(s).sacData(:,3), 'o');
    title('anglular orientation')
    
end % for s

save tmpdata.mat subj

%% group statistics
y = [];
y.meanSacRate = [];
y.meanBlkRate = [];
y.meanPupPBTH = [];
y.meanPupPSTH{1} = [];
y.meanPupPSTH{2} = [];
y.meanSacPSTH{1} = [];
y.meanSacPSTH{2} = [];
y.meanBlkPSTH{1} = [];
y.meanBlkPSTH{2} = [];
y.medianDuration = [];

for s = 1:length(subj)
    
    y.meanSacRate = [y.meanSacRate; nanmean(subj(s).sacRate,1)];
    y.meanBlkRate = [y.meanBlkRate; nanmean(subj(s).blkRate,1)];
    y.meanPupPBTH = [y.meanPupPBTH; nanmean(subj(s).pupPBTH,1)];
    
    if iscell(subj(s).percepts)
        perceptData = vertcat(subj(s).percepts{:});
    else
        perceptData = subj(s).percepts;
    end
    
    perceptDurations = perceptData(:,2)-perceptData(:,1);
    perceptIds = perceptData(:,3);
    
    y.medianDuration = [y.medianDuration;...
        [median(perceptDurations(perceptIds==subj(s).perceptKeys(1)))...
        median(perceptDurations(perceptIds==subj(s).perceptKeys(2)))]];
    
    for c = 1:length(y.meanPupPSTH)
        y.meanPupPSTH{c} = [y.meanPupPSTH{c}; nanmean(subj(s).pupPSTH{c},1)];
        y.meanSacPSTH{c} = [y.meanSacPSTH{c}; nanmean(subj(s).sacPSTH{c},1)];
        y.meanBlkPSTH{c} = [y.meanBlkPSTH{c}; nanmean(subj(s).blkPSTH{c},1)];
    end
    
end
%%
[~,~,~,stats] = ttest(y.meanSacRate(:,1), y.meanSacRate(:,2));
[~,~,~,stats] = ttest(y.meanBlkRate(:,1), y.meanBlkRate(:,2));

%%
% mean psth
measuresOfInterest = fieldnames(y);
for i = 1:length(measuresOfInterest)
    Y = y.(measuresOfInterest{i});
    figure;
        
    if ~isempty(strfind(measuresOfInterest{i}, 'Rate'))|| ~isempty(strfind(measuresOfInterest{i}, 'Duration'))
        [~,p,ci,stats] = ttest(Y(:,1), Y(:,2));
        bar(1:size(Y,2), mean(Y,1)); hold on; 
        my_errorbar(1:size(Y,2), mean(Y,1), sem(Y));
        text(1,0, sprintf('p = %f, t = %f', p, stats.tstat))
    end
    
    
    if ~isempty(strfind(measuresOfInterest{i}, 'PSTH'))
        for c = 1:length(Y)
            plot(thex, mean(Y{c}), conditionColors(c)); hold on;
            shadedErrorBar(thex, mean(Y{c}), sem(Y{c}), conditionColors(c));
            grid on
        end
    end
    
    if ~isempty(strfind(measuresOfInterest{i}, 'PBTH'))
        plot(thex, mean(Y), conditionColors); hold on;
        shadedErrorBar(thex, mean(Y), sem(Y));
        grid on
    end
    title(measuresOfInterest{i});
end

