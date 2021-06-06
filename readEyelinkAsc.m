function eye = readEyelinkAsc(varargin)
%
% function varargout = readEyelinkAsc(varargin)
% reads Eyelink data file
%
% input: data file name(s)
%
% output:
% d - data: left x y pupil, right x y pupil
% b - blink start and end samples for each eye separately
% t - trial information ntrialsx2: first column - trial onset, second -
%    trial type
% k - key press information of size npressesx2: first column sample,
%   second - keyid
% o - other messages (not key and not trial)
%
% function also works in MONO mode, with corresponding lines being empty
%
% Natalia Zaretskaya 27.04.14
% last modified 14.03.2019
% extended to collect saccade information 04.12.2019

if (nargin < 1)
    
    [fileName,filePath] = uigetfile('*.asc','choose file(s)', 'MultiSelect', 'on');
    
elseif nargin < 2
    
    error('please specify file name and path')
    
else
    
    fileName = varargin{1};
    filePath = varargin{2};
    
end

if ischar(fileName)
    fileName = {fileName};
end


% defaults:
plotFlag = 1;

for i = 1:length(varargin)
    if strcmp(varargin{i}, 'plotFlag')
        plotFlag = varargin{i+1};
    end
end


for i = 1:length(fileName)
    
    fprintf('reading file %s \n', fileName{i});
    
    % read everything line-wise
    fid = fopen(fullfile (filePath, deblank(fileName{i})));
    [sd] = textscan(fid,'%s',...
        'headerlines', 0, 'returnOnError',0,...
        'Delimiter', '\n', 'CollectOutput', false);
    fclose(fid);
    sd = sd{1};
    
    % determine whether a line is data, mesage or other
    isEvent = false(size(sd,1),1);
    isData = false(size(sd,1),1);
    isMessage = false(size(sd,1),1);
    for j = 1:length(sd) % read each line
        if ~isempty(sd{j})
            if ~isempty(strfind('1234567890', sd{j}(1)))
                isData(j) = true;
            elseif ~isempty(strfind('M', sd{j}(1)))
                isMessage(j) = true;
            elseif ~isempty(strfind('E', sd{j}(1)))
                isEvent(j) = true;
            elseif ~isempty(strfind(sd{j}, 'START'))
                startSample = j;
            end
        end
    end
    isData(1:startSample-1) = false;
    isEvent(1:startSample-1) = false;
    isMessage(1:startSample-1) = false;
    dataIdx = find(isData);
    messageIdx = find(isMessage);
    eventIdx = find(isEvent);
    
    % collect data: time leftx lefty leftpupil rightx righty rightpupil ???
    d = zeros(length(dataIdx),8);
    for j = 1:length(dataIdx)
        stringData = sd{dataIdx(j)};
        tmp = textscan(stringData, '%f%f%f%f%f%f%f%f', 'TreatAsEmpty', '.', 'CollectOutput', true);
        d(j,:) = tmp{1};
    end
    
    
    % checke eye events
    % (only blinks implemented)
    bl = [];
    br = [];
    sl = [];
    sr = [];
    for j = 1:length(eventIdx)
        stringData = sd{eventIdx(j)};
%         tmp = textscan(stringData, '%s%s%f%f%*f');
        tmp = textscan(stringData, '%s');
        
        if strcmp(tmp{1}{1}, 'EBLINK')
            if strcmp(tmp{1}{2}, 'L')
                bl = [bl; str2double(tmp{1}{3}) str2double(tmp{1}{4})];
            elseif strcmp(tmp{1}{2}, 'R')
                br = [br; str2double(tmp{1}{3}) str2double(tmp{1}{4})];
            end
        elseif strcmp(tmp{1}{1}, 'ESACC')
            if strcmp(tmp{1}{2}, 'L')
                sl = [sl; str2double(tmp{1}{3}) str2double(tmp{1}{4})...
                    str2double(tmp{1}{5}) str2double(tmp{1}{6})...
                    str2double(tmp{1}{7}) str2double(tmp{1}{8})...
                    str2double(tmp{1}{9})];
            elseif strcmp(tmp{1}{2}, 'R')
                sr = [sr; str2double(tmp{1}{3}) str2double(tmp{1}{4})...
                    str2double(tmp{1}{5}) str2double(tmp{1}{6})...
                    str2double(tmp{1}{7}) str2double(tmp{1}{8})...
                    str2double(tmp{1}{9})];
            end
        end
    end
    b.l = bl;
    b.r = br;
    s.l = sl;
    s.r = sr;
    
    % check custom messages to get trials and key presses
    t = [];
    k = [];
    o = [];
    for j = 1:length(messageIdx)
        stringData = sd{messageIdx(j)};
        tmp = textscan(stringData, '%*s%f%s%f');
        if strcmp(tmp{2}, 'TRIAL')
            t = [t; tmp{1} tmp{3}];
        elseif strcmp(tmp{2}, 'KEY')
            k = [k; tmp{1} tmp{3}];
        elseif strcmp(tmp{2}, 'OFFS')
            o = [o; tmp{1} tmp{3}];
        end
    end
    
    
end

eye.d = d;
eye.b = b;
eye.s = s;
eye.k = k;
eye.t = t;
eye.o = o;

if plotFlag
    figure;
    plot(d(:,1), d(:,2:end-1));
    hold on
    if ~isempty(k)
        plot(k(:,1), k(:,2), 'k.')
    end
    if ~isempty(t)
        plot(t(:,1), t(:,2), 'b.')
    end
    if ~isempty(s)
        plot(s.l(:,1), ones(size(s.l(:,1))), 'r.')
    end
    if ~isempty(b)
        plot(b.l(:,1), ones(size(b.l(:,1))), 'g.')
    end
    xlim([d(1,1) d(end,1)])
    title(fileName)
    xlabel('time (ms)')
    legend({'xl','yl','pl','xr','yr','pr', 'key', 'trial', 'lsacc', 'lblink'})
    grid on
end

end


