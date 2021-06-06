function [ epochs ] = defineEpochsFromPercepts( percepts, c, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% c.alignTo = 'onset';
% c.limitType = 'constant';
% c.limits = [-1000 1000];

% c.alignTo = 'midPoint';
% c.limitType = 'constant';
% c.limits = [-1000 1000];
% 
% c.alignTo = 'onset';
% c.limitType = 'scaled';
% c.limits = [-0.5 0.5]; % half of the previous percept, half of the next
% 
% c.alignTo = 'midPoint';
% c.limitType = 'scaled';
% c.limits = [-0.5 0.5];
% 
% c.alignTo = 'onset';
% c.limitType = 'constantOffset';
% c.limits = [-1000 1000]; % half of the previous percept, half of the next
% 
% c.alignTo = 'midPoint';
% c.limitType = 'constantOffset';
% c.limits = [-1000 1000];

plotFlag = 1;

% check optional arguments
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'plotFlag')
        plotFlag = varargin{i+1};
    end
    
end

fprintf('=> Defining epochs from percepts:\n')

for trl = 1:length(percepts)
    
    if strcmp(c.alignTo, 'onset')
        triggers = percepts{trl}(:,1);
    elseif strcmp(c.alignTo, 'midPoint')
        triggers = percepts{trl}(:,1)+(percepts{trl}(:,2)-percepts{trl}(:,1))/2;
    else
        error('Undefined alignment method')
    end
    
    
    durations = percepts{trl}(:,2)-percepts{trl}(:,1);
    
    if strcmp(c.limitType, 'constant');
        
        epochs{trl} = [triggers+c.limits(1) triggers+c.limits(2) percepts{trl}(:,3)];
        
    elseif strcmp(c.limitType, 'scaled')
        
        if strcmp(c.alignTo, 'midPoint')
            epochs{trl} = [triggers+durations*c.limits(1)...
                triggers+durations*c.limits(2) ...
                percepts{trl}(:,3)];
            
        elseif strcmp(c.alignTo, 'onset')
            epochs{trl} = [triggers+[NaN; durations(1:end-1)]*c.limits(1)...
                triggers+durations*c.limits(2) ...
                percepts{trl}(:,3)];
        end
        
        
    elseif strcmp(c.limitType, 'constantOffset')
        
        if strcmp(c.alignTo, 'midPoint')
            epochs{trl} = [triggers-durations/2-c.limits(1)...
                triggers+durations/2-c.limits(2) ...
                percepts{trl}(:,3)];
            
        elseif strcmp(c.alignTo, 'onset')
            epochs{trl} = [triggers-[NaN; durations(1:end-1)]-c.limits(1)...
                           triggers+durations-c.limits(2) ...
                           percepts{trl}(:,3)];
        end
        
    end
    
    
    % clean the resulting samples
    % delete negative durations, if any;
    isNegative = (epochs{trl}(:,2) - epochs{trl}(:,1))<0;
    epochs{trl}(isNegative,:) = [];
    fprintf('%u negative durations deleted\n', numel(find(isNegative)));
    
    [nanRows, nanColumns] = find(isnan(epochs{trl}));
    epochs{trl}(nanRows,:) = [];
    fprintf('%u NaN-containing rows deleted\n', numel(nanRows));
    
end % trl




if plotFlag
    
    figure;
    
    for trl = 1:length(epochs)
        
        subplot(length(epochs), 1, trl);
        
        perceptTypes = unique(percepts{trl}(:,3));
        
        for p = 1:length(perceptTypes)
            
            % plot original data
            idx = find(percepts{trl}(:,3) == perceptTypes(p));
            for i = 1:length(idx)
                plot([percepts{trl}(idx(i),1) percepts{trl}(idx(i),2)], [p p], 'LineWidth', 3, 'Color', 'b')
                hold on
            end
            ylim([0 length(perceptTypes)+1]);
            
            
            % plot clean data
            if ~isempty(epochs{trl})
                
                idx = find(epochs{trl}(:,3) == perceptTypes(p));
                
                for i = 1:length(idx)
                    plot([epochs{trl}(idx(i),1) epochs{trl}(idx(i),2)], [p p]+0.25, 'LineWidth', 3, 'Color', 'r')
                end
                
                hold on
                ylim([0 length(perceptTypes)+1]);
                
            end
            
        end
        
    end
    
end





end

