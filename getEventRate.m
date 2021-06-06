function [ rate ] = getEventRate( events, epochs, samplingRate )
% computes event rate withing given epochs
%
% INPUTS
% events - matrix of n_events x 2 with indices of event onsets and offsets (e.g. eye blinks, saccades)
% epochs - matrix of n_epochs x 3 with indices of onsets abd offsets of epochs
%  plus condition ID in the 3d column
% samplingRate
%
% OUTPUTS
% rate n_conditions x 1 - rate of events in each condition type
%
% Natalia Zateskaya 07.2014
% update 26.09.2014

eventCount = [];
timeCount = [];


epochTypes=unique(epochs(:,3));

for i = 1:length(epochTypes)
    
    thisEpochType = epochs(:,3) == epochTypes(i);
    
    [aidx, ~] = epochOverlap(events(:,[1 2]), epochs(thisEpochType, [1 2]) );
%     eventCount(trl,i) = numel(aidx);
%     timeCount(trl,i) = sum(epochs(thisEpochType, 2)-epochs(thisEpochType, 1));
    eventCount(i) = numel(aidx);
    timeCount(i) = sum(epochs(thisEpochType, 2)-epochs(thisEpochType, 1));
    
end




% sum across tirals
eventCount = sum(eventCount,1);
timeCount = sum(timeCount,1);

% transform to Hz
rate = eventCount./(timeCount/samplingRate);

end

