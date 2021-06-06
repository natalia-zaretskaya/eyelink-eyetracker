function idata = interpolateEvents(x, data, events, addsamples)
% function interploateEvents interpolates samples that are within a certain
% "event". An event is defined by smaple at the start and sample at the
% end.
% Usage:  idata = interpolateEvents(x, data, events, addsamples)
% INPUTS: 
% x - vector of the x-axis (e.g. time, sample nubmer)
% data - vector time course, e.g. left-eye x coordinate
% events - n_events x 2 matrix with onset and offset time of each event
% addsamples - integer number of samples to be added before and after the
% identified event.
% OUTPUTS:
% idata - vector with interplated samples
%
% Natalia Zaretskaya 04.14

valid = true(size(x));
events(:,1) = events(:,1)-addsamples;
events(:,2) = events(:,2)+addsamples;

valid(isnan(data)) = false; % also set Nan to invalid

for i = 1:size(events,1)
    valid(find(x==events(i,1)):find(x==events(i,2))) = false; % invalid
end

ipoints = interp1(find(valid), data(valid,:), find(~valid),  'linear', 'extrap');

idata = data;

idata(~valid,:) = ipoints;

end