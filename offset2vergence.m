function [ vergenceAngle, y_intersect ] = offset2vergence( offset, distance_eyeToEye, distance_eyeToScreen )
% function [ vergenceAngle ] = offset2vergence( offset, distance_eyeToEye, distance_eyeToScreen )
% Converts offset between the two images on the screen to a perceived
% distance and ideal vergence angle for that distance. Based on some basic
% school geometry.
% INPUTS:
% offset - a vector with offset values. Positive offset - displacement of
% the left-eye image to the right, objects become closer. Negative offset -
% displacement of left-eye image to the LEFT, objects become farther away.
% distance_eyeToEye - distance between the two pupils
% distance_eyeToScreen
% all distances and the offset have to have the same units, e.g. all in cm
% OUTPUTS: 
% vergenceAngle - vector with ideal vergence angles, one for every offset
% value
% y_intersect - theoretical distance of focal point from observer

vergenceAngle = zeros(size(offset,1), size(offset,2));

for i = 1:length(offset)
    % define two lines of sight
    % left eye
    x1 = [-distance_eyeToEye/2 offset(i)];
    y1 = [0 distance_eyeToScreen];
    % right eye
    x2 = [distance_eyeToEye/2 -offset(i)];
    y2 = [0 distance_eyeToScreen];
    
    % find intersection point of two gaze lines
    p1 = polyfit(x1,y1,1);
    p2 = polyfit(x2,y2,1);
    %calculate intersection
    x_intersect = fzero(@(x) polyval(p1-p2,x),3);
    y_intersect = polyval(p1,x_intersect);
    
    % determine vergence angle using known y-distance
    vergenceAngle(i) = atand((distance_eyeToEye/2)/y_intersect)*2;
    
    % plot this for debugging
    %     figure
    %     line(x1,y1);
    %     hold on;
    %     line(x2,y2);
    %     plot(x_intersect,y_intersect,'r*')
    %     refline(0, distance_eyeToScreen)
    
end % for i
end

