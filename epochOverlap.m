function [aidx, bidx] = epochOverlap(A, B)
%
% detect overlap between two set of epoch, e.g. artifact
% and trial which contains artifactual samples
% solution from
% http://stackoverflow.com/questions/10301718/calculating-overlap-in-mx2-and-nx2-matrices
%
% inputs:
% A - fieldtrip format structure with cfg.trl matrix with trial information
% B - matrix with start and end indices of artifacts
%
% outputs:
% aidx - indices to A that overlap with B
% bdix - indices to B that overlap with A
% Natalia 02.2012
% update 06.2014 - more general version


overlap = @(x, y)y(:, 1) < x(:, 2) & y(:, 2) > x(:, 1);

[tmp1, tmp2] = meshgrid(1:size(A, 1), 1:size(B, 1));

M = reshape(overlap(A(tmp1, :), B(tmp2, :)), size(B, 1), [])';

[aidx, bidx] = find(M);


end