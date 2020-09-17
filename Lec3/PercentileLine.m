function line = PercentileLine(xvals, Eta, cdf, percentile)
% PercentileLine Returns percentile line given density diffusion path.
% Inputs:
% xvals: a N1*1 array, the diffusion path steps.
% Eta: a N2*1 array, the Eta grid
% cdf: a N2*N1 matrix, the cumulative density along path
% percentile - percentile to be calculated

% Note we interpolate the closet 5 point to the defined percentile, hence
% require N1>5.

[~,Ix] = sort(abs(cdf-0.01*percentile));
line = zeros(length(xvals),1);
for k = 1:length(xvals)
    line(k) = interp1(cdf(Ix(1:5,k),k),Eta(Ix(1:5,k)+1),0.01*percentile);
end
