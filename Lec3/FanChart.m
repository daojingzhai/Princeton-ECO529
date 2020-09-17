function [lh, ph] = FanChart(T, Eta, pdfORcdf, type, prctiles)
% FANCHART  Create a fan chart visualization
% A fan chart is a plot of time-varying distribution percentiles shown as
% shaded bands around a central (median/mean) line. 
% REQUIRED INPUTS:
% T: a N1*1 array, the diffusion path steps.
% Eta: a N2*1 array, the Eta grid
% pdfORcdf: a N2*N1 matrix, the density or cumulative density along path.
% type: 'cdf' or 'pdf'

% OPTIONAL INPUTS:
% prctiles - an array of percentiles to calculate (default: 5:5:95)

% Parse inputs
if nargin < 5 || isempty(prctiles)
    prctiles = 5:5:95;
end
ip = inputParser;
addParamValue(ip, 'parent', gca, @(x)ishandle(x)&&strcmp(get(x,'type'),'axes')); %#ok<*NVREPL>
addParamValue(ip, 'alpha', 1, @(x)isscalar(x)&&isnumeric(x));
addParamValue(ip, 'colormap', @boeRedMap, @(x)isa(x,'function_handle')||ischar(x)||iscell(x)&&ischar(x{1}));
parse(ip);
results = ip.Results;
parent = results.parent;
alpha = results.alpha;
cmapFun = results.colormap;

% Calculate fan chart bands
if strcmp(type, 'pdf')
    % cdf for all times after shock
    cdf = cumsum((pdfORcdf(:,2:end)').*((Eta(2:end)-Eta(1:end-1))*ones(1,size(pdfORcdf,1))));
else
    cdf = pdfORcdf;
end

bands = zeros(length(T),length(prctiles));
for i = 1:length(prctiles)
    band_temp = PercentileLine(T, Eta, cdf, prctiles(i));
    bands(:,i) = band_temp;
end

% Calculate centerline
% centerline, choose median
centerline = PercentileLine(T, Eta, cdf, 50);
% centerline, choose 25% percentile
line25 = PercentileLine(T, Eta, cdf, 25);
% centerline, choose 75% percentile
line75 = PercentileLine(T, Eta, cdf, 75);

% Xvalues for bands
xplot = T([1:end end:-1:1]);

% Create colormap
ncolors = floor(length(prctiles)/2);
if iscell(cmapFun)
    col = feval(cmapFun{1}, cmapFun{2:end}, ncolors);
else
    col = feval(cmapFun, ncolors);
end

% Create plot
if verLessThan('matlab', '8.4')
    ph = zeros(ncolors,1);
else
    ph = gobjects(ncolors,1);
end
for i = 1:ncolors
    ph(i) = patch(xplot, [bands(:,0+i);  flipud(bands(:,end-i+1))]', col(i,:) ,...
        'EdgeColor', col(i,:), 'Parent', parent, 'FaceAlpha', alpha);
end
lh = line(T, centerline, 'LineWidth', 2, 'Color', 'k', 'Parent', parent);
line(T, line25, 'LineStyle', '--', 'LineWidth', .5, 'Color', 'k', 'Parent', parent);
line(T, line75, 'LineStyle', '--', 'LineWidth', .5, 'Color', 'k', 'Parent', parent);
ph = flipud(ph);


function map = boeRedMap(ncolors)

colors = [254 230 222
    252 211 196
    250 190 171
    248 170 148
    246 151 127
    244 132 108
    243 113 92
    241 91 75
    237 27 46];
map = colors/256;
if nargin > 0
    map = interp1(1:size(map,1), map, 1:ncolors);
end
