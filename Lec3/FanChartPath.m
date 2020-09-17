function [lh, ph] = FanChartPath(xvals, paths, centerType, prctiles, varargin)
% FANCHART  Create a fan chart visualization
%
% A fan chart is a plot of time-varying distribution percentiles shown as
% shaded bands around a central (median/mean) line. It is useful for plotting paths of a
% Monte-Carlo simulation of time series
% 
% Syntax:
%
% [centObj, bandsObj] = fanChart(xvals, paths, centerType, prctiles, 'param', value, ...)
%
% Inputs:
% xvals - the X-axis values for the data to be visualized. If empty, xvals = 1:size(paths,1)
% paths - a nSteps-by-nTrials matrix of simulated paths. 
% centerType - string 'mean' or 'median' (default: 'median')
% prctiles - an array of percentiles to calculate (default: 5:5:95)
%
% Optional parameter-value inputs:
% parent - the axes object to draw the chart in (default: gca)
% alpha -  the transparency coefficient for the bands (default: 1 (no
%          transparency))
% colormap - a string, cell array or function handle to a colormap function (default: boeRedMap).
%            the cell array syntax should contain a string name of the colormap function as the first
%            element and optional arguments to that function as additional elements
%            Available colormap functions in fanChart are boeRedMap,
%            shadesOfRed and shadesOfColor
%
% See also fanChart_examples
% Copyright 2014 The MathWorks, Inc.
% Parse inputs
if nargin < 3 || isempty(centerType)
    centerType = 'median';
end
if nargin < 4 || isempty(prctiles)
    prctiles = 5:5:95;
end
ip = inputParser;
addParamValue(ip, 'parent', gca, @(x)ishandle(x)&&strcmp(get(x,'type'),'axes')); %#ok<*NVREPL>
addParamValue(ip, 'alpha', 1, @(x)isscalar(x)&&isnumeric(x));
addParamValue(ip, 'colormap', @boeRedMap, @(x)isa(x,'function_handle')||ischar(x)||iscell(x)&&ischar(x{1}));
parse(ip, varargin{:});
results = ip.Results;
parent = results.parent;
alpha = results.alpha;
cmapFun = results.colormap;
% Calculate fan chart bands
% Initially we will do this with percentiles
bands = prctile(paths, prctiles, 2);
% Calculate centerline
centerline = feval(lower(centerType), paths, 2);
% Xvalues for bands
xplot = xvals([1:end end:-1:1]);
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
    ph(i) = patch(xplot, [bands(:,0+i);  flipud(bands(:,end-i+1))]', col(i,:) ,'EdgeColor', col(i,:), 'Parent', parent, 'FaceAlpha', alpha);
end
lh = line(xvals, centerline, 'LineWidth', 2, 'Color', 'k', 'Parent', parent);
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
function map = shadesOfRed(ncolors)
%shades = linspace(.8,0,ncolors);
%shades = log10(linspace(10^(.8),10^(0),ncolors));
shades = log(linspace(exp(.88),exp(.11),ncolors));
%col = zeros(ncolors, 3);
%col(:,1) = shades;
map = repmat(shades(:),1,3);
map(:,1) = 1;
function map = shadesOfColor(color, ncolors)
HSV = rgb2hsv(color);
shades = log(linspace(exp(.11),exp(.88),ncolors));
%shades = linspace(.11,.88,ncolors);
colors = zeros(ncolors, 3);
colors(:,1) = HSV(1);
colors(:,2) = shades;
colors(:,3) = HSV(3);
map = hsv2rgb(colors);