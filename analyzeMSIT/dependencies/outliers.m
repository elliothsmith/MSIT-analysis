function [idx,iprange,fence]=outliers(data,varargin)
% [idx,iprange,inner_fence]=outliers(data,tile,thresh)
% 
% Find outliers in a vector.
% skellis 12/6/2010
%
% method: http://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm
%
% Inputs:
%   data   - vector of observations from a single variable
%   tile   - percentiles used (optional; default [25 75])
%   thresh - multiplied by the inter-percentile range to get the inner 
%            fence (optional; default 1.5)
%
% Outputs:
%   idx         - index of outliers in the input data vector
%   iprange     - the inner-percentile range (see "help iqr")
%   inner_fence - the thresholds calculated for outlying values

% defaults
tile=[25 75];
thresh=1.5;
sides='both';

fprintf('\nIndexing %s sided outliers...',sides)

% user-specified values
if(nargin>3)
    sides=varargin{3};
end
if(nargin>2)
    thresh=varargin{2};
end
if(nargin>1)
    tile=varargin{1};
end

% find the the percentiles specified by "tile" and inter-percentile range
percentiles=prctile(data,tile); 
iprange=diff(percentiles,1);

% calculate fences (quartile +/- thresh*iqr)
fence=[percentiles(1)-thresh*iprange; percentiles(2)+thresh*iprange];

% find points outside the inner fences
if(strcmpi(sides,'both'))
    idx=(data<fence(1)|data>fence(2));
elseif(strcmpi(sides,'left'))
    idx=(data<fence(1));
elseif(strcmpi(sides,'right'))
    idx=(data>fence(2));
end