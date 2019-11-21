function hout=suptitle(varargin)
%SUPTITLE Puts a title above all subplots.
%	SUPTITLE('text') adds text to the top of the figure
%	above all subplots (a "super title"). Use this function
%	after all subplot commands.

% Drea Thomas 6/15/95 drea@mathworks.com

% Warning: If the figure or axis units are non-default, this
% will break.

str = varargin{1};
if length(varargin)>1
	params = varargin(2:2:end);
	vals = varargin(3:2:end);
	if length(params)~=length(vals)
		error('Parameter and value pairs are not balanced');
	end
else
	params = {};
	vals = {};
end

% Parameters used to position the supertitle.

% Amount of the figure window devoted to subplots
plotregion = .92;

% Y position of title in normalized coordinates
titleypos  = .95;

% Fontsize for supertitle
fs = get(gcf,'defaultaxesfontsize')+4;

% Fudge factor to adjust y spacing between subplots
fudge=1;

haold = get(gcf,'currentAxes');
figunits = get(gcf,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.
