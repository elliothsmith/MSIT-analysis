function axisHandle = plotmultipleaxes(k,plotDIMx,plotDIMy,plotBorder,figHandle)
%PLOTMULTIPLEAXES Create axes in tiled positions
%   H = plotMultipleAxes(k,M,N,Border,figHandle), breaks the Figure window
%   specified by figHandle into an N-by-M matrix of axes, selects 
%   the k-th axes for the current plot, and returns the axes handle.
%   
%   The Border variable allows the user to customize the 
%   size of subplots and the gaps between subplots. Try 0.002.


% Version Date: 20110503
% Author: Tyler Davis

% updated: 20130828 (Elliot Smith)
% - added help text, added plot border to input arguments.

% updated: 20130906 (Elliot Smith)
% - fixed dimensions in help text, changed capitalization

plotIndx = plotBorder:plotBorder+(1-(plotBorder*(plotDIMx+1)))/plotDIMx:1;
plotIndy = plotBorder:plotBorder+(1-(plotBorder*(plotDIMy+1)))/plotDIMy:1;
plotIndx = plotIndx(plotIndx<1);
plotIndy = plotIndy(plotIndy<1);
[xIdx,yIdx] = meshgrid(plotIndx,plotIndy);
yIdx = flipud(yIdx);

axisHandle = subplot('Position',[xIdx(k) yIdx(k) (1-(plotBorder*(plotDIMx+1)))/plotDIMx (1-(plotBorder*(plotDIMy+1)))/plotDIMy],'Parent',figHandle);