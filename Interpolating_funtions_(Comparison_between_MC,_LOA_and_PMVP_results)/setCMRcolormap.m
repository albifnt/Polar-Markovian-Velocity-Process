% Daniel W. Meyer
% Institute of Fluid Dynamics, ETH Zurich
% March 2016
%
% use colormap that shows nicely in color and b/w
% cflag: set true for continuous colormap
function setCMRcolormap(cflag)
% discrete CMR colormap with 9 values
% cmapCMR =[.95 .95 0.95; .9 .9 .5; .9 .75 .1; .9 .5 0; 1 .25 .15; .6 .2 .50; .3 .15 .75; .15 .15 .5; 0 0 0];
cmapCMR =[0 0 0; .15 .15 .5; .3 .15 .75; .6 .2 .50; 1 .25 .15; .9 .5 0; .9 .75 .1; .9 .9 .5; .95 .95 0.95];
if ((nargin == 1) && (~cflag)) % discrete CMR colormap with 9 values
    colormap(cmapCMR)
else % continous CMR colormap with 64 base values
    % interpolate to 64 'smooth' values
    xi = 1:8/63:9; % 64 color levels instead of 9
    x = 1:9;
    for i = 1:3
      cmapCMRs(:,i) = spline(x,cmapCMR(:,i),xi)'; % spline fit intermediate values
    end
    % eliminate spurious values outside of rangecolormap(cmapCMR)
    cmapCMRs = cmapCMRs - min(min(cmapCMRs)); % lower bound is zero
    cmapCMRs = cmapCMRs / max(max(cmapCMRs)); % upper bound is one
    colormap(cmapCMRs)
end
