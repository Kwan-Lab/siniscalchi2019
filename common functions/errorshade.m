function errorshade(x,h,l,c,t)
% % errorshade(x,h,l,c) %
%PURPOSE:   Plot shaded area, e.g. for confidence intervals
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   x:  the independent variable
%   h:  the upper bound for the dependent variable
%   l:  the lower bound for the dependent variable
%   c:  the color of the shading, e.g., [0.7 0.7 0.7] for gray
%   t:  fraction solid/translucent (0 to 1)

%%
if nargin < 5
    facealpha = 1; %not translucent
else
    facealpha = t; %make the shading translucent
end

x=x(:)';
h=h(:)';
l=l(:)';

%remove the NaN values
badIdx = isnan(h) | isnan(l);
x=x(~badIdx);
h=h(~badIdx);
l=l(~badIdx);

%plot the shading
h=fill([x fliplr(x)],[h fliplr(l)],c,'LineStyle','none');
set(h,'facealpha',facealpha);