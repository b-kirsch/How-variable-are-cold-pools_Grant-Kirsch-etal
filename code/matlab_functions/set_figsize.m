function set_figsize( ifig, figsize )
%
% author: Leah Grant (leah.grant@colostate.edu), Aug. 2013
%
% function set_figsize( ifig, figsize )
%   sets figure print size/properties, and moves/resizes
%   Matlab figure display
% 
% ifig is the figure handle
%
% figsize is a 2-element array with the paper size 
%   [width height] given in inches
% if figsize is not passed in, the default paper size 
%   will be used (8.5x11) and the function will simply
%   move the figure to the lower-left corner of the screen
%
% also sets paper position s.t. a half-inch border exists around
%   the paper size (for printing)
%
% example: set_figsize( gcf, [15 10] )
%          set_figsize( gcf, [] )   % to use default 8.5x11 size



% set paper units to inches (this should be the default)
set(ifig,'PaperUnits','Inches')

% if the figsize wasn't passed in, get the current fig size
if isempty(figsize)
    figsize=get(ifig,'PaperSize');
end

% move figure to lower left corner and scale to figsize
% set units to inches
set(ifig,'Units','Inches')
set(ifig,'OuterPosition', [0 0 figsize] )

% resize figure so it prints correctly
set(ifig,'PaperSize',figsize)
% set paper position with half inch margins
% [up, left from corner, length, width] units in inches
set(ifig,'PaperPosition',[0.5 0.5 figsize-1])
