function ax = multipan( ifig, nr, nc, panel, varargin )
% function multipan( ifig, nr, nc, panel, varargin )
%
% author: Leah Grant (leah.grant@colostate.edu), Aug 2017
%
% ifig = figure handle
% nr: number of rows
% nc: number of columns
% panel: panel number on which to plot
% 
%   example:
% multipan( gcf, 2, 3, 5 )
%   will set up axes for the fifth panel in a 2x3 plot like so:
%  -----------
% | 1 | 2 | 3 |
%  -----------
% | 4 | 5 | 6 |
%  -----------
% 
% 
% possible pairs for varargin:
% 
% 'OuterMargin', [left bottom right top]
% 'om', [left bottom right top]
%    space to leave around left, bottom, right, top of full 
%    plot area, in normalized units
%    by default, the values are set to [0.05 0.05 0.05 0.05]
%
% 'InnerMargin', [horiz vert]
% 'im', [horiz vert]
%    horiz is the space between plots horizontally (width)
%    vert is the space between plots vertically (height)
%    both in normalized units
%    by default, these are set to [0.05 0.05]


% Initialize outermargin and innermargin values
OutM=[0.05 0.05 0.05 0.05];
InM=[0.05 0.05];

% loop through arg pairs
nArgExtra=length(varargin);
if nArgExtra>0
    for i=1:2:nArgExtra
        if strcmpi(varargin{i},'outermargin')
            OutM=varargin{i+1};
        elseif strcmpi(varargin{i},'om')
            OutM=varargin{i+1};
        elseif strcmpi(varargin{i},'innermargin')
            InM=varargin{i+1};
        elseif strcmpi(varargin{i},'im')
            InM=varargin{i+1};
        end
    end
end           

% set figure units to normalized
set(ifig,'units','normalized')

% determine width of plots in normalized units
pwid = ( 1. - OutM(3) - OutM(1) - (nc-1)*InM(1) )/nc;
% determine height of plots in normlized units
phgt = ( 1. - OutM(4) - OutM(2) - (nr-1)*InM(2) )/nr;


% figure out current row and column for panel number
% figure out columns first, since panel numbers increase 
%   left to right first, then top to bottom
CurCol = mod(panel,nc); 
  if CurCol==0; CurCol=nc; end
CurRow = mod( (panel-CurCol)/nc + 1, nr);
  if CurRow==0; CurRow=nr; end

  
% set up position
left = OutM(1) + (CurCol-1)*pwid + (CurCol-1)*InM(1);
bot = 1. - OutM(4) - CurRow*phgt - (CurRow-1)*InM(2);

pos = [left bot pwid phgt];

  
% set up axes (this handle will be returned from the function)
ax = axes( 'Parent', ifig, 'Units', 'Normalized', 'Position', pos );
