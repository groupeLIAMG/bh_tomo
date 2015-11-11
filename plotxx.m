function [axout,h1,h2]=plotxx(varargin)
%PLOTXX Graphs with x tick labels on top and bottom
%   PLOTXX(X1,Y1,X2,Y2) plots Y1 versus X1 with x-axis labeling
%   on bottom and plots Y2 versus X2 with x-axis labeling on
%   top.
%
%   PLOTXX(X1,Y1,X2,Y2,FUN) uses the plotting function FUN 
%   instead of PLOT to produce each graph. FUN can be a 
%   function handle or a string that is the name of a plotting 
%   function, e.g. plot, semilogx, semilogy, loglog, stem, 
%   etc. or any function that accepts the syntax H = FUN(X,Y).
%   For example
%      PLOTXX(X1,Y1,X2,Y2,@loglog)  % Function handle
%      PLOTXX(X1,Y1,X2,Y2,'loglog') % String
%   
%   PLOTXX(X1,Y1,X2,Y2,FUN1,FUN2) uses FUN1(X1,Y1) to plot the data for
%   the left axes and FUN2(X2,Y2) to plot the data for the right axes.
%
%   PLOTXX(AX,...) plots into AX as the main axes, instead of GCA.  If AX
%   is the vector of axes handles returned by a previous call to PLOTXX,
%   then the secondary axes will be ignored.
%
%   [AX,H1,H2] = PLOTXX(...) returns the handles of the two axes created in
%   AX and the handles of the graphics objects from each plot in H1
%   and H2. AX(1) is the left axes and AX(2) is the right axes.  
%
%   See also PLOT, function_handle

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.18.4.9 $  $Date: 2005/06/21 19:29:44 $

%   Made by Bernard Giroux from the plotyy.m file by MathWorks
%

% Parse possible Axes input
%   First see if we have a vector of associated PLOTXX axes
if length(varargin{1}(:))==2 && ...
        all(ishandle(varargin{1}(:))) && ...        
        isequal(lower(strvcat(get(varargin{1},'type'))),['axes';'axes']) && ...
        any(findall(varargin{1}(1),'tag','PlotyyDeleteProxy')) && ...
        any(findall(varargin{1}(2),'tag','PlotyyDeleteProxy'))
    varargin{1} = varargin{1}(1);
end
[cax,args,nargs] = axescheck(varargin{:});
caxspecified = ~isempty(cax);

cax = newplot(cax);
fig = get(cax,'Parent');

if nargs<4, error('Not enough input arguments.'); end
[x1,y1,x2,y2] = deal(args{1:4});
if nargs<5, fun1 = @plot; else fun1 = args{5}; end
if nargs<6, fun2 = fun1;  else fun2 = args{6}; end

hold_state = ishold(cax);
co = get(fig,'DefaultAxesColorOrder');

% Plot first plot
ax(1) = cax;
set(fig,'NextPlot','add')
[h1,ax(1)] = fevalfun(fun1,ax(1),x1,y1,caxspecified);
set(ax(1),'Box','on','ActivePositionProperty','position')
xlim1 = get(ax(1),'XLim');
ylim1 = get(ax(1),'YLim');


% Try to produce different colored lines on each plot by
% shifting the colororder by the number of different colors
% in the first plot.
colors = get(findobj(h1,'Type','line'),{'Color'});
if ~isempty(colors),
  colors = unique(cat(1,colors{:}),'rows');
  n = size(colors,1);
  % Set the axis color to match the single colored line
  if n==1, set(ax(1),'XColor',colors), end
  set(fig,'DefaultAxesColorOrder',co([n+1:end 1:n],:))
end

% Plot second plot
ax1hv = get(ax(1),'HandleVisibility');
ax(2) = axes('HandleVisibility',ax1hv,'Units',get(ax(1),'Units'), ...
    'Position',get(ax(1),'Position'),'Parent',fig);
[h2,ax(2)] = fevalfun(fun2,ax(2),x2,y2,caxspecified);
set(ax(2),'XAxisLocation','top','Color','none', ...
          'XGrid','off','YGrid','off','Box','off');
xlim2 = get(ax(2),'XLim');
ylim2 = get(ax(2),'YLim');

% Check to see if the second plot also has one colored line
colors = get(findobj(h2,'Type','line'),{'Color'});
if ~isempty(colors),
  colors = unique(cat(1,colors{:}),'rows');
  n = size(colors,1);
  if n==1, set(ax(2),'XColor',colors), end
end

islog1 = strcmp(get(ax(1),'XScale'),'log');
islog2 = strcmp(get(ax(2),'XScale'),'log');

if islog1, xlim1 = log10(xlim1); end
if islog2, xlim2 = log10(xlim2); end

% Find bestscale that produces the same number of y-ticks for both
% the left and the right.
[low,high,ticks] = bestscale(xlim1(1),xlim1(2),xlim2(1),xlim2(2),islog1 | islog2);

if ~isempty(low)
    if islog1,
      xticks1 = logsp(low(1),high(1),ticks(1));
      low(1) = 10.^low(1);
      high(1) = 10.^high(1);
      xlim1 = 10.^xlim1;
    else
      xticks1 = linspace(low(1),high(1),ticks(1));
    end

    if islog2,
      xticks2 = logsp(low(2),high(2),ticks(2));
      low(2) = 10.^low(2);
      high(2) = 10.^high(2);
      xlim2 = 10.^xlim2;
    else
      xticks2 = linspace(low(2),high(2),ticks(2));
    end

    % Set ticks on both plots the same
    set(ax(1),'XLim',[low(1) high(1)],'XTick',xticks1);
    set(ax(2),'XLim',[low(2) high(2)],'XTick',xticks2);
    set(ax,'YLim',[min(ylim1(1),ylim2(1)) max(ylim1(2),ylim2(2))])

    % Set tick labels if axis ticks aren't at decade boundaries
    % when in log mode
    if islog1
      decade1 =  abs(floor(log10(xticks1)) - log10(xticks1));
    end
    if islog2
      decade2 =  abs(floor(log10(xticks2)) - log10(xticks2));
    end
    if islog1 & any(decade1 > 0.1)
      for i=length(xticks1):-1:1
        xtickstr1{i} = sprintf('%3g',xticks1(i));
      end
      set(ax(1),'XTickLabel',xtickstr1)
    end

    if islog2 & any(decade2 > 0.1)
      for i=length(xticks2):-1:1
        xtickstr2{i} = sprintf('%.3g',xticks2(i));
      end
      set(ax(2),'XTickLabel',xtickstr2)
    end
    
else
    % Use the default automatic scales and turn off the box so we
    % don't get double tick marks on each side.  We'll still get
    % the grid from the left axes though (if it is on).
    set(ax,'Box','off')
end

% create DeleteProxy objects (an invisible text object in
% the first axes) so that the other axes will be deleted
% properly.
DeleteProxy(1) = text('Parent',ax(1),'Visible','off',...
                      'Tag','PlotyyDeleteProxy',...
                      'HandleVisibility','off',...
        'DeleteFcn','try;delete(get(gcbo,''userdata''));end');
DeleteProxy(2) = text('Parent',ax(2),'Visible','off',...
                       'Tag','PlotyyDeleteProxy',...
                       'HandleVisibility','off',...
        'DeleteFcn','try;delete(get(gcbo,''userdata''));end');
set(DeleteProxy(1),'UserData',ax(2));
set(DeleteProxy(2),'UserData',DeleteProxy(1));

if ~hold_state, hold(cax,'off'), end

% Reset colororder
set(fig,'DefaultAxesColorOrder',co,'CurrentAxes',ax(1))

if nargout>0, axout = ax; end



%---------------------------------------------------
function [low,high,ticks] = bestscale(umin,umax,vmin,vmax,islog)
%BESTSCALE Returns parameters for "best" yy scale.

penalty = 0.02;

% Determine the good scales
[ulow,uhigh,uticks] = goodscales(umin,umax);
[vlow,vhigh,vticks] = goodscales(vmin,vmax);

% Find good scales where the number of ticks match
[u,v] = meshgrid(uticks,vticks);
[j,i] = find(u==v);

if islog % Filter out the cases where power of ten's don't match
  for k=length(i):-1:1
    utest = logsp(ulow(i(k)),uhigh(i(k)),uticks(i(k)));
    vtest = logsp(vlow(j(k)),vhigh(j(k)),vticks(j(k)));
    upot = abs(log10(utest)-round(log10(utest))) < 10*eps*log10(utest);
    vpot = abs(log10(vtest)-round(log10(vtest))) < 10*eps*log10(vtest);
    if ~isequal(upot,vpot),
       i(k) = [];
       j(k) = [];
    end
  end
end

if ~isempty(i)
  udelta = umax-umin;
  vdelta = vmax-vmin;
  ufit = ((uhigh(i)-ulow(i)) - udelta)./(uhigh(i)-ulow(i));
  vfit = ((vhigh(j)-vlow(j)) - vdelta)./(vhigh(j)-vlow(j));

  fit = ufit + vfit + penalty*(max(uticks(i)-6,1)).^2;

  % Choose base fit
  k = find(fit == min(fit)); k=k(1);
  low = [ulow(i(k)) vlow(j(k))];
  high = [uhigh(i(k)) vhigh(j(k))];
  ticks = [uticks(i(k)) vticks(j(k))];
else
  % Return empty to signal calling routine that we weren't able to
  % find matching scales.
  low = [];
  high = [];
  ticks = [];
end



%------------------------------------------------------------
function [low,high,ticks] = goodscales(xmin,xmax)
%GOODSCALES Returns parameters for "good" scales.
%
% [LOW,HIGH,TICKS] = GOODSCALES(XMIN,XMAX) returns lower and upper
% axis limits (LOW and HIGH) that span the interval (XMIN,XMAX) 
% with "nice" tick spacing.  The number of major axis ticks is 
% also returned in TICKS.

BestDelta = [ .1 .2 .5 1 2 5 10 20 50 ];
penalty = 0.02;

% Compute xmin, xmax if matrices passed.
if length(xmin) > 1, xmin = min(xmin(:)); end
if length(xmax) > 1, xmax = max(xmax(:)); end
if xmin==xmax, low=xmin; high=xmax+1; ticks = 1; return, end

% Compute fit error including penalty on too many ticks
Xdelta = xmax-xmin;
delta = 10.^(round(log10(Xdelta)-1))*BestDelta;
high = delta.*ceil(xmax./delta);
low = delta.*floor(xmin./delta);
ticks = round((high-low)./delta)+1;



%---------------------------------------------
function  y = logsp(low,high,n)
%LOGSP Generate nice ticks for log plots
%   LOGSP produces linear ramps between 10^k values.

y = linspace(low,high,n);

k = find(abs(y-round(y))<=10*eps*max(y));
dk = diff(k);
p = find(dk > 1);

y = 10.^y;

for i=1:length(p)
  r = linspace(0,1,dk(p(i))+1)*y(k(p(i)+1));
  y(k(p(i))+1:k(p(i)+1)-1) = r(2:end-1);
end



%---------------------------------------------
function [h,ax] = fevalfun(func,ax,x,y,caxspecified)
%FEVALFUN Evaluate the function with or without a specified axes
%    FEVALFUN returns an appropriate error if PLOTYY cannot determine
%    the syntax for calling it with or without a specified axes.

if ~caxspecified
    % If no axes input was specified, then assume we are dealing with
    %   a handlevisibility/on axes returned by NEWPLOT
    h = feval(func,x,y);
else
    % If an axes input was specified, then first try to call func using
    %   the axes input as the first input argument
    try
        h = feval(func,ax,x,y);
    catch
        % If func won't accept an axes a the first input arguemnt, try
        %   to call func with a parent/axes pair
        try
            h = feval(func,x,y,'Parent',ax);
        catch
            % Make a string out of the function name
            if ~ischar(func)
                funcstr = func2str(func);
            else
                funcstr = func;
            end
            % Finally throw the error
            warning(['Function ''%s'' might not allow the common syntax ' ...
                'for parent\n         axes specification.  The axes ' ...
                'handle either must be the first\n         input argument, ' ...
                'or must be specifed as a property/value pair.\n         ' ...
                'See the following error for more information.\n'], ...
                funcstr)
            rethrow(lasterror)
        end
    end
end
if ~isempty(h)   %%%YH
ax = get(h(1),'Parent');
end
