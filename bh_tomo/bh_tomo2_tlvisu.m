function bh_tomo2_tlvisu( varargin )
%BH_TOMO2_TLVISU Visualize time-lapse tomographic inversion results

rep='';
file='';
if nargin>=2
    rep=varargin{1};
    file=varargin{2};
end

cminAmp = 1;
cmaxAmp = 3;
cminTT = 0.06;
cmaxTT = 0.12;

fs = 11;
if nargin>=3
    fs = varargin{3};
elseif ispc
    fs = 9;
end
vScale = 1;
if ispc
    vScale = 0.81;
end

width = 1400;
height = 875*vScale;
% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_tlinv',...
    'Name','bh_tomo_tlinv',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);




f.Visible = 'on';


    function resizeUI(varargin)

    end

    function closeWindow(varargin)
        quitUI()
    end
    function quitUI(varargin)
        delete(f)
    end


end