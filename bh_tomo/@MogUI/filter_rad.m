function filter_rad(obj,varargin)

fs = 11;
if isa(varargin{1},'matlab.ui.control.UIControl')
    fs = varargin{1}.FontSize;
elseif ispc
    fs = 9;
end

no = obj.handles.listMogs.Value;
if no<=0 || no>length(obj.mogs)
    warndlg('No MOG selected')
    return
end

mog = obj.mogs(no);

width = 900;
height = 600;

% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

str = get_str_locale();

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_spectra',...
    'Name','Filter Traces',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

haxes1 = axes('Units','points','Parent',f);    

hcontrol = uipanel(f,'Title','Control',...
    'Units','points',...
    'FontSize',fs);

hnumberl = uicontrol('Style','text',...
    'String','Trace no',...
    'Units','points',...
    'HorizontalAlignment','right',...
    'FontSize',fs,...
    'Parent',hcontrol);
hnumber = uicontrol('Style','edit',...
    'String','1',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@trace_no,...
    'Parent',hcontrol);
hnext = uicontrol('Style','pushbutton',...
    'String','Next Trace',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@next,...
    'Parent',hcontrol);
hprevious = uicontrol('Style','pushbutton',...
    'String','Previous Trace',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@previous,...
    'Parent',hcontrol);

hdetrend = uicontrol('Style','checkbox',...
    'String','Linear Detrending',...
    'Value', mog.processParams.detrend, ...
    'Units','points',...
    'Callback',@process,...
    'FontSize',fs,...
    'Parent',hcontrol);
hbandpass = uicontrol('Style','checkbox',...
    'String','Apply Bandpass Filter',...
    'Value', mog.processParams.bandpass, ...
    'Units','points',...
    'Callback',@process,...
    'FontSize',fs,...
    'Parent',hcontrol);
hlowcutl = uicontrol('Style','text',...
    'String','Low-Cut Frequency',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',hcontrol);
hlowcut = uicontrol('Style','edit',...
    'String',num2str(mog.processParams.lowcut),...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@process,...
    'Parent',hcontrol);
hhighcutl = uicontrol('Style','text',...
    'String','High-Cut Frequency',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',hcontrol);
hhighcut = uicontrol('Style','edit',...
    'String',num2str(mog.processParams.highcut),...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@process,...
    'Parent',hcontrol);

hstore = uicontrol('Style','pushbutton',...
    'String','Store Filtered Traces',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@store,...
    'Parent',hcontrol);
hreset = uicontrol('Style','pushbutton',...
    'String','Reset Traces',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@reset,...
    'Parent',hcontrol);

filt_traces = mog.traces;
update_trace()

f.Visible = 'on';

    function resizeUI(varargin)
        f.Visible = 'off';
        
        hBorder = 15;
        width = f.Position(3);
        height = f.Position(4);
        
        axBorder = haxes1.Position(2);
        inWidth = width-axBorder-2*hBorder;
        inHeight = height - axBorder-2*hBorder;
        cWidth = width-2*hBorder;
        cHeight =  0.4*inHeight;
        
        haxes1.Position = [axBorder axBorder+cHeight+hBorder inWidth 0.6*inHeight-hBorder];
        hcontrol.Position = [hBorder hBorder cWidth cHeight];
        
        vSize = hnext.Position(4);
        hSize = (cWidth-7*hBorder)/6;
        
        hnumberl.Position = [hBorder cHeight-vSize-hBorder hSize vSize];
        hnumber.Position = [2*hBorder+hSize cHeight-vSize-hBorder hSize vSize];
        hprevious.Position = [hBorder cHeight-2*vSize-2*hBorder hSize vSize];
        hnext.Position = [2*hBorder+hSize cHeight-2*vSize-2*hBorder hSize vSize];

        hdetrend.Position  = [4*hBorder+3*hSize cHeight-vSize-hBorder 2*hSize+hBorder vSize];
        hbandpass.Position = [4*hBorder+3*hSize cHeight-2*vSize-2*hBorder 2*hSize+hBorder vSize];
        hlowcutl.Position  = [4*hBorder+3*hSize cHeight-3*vSize-3*hBorder hSize vSize];
        hlowcut.Position   = [5*hBorder+4*hSize cHeight-3*vSize-3*hBorder hSize vSize];        
        hhighcutl.Position = [4*hBorder+3*hSize cHeight-4*vSize-4*hBorder hSize vSize];
        hhighcut.Position  = [5*hBorder+4*hSize cHeight-4*vSize-4*hBorder hSize vSize];        
        
        hstore.Position = [6*hBorder+5*hSize hBorder hSize vSize];        
        hreset.Position = [5*hBorder+4*hSize hBorder hSize vSize];        
        
        f.Visible = 'on';
    end
    function closeWindow(varargin)
        delete(f)
    end
    function trace_no(varargin)
        trc_nb = str2double(hnumber.String);
        if trc_nb < 1
            trc_nb = 1;
        end
        if trc_nb > mog.data.ntrace
            trc_nb = mog.data.ntrace;
        end
        hnumber.String = num2str(trc_nb);
        update_trace()
    end
    function next(varargin)
        trc_nb = str2double(hnumber.String)+1;
        if trc_nb < 1
            trc_nb = 1;
        end
        if trc_nb > mog.data.ntrace
            trc_nb = mog.data.ntrace;
        end
        hnumber.String = num2str(trc_nb);
        update_trace()
    end
    function previous(varargin)
        trc_nb = str2double(hnumber.String)-1;
        if trc_nb < 1
            trc_nb = 1;
        end
        if trc_nb > mog.data.ntrace
            trc_nb = mog.data.ntrace;
        end
        hnumber.String = num2str(trc_nb);
        update_trace()
    end
    function process(varargin)
        traces = mog.data.rdata;
        if hdetrend.Value == 1
            traces = detrend_rad(mog.data.rdata);
        end
        if hbandpass.Value == 1
            dt = mog.data.timec * mog.fac_dt;

            flc = 1e-3*str2double(hlowcut.String);
            fhc = 1e-3*str2double(hhighcut.String);
            Fs = 1/dt;
            HalfFs = Fs/2;
            Wp = [flc/HalfFs fhc/HalfFs];
            Ws = [0.7 * Wp(1) 1.5*Wp(2)];
            Rp = 3; Rs = 40;
            [n,Wn] = cheb1ord(Wp,Ws,Rp,Rs);

            [b,a] = cheby1(n,0.5,Wn);

            h=waitbar(0,'Filtering ...');
            for n=1:mog.data.ntrace
                filt_traces(:,n) = filtfilt(b, a, traces(:,n));
                waitbar(n/mog.data.ntrace,h)
            end
            close(h)
        else
            filt_traces = traces;
        end
        update_trace()
    end
    function update_trace(varargin)
        trc_nb = str2double(hnumber.String);
        trace = mog.data.rdata(:,trc_nb);
        ftrace = filt_traces(:,trc_nb);
        snr1 = computeSNR(trace);
        snr2 =  computeSNR(ftrace);
        
        if hbandpass.Value==1 || hdetrend.Value==1
            c = [0.8 0.8 0.8];
            plot(haxes1, mog.data.timestp, trace, 'color', c)
            hold on
            c = [0 0.4470 0.7410];
            plot(haxes1, mog.data.timestp, ftrace, 'color', c)
            title(['S/N (orig): ',num2str(snr1), 'S/N (filt): ',num2str(snr2)])
            hold off
            legend(haxes1, 'Original', 'Filtered')
        else
            c = [0 0.4470 0.7410];
            plot(haxes1, mog.data.timestp, trace, 'color', c)
            title(['S/N: ',num2str(snr1)])
        end
        xlabel(haxes1, [str.s22,' [',mog.data.tunits,']'])
        ylabel(haxes1, str.s21)
    end
    function store(varargin)
        mog.processParams.detrend = hdetrend.Value;
        mog.processParams.bandpass = hbandpass.Value;
        mog.processParams.lowcut = str2double(hlowcut.String);
        mog.processParams.highcut = str2double(hhighcut.String);
        mog.traces = filt_traces;
        obj.notify('mogEdited')
    end
    function reset(varargin)
        hdetrend.Value = 0;
        hbandpass.Value = 0;
        mog.processParams.detrend = hdetrend.Value;
        mog.processParams.bandpass = hbandpass.Value;
        mog.traces = mog.data.rdata;
        update_trace()
        obj.notify('mogEdited')
    end
end