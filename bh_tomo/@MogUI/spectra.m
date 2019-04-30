function spectra(obj,varargin)

if exist('pwelch','file')~=2
    warndlg('You need the Signal Processing toolbox to get MOG Spectra')
    return
end
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


f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_spectra',...
    'Name','MOG Spectra',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

haxes1 = axes('Units','points','Parent',f);
haxes2 = axes('Units','points','Parent',f);
haxes3 = axes('Units','points','Parent',f);


Tx = unique(mog.data.Tx_z);
rdata = mog.traces;
Amax = kron(max(abs(rdata)), ones(size(rdata,1),1));
rdata = rdata./Amax;


% find number of Tx positions
uTx = unique([mog.data.Tx_x' mog.data.Tx_y' mog.data.Tx_z'],'rows');
text = cell(size(uTx,1),1);
for nn=1:length(text)
    text{nn} = num2str(nn);
end
htxNumber = uicontrol('Style','text',...
    'String','Tx Number',...
    'HorizontalAlignment','center',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',f);
htxNumberPopup = uicontrol('Style','popupmenu',...
    'String',text,...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',f);
htxElev = uicontrol('Style','text',...
    'String','Tx Elevation :',...
    'HorizontalAlignment','center',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',f);

hpsd = uicontrol('Style','text',...
    'String','PSD Estimation Method',...
    'HorizontalAlignment','center',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',f);
hpsdPopup = uicontrol('Style','popupmenu',...
    'String',{'Welch','Burg - order 2','Burg - order 3','Burg - order 4'},...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',f);

hlowpass = uicontrol('Style','checkbox',...
    'String','Apply Low Pass Filter',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',f);

hfmax = uicontrol('Style','text',...
    'String','F Max',...
    'HorizontalAlignment','center',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',f);
hfmaxEdit = uicontrol('Style','edit',...
    'String','400',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',f);

hsnrscale = uicontrol('Style','text',...
    'String','SNR Scale',...
    'HorizontalAlignment','center',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',f);
hsnrscalePopup = uicontrol('Style','popupmenu',...
    'String',{'Linear','Log'},...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@snrscale,...
    'Parent',f);

hdfreq = uipanel(f,'Title','Dominant Frequency',...
    'Units','points',...
    'FontSize',fs);

hdfreqCheck = uicontrol('Style','checkbox',...
    'String','Compute & Show',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',hdfreq);
hdfmin = uicontrol('Style','text',...
    'String','F min',...
    'HorizontalAlignment','right',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',hdfreq);
hdfminEdit = uicontrol('Style','edit',...
    'String','0',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',hdfreq);
hdfmax = uicontrol('Style','text',...
    'String','F max',...
    'HorizontalAlignment','right',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',hdfreq);
hdfmaxEdit = uicontrol('Style','edit',...
    'String','400',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFigs,...
    'Parent',hdfreq);

updateFigs()

f.Visible = 'on';

    function resizeUI(varargin)
        f.Visible = 'off';
        
        hBorder = 15;
        width = f.Position(3);
        height = f.Position(4);
        
        axBorder = haxes1.Position(2);
        
        axWidth = 9/40*width;
        haxes1.Position = [axBorder axBorder axWidth height-2.5*axBorder];
        haxes2.Position = [1.5*axBorder+axWidth axBorder axWidth height-2.5*axBorder];
        haxes3.Position = [2*axBorder+2*axWidth axBorder 1/6*width height-2.5*axBorder];
        
        hSize = width/22;
        
        vFac = 0.8*height/600;
        if vFac<1
            vFac = 1;
        end
        vSize = 22*vFac;
        vSpace = 5*vFac;
        
        htxNumber.Position = [width-hBorder-3*hSize height-axBorder 2*hSize vSize];
        htxNumberPopup.Position = [width-hBorder-3*hSize height-axBorder-vSize 2*hSize vSize];
        htxElev.Position = [width-hBorder-3.5*hSize height-axBorder-2*vSize-vSpace 3*hSize vSize];
        
        hpsd.Position = [width-hBorder-3.5*hSize height-axBorder-3*vSize-4*vSpace 3*hSize vSize];
        hpsdPopup.Position = [width-hBorder-3*hSize height-axBorder-4*vSize-5*vSpace 2*hSize vSize];
        
        hlowpass.Position = [width-hBorder-4*hSize height-axBorder-5*vSize-8*vSpace 4*hSize vSize];
        
        hfmax.Position = [width-hBorder-3.5*hSize height-axBorder-6*vSize-11*vSpace 3*hSize vSize];
        hfmaxEdit.Position = [width-hBorder-3*hSize height-axBorder-7*vSize-12*vSpace 2*hSize vSize];
        hsnrscale.Position = [width-hBorder-3.5*hSize height-axBorder-8*vSize-15*vSpace 3*hSize vSize];
        hsnrscalePopup.Position = [width-hBorder-3*hSize height-axBorder-9*vSize-16*vSpace 2*hSize vSize];
        
        hSize = 4*hSize;
        hdfreq.Position = [width-hBorder-hSize axBorder hSize 6*vSize];
        
        hSize = hSize/4;
        hdfmax.Position = [hBorder hBorder hSize vSize];
        hdfmin.Position = [hBorder hBorder+vSize+vSpace hSize vSize];
        hSize = 1.2*hSize;
        hdfmaxEdit.Position = [hBorder+hSize hBorder hSize vSize];
        hdfminEdit.Position = [hBorder+hSize hBorder+vSize+vSpace hSize vSize];
        
        hdfreqCheck.Position = [hBorder hBorder+2*vSize+3*vSpace 3*hSize vSize];
        f.Visible = 'on';
    end
    function closeWindow(varargin)
        delete(f)
    end
    function snrscale(varargin)
        if hsnrscalePopup.Value==1
            haxes3.XScale = 'linear';
        else % log
            haxes3.XScale = 'log'; 
        end
    end
    function updateFigs(varargin)
        
        n = htxNumberPopup.Value;
        
        dt = mog.data.timec * mog.fac_dt;
        
        ind =  Tx(n) == mog.data.Tx_z ;
        traces = rdata(:,ind);
        nfft = 2^(1+nextpow2(size(traces,1)));
        
        fac_f = 1.0;
        fac_t = 1.0;
        if strcmp(mog.data.tunits,'ns')
            % radar: freq nominale donnée en MHz
            fac_f = 10^6;
            fac_t = 10^-9;
        elseif strcmp(mog.data.tunits,'ms')
            % sismique: on assume que f dominante est en kHz
            warndlg('Assuming Tx nominal frequency in _kHz_')
            fac_f = 10^3;
            fac_t = 10^-3;
        else
            warndlg({'Assuming Tx nominal frequency in _Hz_',...
                'Assuming time step in _s_'})
        end
        f0 = mog.data.rnomfreq * fac_f;
        dt = dt * fac_t;
        Fs = 1/dt;
        
        % bruit pris sur dernière 20 ns
        win_snr = round(20/mog.data.timec);
        snr = data_select(traces, f0, dt, win_snr);
        
        if hlowpass.Value == 1
            HalfFs = Fs/2;
            Wp = 1.4*f0/HalfFs;
            Ws = 1.6*f0/HalfFs;
            Rp = 3;
            Rs = 40;
            [nc,Wn] = cheb1ord(Wp,Ws,Rp,Rs);
            
            [b,a] = cheby1(nc,0.5,Wn);
            
            %[b,a] = ellip(10, Rp, Rs, Wp);
            
            for nt=1:size(traces,2)
                traces(:,nt) = filtfilt(b,a,traces(:,nt));
            end
        end
        if hdfreqCheck.Value == 1
            [tmp,freq] = pburg(traces(:,1),2,nfft,Fs);
            fdom = zeros(size(traces,2), 1);
            Pxx = zeros(length(tmp), size(traces,2));
            
            min_fdom = str2double(hdfminEdit.String) * fac_f;
            max_fdom = str2double(hdfmaxEdit.String) * fac_f;
            
            for nt=1:size(traces,2)
                [fdom(nt), Pxx(:,nt)] = get_fdom(traces(:,nt), dt, f0, min_fdom, max_fdom);
            end
        else
            method = hpsdPopup.Value;
            if method==1
                [tmp, freq] = pwelch(traces(:,1),[],[],[],Fs);
                Pxx = zeros(length(tmp), size(traces,2));
                Pxx(:,1) = tmp;
                for nt=2:size(traces,2)
                    Pxx(:,nt) = pwelch(traces(:,nt),[],[],[],Fs);
                end
            else
                ordre = method;
                [tmp,freq] = pburg(traces(:,1),ordre,nfft,Fs);
                Pxx = zeros(length(tmp), size(traces,2));
                Pxx(:,1) = tmp;
                for nt=2:size(traces,2)
                    Pxx(:,nt) = pburg(traces(:,nt),ordre,nfft,Fs);
                end
            end
        end
        z = mog.data.Rx_z(ind);
        timestp = mog.data.timestp * mog.fac_dt;
        
        imagesc(timestp, z, traces','Parent',haxes1)
        max(timestp)
        xlabel(haxes1,['Time [',mog.data.tunits,']'],'FontSize',fs)
        ylabel(haxes1,['Rx elevation [',mog.data.cunits,']'],'FontSize',fs)
        title(haxes1,'Normalized amplitudes','FontSize',fs+1)
        
        Pxx(isnan(Pxx)) = 0;
        
        imagesc(freq/fac_f, z, log10(Pxx'),'Parent',haxes2)
        xlabel(haxes2,'Frequency [MHz]','FontSize',fs)
        title(haxes2,'log_{10} Power spectra','FontSize',fs+1)
        if hdfreqCheck.Value == 1
            hold(haxes2,'on')
            ind = fdom~=-1;
            plot(haxes2, fdom(ind)/fac_f,z(ind),'ko')
            hold(haxes2,'off')
        end
        haxes2.XLim = [0 str2double(hfmaxEdit.String)];
        
        plot(haxes3, snr, z)
        xlabel(haxes3,'SNR','FontSize',fs)
        title(haxes3,'Signal-to-Noise ratio','FontSize',fs+1)
        if hsnrscalePopup.Value==1
            haxes3.XScale = 'linear';
        else % log
            haxes3.XScale = 'log';
        end
        haxes3.YLim = haxes1.YLim;
        haxes2.YDir = 'reverse';
        
        htxElev.String = ['Tx elevation: ',num2str(Tx(n)),' m'];
    end
end
