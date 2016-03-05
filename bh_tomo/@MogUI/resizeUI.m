function resizeUI(obj,varargin)

oldUnits = obj.handles.hp.Units;
obj.handles.hp.Units = 'points';
obj.handles.hp.Visible = 'off';

p = obj.handles.hp.Position;
width = p(3);  % prefered: 630
height = p(4); % prefered: 380

hSize = width/7;
hSpace = width/120;
hBorder = width/30;

vFac = 0.8*height/360;
if vFac<1
    vFac = 1;
end
vSize = 22*vFac;
vSpace = 5*vFac;
vBorderTop = 45*vFac;
vBorder = 15*vFac;

obj.handles.addMog.Position = [hBorder height-vBorderTop 1.5*hSize vSize];
obj.handles.removeMog.Position = [hBorder+2*hSpace+1.5*hSize height-vBorderTop 1.5*hSize vSize];

obj.handles.listMogs.Position = [1.5*hBorder vBorder+5*vSpace+5.5*vSize 2*hSpace+3*hSize-hBorder height-2*vBorder-vBorderTop-6.5*vSize];

obj.handles.textType.Position = [hBorder vBorder+4*vSpace+4*vSize 1.25*hSize vSize];
obj.handles.textTx.Position = [hBorder vBorder+3*vSpace+3*vSize 1.25*hSize vSize];
obj.handles.textRx.Position = [hBorder vBorder+2*vSpace+2*vSize 1.25*hSize vSize];
obj.handles.popupType.Position = [hBorder+2*hSpace+1.5*hSize vBorder+4*vSpace+4*vSize 1.5*hSize vSize];
obj.handles.popupTx.Position = [hBorder+2*hSpace+1.5*hSize vBorder+3*vSpace+3*vSize 1.5*hSize vSize];
obj.handles.popupRx.Position = [hBorder+2*hSpace+1.5*hSize vBorder+2*vSpace+2*vSize 1.5*hSize vSize];

obj.handles.airBefore.Position = [hBorder vBorder+vSpace+vSize 1.5*hSize vSize];
obj.handles.airAfter.Position = [hBorder vBorder 1.5*hSize vSize];
obj.handles.textAirBefore.Position = [hBorder+2*hSpace+1.5*hSize vBorder+vSpace+vSize 1.5*hSize vSize];
obj.handles.textAirAfter.Position = [hBorder+2*hSpace+1.5*hSize vBorder 1.5*hSize vSize];

obj.handles.renameMog.Position = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-vSize hSize vSize];
obj.handles.importMog.Position = [width-hBorder-2*hSize-hSpace   height-vBorderTop-vSize hSize vSize];
obj.handles.rawData.Position   = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-vSpace-2*vSize hSize vSize];
obj.handles.zop.Position       = [width-hBorder-2*hSize-hSpace   height-vBorderTop-vSpace-2*vSize hSize vSize];
obj.handles.spectra.Position   = [width-hBorder-hSize            height-vBorderTop-vSpace-2*vSize hSize vSize];
obj.handles.statsTt.Position   = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-2*vSpace-3*vSize hSize vSize];
obj.handles.statsAmp.Position  = [width-hBorder-2*hSize-hSpace   height-vBorderTop-2*vSpace-3*vSize hSize vSize];
obj.handles.rays.Position      = [width-hBorder-hSize            height-vBorderTop-2*vSpace-3*vSize hSize vSize];
obj.handles.exportTt.Position  = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-3*vSpace-4*vSize hSize vSize];
obj.handles.exportTau.Position = [width-hBorder-2*hSize-hSpace   height-vBorderTop-3*vSpace-4*vSize hSize vSize];
obj.handles.prune.Position     = [width-hBorder-hSize            height-vBorderTop-3*vSpace-4*vSize hSize vSize];

obj.handles.textDate.Position = [width-hBorder-2.5*hSize-hSpace vBorder hSize vSize];
obj.handles.editDate.Position = [width-hBorder-1.5*hSize vBorder 1.5*hSize vSize];

obj.handles.textFreq.Position = [width-hBorder-3*hSize-hSpace vBorderTop+5*(vSpace+vSize) 2*hSize vSize];
obj.handles.editFreq.Position = [width-hBorder-hSize vBorderTop+5*(vSpace+vSize) hSize vSize];
obj.handles.textFeedRx.Position = [width-hBorder-3*hSize-hSpace vBorderTop+4*(vSpace+vSize) 2*hSize vSize];
obj.handles.editFeedRx.Position = [width-hBorder-hSize vBorderTop+4*(vSpace+vSize) hSize vSize];
obj.handles.textFeedTx.Position = [width-hBorder-3*hSize-hSpace vBorderTop+3*(vSpace+vSize) 2*hSize vSize];
obj.handles.editFeedTx.Position = [width-hBorder-hSize vBorderTop+3*(vSpace+vSize) hSize vSize];
ext = obj.handles.checkDtCorr.Extent;
obj.handles.checkDtCorr.Position = [width-hBorder-1.15*ext(3)-hSize-hSpace vBorderTop+2*(vSpace+vSize) 1.15*ext(3) vSize];
obj.handles.editDtCorr.Position = [width-hBorder-hSize vBorderTop+2*(vSpace+vSize) hSize vSize];
obj.handles.textMultFac.Position = [width-hBorder-3*hSize-hSpace vBorderTop+vSpace+vSize 2*hSize vSize];
obj.handles.editMultFac.Position = [width-hBorder-hSize vBorderTop+vSpace+vSize hSize vSize];

obj.handles.hp.Visible = 'on';
obj.handles.hp.Units = oldUnits;
end
