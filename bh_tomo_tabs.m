function bh_tomo_tabs(hObject, eventdata, handles, varargin)
%created  2012-10-11  YH
%Last modified 2012-  YH

%Create tabs

handles.hTabGroup = uitabgroup('Parent',handles.bh_tomo,'Position',[.1 .99 .9 .005]);drawnow; 
handles.tabs.home = uitab('parent',handles.hTabGroup, 'title','      Home     ');
handles.tabs.db = uitab('parent',handles.hTabGroup, 'title','      Data base     ');
handles.tabs.pick = uitab('parent',handles.hTabGroup, 'title','      Semi-automatic picking (x-corr)      ');
handles.tabs.phasePick = uitab('parent',handles.hTabGroup, 'title','      Automatic picking (AIC-CWT)      ');
handles.tabs.tt = uitab('parent',handles.hTabGroup, 'title','      Process travel time      ');
handles.tabs.amp = uitab('parent',handles.hTabGroup, 'title','     Process amplitudes     ');
handles.tabs.fitCovar = uitab('parent',handles.hTabGroup, 'title','      Covariance model      ');
handles.tabs.inv = uitab('parent',handles.hTabGroup, 'title','      Inversion      ');
handles.tabs.interp = uitab('parent',handles.hTabGroup, 'title','      Interpretation      ');
handles.tabs.timelaps = uitab('parent',handles.hTabGroup, 'title','      Time-lapse data processing      ');
%handles.hpTab1 = uipanel('Parent',handles.tab1,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab2 = uipanel('Parent',handles.tab2,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab3 = uipanel('Parent',handles.tab3,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab4 = uipanel('Parent',handles.tab4,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab5 = uipanel('Parent',handles.tab5,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab6 = uipanel('Parent',handles.tab6,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab7 = uipanel('Parent',handles.tab7,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab8 = uipanel('Parent',handles.tab8,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
%handles.hpTab9 = uipanel('Parent',handles.tab9,'Title','','FontSize',8,...
%              'Tag', 'tabPanel','Position',[.01 .01 .98 .98]);
          
%TabsVisibleOff(handles);
%guidata(handles.bh_tomo, handles);
guidata(hObject, handles);
set(handles.hTabGroup,'SelectedIndex',1,'SelectionChangeFcn',...  
    @(hObject, eventdata) tab_selection_change_Callback(hObject, eventdata, handles));

handles.hpTab1 = findobj(handles.tabs.home,'Tag', 'tabPanelHome'); 
if isempty(handles.hpTab1)
    handles.hpTab1 = uipanel('Parent',handles.tabs.home,'Title','','FontSize',8,...
          'Tag', 'tabPanelHome','Position',[.01 .01 .98 .98]);    
    HomePanel(hObject, eventdata, handles);
end
%guidata(handles.bh_tomo, handles);
guidata(hObject, handles);
%set(handles.hpTab1, 'VIsible', 'on');
% set (handles.Menu_File_Import_Data, 'Enable', 'on');
% set (handles.Menu_File_Input_Data, 'Enable', 'on');
% set (handles.Menu_File_Save_Project_As, 'Enable', 'on');
% set (handles.Menu_File_Save_Project, 'Enable', 'on');
%guidata(hObject, handles);   


