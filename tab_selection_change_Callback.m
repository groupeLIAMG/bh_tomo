function tab_selection_change_Callback(hObject, eventdata,handles)
%Created 2012-10-11  YH
%Last modified 2012-10-  YH
%Create tabpanels
%Another way to switch between the tabpanels
%handles = guidata(handles.figure1);

tabs = get(hObject,'Children');    
oldTab = eventdata.OldValue;              
oldName = get(tabs(oldTab),'title');         
newTab = eventdata.NewValue;               
newName = get(tabs(newTab),'title');         
hpNewTab = tabs(newTab);
if newTab == 2
    handles.hpTab2 = findobj(handles.tabs.db,'Tag', 'tabPanelDB'); 
    if isempty(handles.hpTab2)
        handles.hpTab2 = uipanel('Parent',handles.tabs.db,'Title','','FontSize',8,...
              'Tag', 'tabPanelDB','Position',[.01 .01 .98 .98]);    
      	guidata( handles.bh_tomo, handles);
        db_file = getappdata(handles.bh_tomo,'db_file');
        bh_tomo_db( 'UserData', db_file );
    end
end
if isempty(Data)
    questdlgData(hObject, eventdata, handles);
else
    switch strtrim(newName)
        case 'Data'
            handles.hpTab2 = findobj(handles.tab2,'Tag', 'tabPanelData');   
            if (isempty(handles.hpTab2) || dataChange(2) == 1) %importData == 1)
                handles.hpTab2 = uipanel('Parent',handles.tab2,'Title','','FontSize',8,...
                          'Tag', 'tabPanelData','Position',[.01 .01 .98 .98]); 
                      guidata( handles.figure1, handles);
                DataPanel( hObject,eventdata,handles,handles.hpTab2);
                dataChange(2)=0;
            %    importData = 0;
            end
        case 'Depth functions'  
            handles.hpTab3 = findobj(handles.tab3,'Tag', 'tabPanelDepth');   
            if (isempty(handles.hpTab3) || dataChange(3) == 1) %importData == 1)
                handles.hpTab3 = uipanel('Parent',handles.tab3,'Title','','FontSize',8,...
                      'Tag', 'tabPanelDepth','Position',[.01 .01 .98 .98]); 
                  guidata( handles.figure1, handles);
                CalculationPanel( hObject,eventdata,handles,handles.hpTab3);
        %        importData = 0;
                dataChange(3) = 0;
            end
        case 'Visualization' 
            handles.hpTab4 = findobj(handles.tab4,'Tag', 'tabPanelVisualization');   
            if (isempty(handles.hpTab4) || dataChange(4) == 1) %importData == 1)
                handles.hpTab4 = uipanel('Parent',handles.tab4,'Title','','FontSize',8,...
                      'Tag', 'tabPanelVisualization','Position',[.01 .01 .98 .98]);  
                  guidata( handles.figure1, handles);
                VisualizationPanel( hObject,eventdata,handles);
                dataChange(4) = 0;
        %        importData = 0;
            end
        case 'Location' 
            handles.hpTab5 = findobj(handles.tab5,'Tag', 'tabPanelLocation');   
            if (isempty(handles.hpTab5) || dataChange(5) == 1) %importData == 1)
                handles.hpTab5 = uipanel('Parent',handles.tab5,'Title','','FontSize',8,...
                      'Tag', 'tabPanelLocation','Position',[.01 .01 .98 .98]); 
                  guidata( handles.figure1, handles);
                LocationParamPanel( hObject, eventdata, handles);
                dataChange(5) = 0;
        %        importData = 0;
            end
        case 'Scale' 
            handles.hpTab6 = findobj(handles.tab6,'Tag', 'tabPanelScale');   
            if (isempty(handles.hpTab6) || dataChange(6) == 1) %importData == 1)
                handles.hpTab6 = uipanel('Parent',handles.tab6,'Title','','FontSize',8,...
                      'Tag', 'tabPanelScale','Position',[.01 .01 .98 .98]);  
                  guidata( handles.figure1, handles);
                ScaleParamPanel( hObject, eventdata, handles);
                dataChange(6) = 0;
        %        importData = 0;
            end
        case 'Skewness'
            handles.hpTab7 = findobj(handles.tab7,'Tag', 'tabPanelSkewness');   
            if (isempty(handles.hpTab7) || dataChange(7) == 1) %importData == 1)
                handles.hpTab7 = uipanel('Parent',handles.tab7,'Title','','FontSize',8,...
                      'Tag', 'tabPanelSkewness','Position',[.01 .01 .98 .98]); 
                  guidata( handles.figure1, handles);
                SkewnessPanel( hObject, eventdata, handles);
                dataChange(7) = 0;
        %        importData = 0;
            end
        case 'Kurtosis' 
            handles.hpTab8 = findobj(handles.tab8,'Tag', 'tabPanelKurtosis');   
            if (isempty(handles.hpTab8) || dataChange(8) == 1) %importData == 1)
                handles.hpTab8 = uipanel('Parent',handles.tab8,'Title','','FontSize',8,...
                      'Tag', 'tabPanelKurtosis','Position',[.01 .01 .98 .98]);
                  guidata( handles.figure1, handles);
                KurtosisPanel( hObject, eventdata, handles);
                dataChange(8) = 0;
        %        importData = 0;
            end
        case 'Outliers' 
            handles.hpTab9 = findobj(handles.tab9,'Tag', 'tabPanelOutliers');   
            if (isempty(handles.hpTab9) || dataChange(9) == 1) %importData == 1)
                handles.hpTab9 = uipanel('Parent',handles.tab9,'Title','','FontSize',8,...
                      'Tag', 'tabPanelOutliers','Position',[.01 .01 .98 .98]); 
                  guidata( handles.figure1, handles);
                OutlierPanel( hObject,eventdata,handles,handles.hpTab9 );
                dataChange(9) = 0;
        %        importData = 0;
            end
    end
end
