function varargout = test_button_ok(varargin)
% TEST_BUTTON_OK MATLAB code for test_button_ok.fig
%      TEST_BUTTON_OK, by itself, creates a new TEST_BUTTON_OK or raises the existing
%      singleton*.
%
%      H = TEST_BUTTON_OK returns the handle to a new TEST_BUTTON_OK or the handle to
%      the existing singleton*.
%
%      TEST_BUTTON_OK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_BUTTON_OK.M with the given input arguments.
%
%      TEST_BUTTON_OK('Property','Value',...) creates a new TEST_BUTTON_OK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_button_ok_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_button_ok_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_button_ok

% Last Modified by GUIDE v2.5 23-May-2013 17:46:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_button_ok_OpeningFcn, ...
                   'gui_OutputFcn',  @test_button_ok_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before test_button_ok is made visible.
function test_button_ok_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test_button_ok (see VARARGIN)

   
% Choose default command line output for test_button_ok
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test_button_ok wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_button_ok_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
currChar = get(hObject,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       delete(hObject);
   end