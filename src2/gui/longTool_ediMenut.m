function varargout = longTool_ediMenut(varargin)
% LONGTOOL_EDIMENUT MATLAB code for longTool_ediMenut.fig
%      LONGTOOL_EDIMENUT, by itself, creates a new LONGTOOL_EDIMENUT or raises the existing
%      singleton*.
%
%      H = LONGTOOL_EDIMENUT returns the handle to a new LONGTOOL_EDIMENUT or the handle to
%      the existing singleton*.
%
%      LONGTOOL_EDIMENUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LONGTOOL_EDIMENUT.M with the given input arguments.
%
%      LONGTOOL_EDIMENUT('Property','Value',...) creates a new LONGTOOL_EDIMENUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before longTool_ediMenut_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to longTool_ediMenut_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help longTool_ediMenut

% Last Modified by GUIDE v2.5 18-Jan-2016 15:24:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @longTool_ediMenut_OpeningFcn, ...
                   'gui_OutputFcn',  @longTool_ediMenut_OutputFcn, ...
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


% --- Executes just before longTool_ediMenut is made visible.
function longTool_ediMenut_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to longTool_ediMenut (see VARARGIN)

% Choose default command line output for longTool_ediMenut
set(hObject, 'Name', 'Lifting Surfaces Menu');
ScreenUnits=get(0,'Units');
set(0,'Units','pixels');
ScreenSize=get(0,'ScreenSize');
set(0,'Units',ScreenUnits);
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
FigPos(2)=1/2*(ScreenSize(4)-FigHeight);
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

handles.output = varargin{1};
handles.model = handles.output;
handles.dir     = varargin{2};

handles = button_layout(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes longTool_ediMenut wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = longTool_ediMenut_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);

function handles = button_layout(handles)
switch handles.model.conf
    case 'convencional'
        set(handles.canardButton,'Visible','off');
        set(handles.horizontalButton,'Visible','on');  
        
    case 'canard'
        set(handles.horizontalButton,'Visible','off');
        set(handles.canardButton,'Visible','on');

    case 'convencional_canard'
        set(handles.horizontalButton,'Visible','on');
        set(handles.canardButton,'Visible','on');
    
    case 'flyWing'
        set(handles.horizontalButton,'Visible','off');
        set(handles.canardButton,'Visible','off');
end


% --- Executes on button press in wingButton.
function wingButton_Callback(hObject, eventdata, handles)
% hObject    handle to wingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.model.ala,~] = wingProperties(handles.model,handles.dir,0);
handles.output = handles.model;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in horizontalButton.
function horizontalButton_Callback(hObject, eventdata, handles)
% hObject    handle to horizontalButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.model.horizontal,~] = hor_canProperties(handles.model, 'horizontal', handles.dir, 0);
handles.output = handles.model;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in canardButton.
function canardButton_Callback(hObject, eventdata, handles)
% hObject    handle to canardButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.model.canard, ~] = hor_canProperties(handles.model, 'canard', handles.dir, 0);
handles.output = handles.model;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in closeButton.
function closeButton_Callback(hObject, eventdata, handles)
% hObject    handle to closeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1,'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    guidata(hObject, handles);
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end
