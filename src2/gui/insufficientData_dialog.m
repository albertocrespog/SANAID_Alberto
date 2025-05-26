function varargout = insufficientData_dialog(varargin)
% INSUFFICIENTDATA_DIALOG MATLAB code for insufficientData_dialog.fig
%      INSUFFICIENTDATA_DIALOG, by itself, creates a new INSUFFICIENTDATA_DIALOG or raises the existing
%      singleton*.
%
%      H = INSUFFICIENTDATA_DIALOG returns the handle to a new INSUFFICIENTDATA_DIALOG or the handle to
%      the existing singleton*.
%
%      INSUFFICIENTDATA_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSUFFICIENTDATA_DIALOG.M with the given input arguments.
%
%      INSUFFICIENTDATA_DIALOG('Property','Value',...) creates a new INSUFFICIENTDATA_DIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before insufficientData_dialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to insufficientData_dialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help insufficientData_dialog

% Last Modified by GUIDE v2.5 26-Mar-2016 19:14:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @insufficientData_dialog_OpeningFcn, ...
                   'gui_OutputFcn',  @insufficientData_dialog_OutputFcn, ...
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


% --- Executes just before insufficientData_dialog is made visible.
function insufficientData_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to insufficientData_dialog (see VARARGIN)

% Choose default command line output for insufficientData_dialog

set(hObject, 'Name', 'Additional Information needs to be defined');
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

handles.dir = varargin{1};
handles.warning_img=imread([handles.dir.images,'warningImage.png']);
axes(handles.axes1);
axis off
imshow(handles.warning_img);
mes = varargin{2};
logData = varargin{3};
set(handles.text2,'String',mes);
set(handles.text3,'String',logData);
% Update handles structure
pos = get(0,'PointerLocation');
set(handles.figure1,'Position',[pos(1)-581/2,pos(2)-212/4,581,212]);
guidata(hObject, handles);

% UIWAIT makes insufficientData_dialog wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = insufficientData_dialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
delete(hObject);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
