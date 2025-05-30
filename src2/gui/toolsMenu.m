function varargout = toolsMenu(varargin)
% TOOLSMENU MATLAB code for toolsMenu.fig
%      TOOLSMENU, by itself, creates a new TOOLSMENU or raises the existing
%      singleton*.
%
%      H = TOOLSMENU returns the handle to a new TOOLSMENU or the handle to
%      the existing singleton*.
%
%      TOOLSMENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOOLSMENU.M with the given input arguments.
%
%      TOOLSMENU('Property','Value',...) creates a new TOOLSMENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toolsMenu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toolsMenu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toolsMenu

% Last Modified by GUIDE v2.5 18-Aug-2015 17:57:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toolsMenu_OpeningFcn, ...
                   'gui_OutputFcn',  @toolsMenu_OutputFcn, ...
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


% --- Executes just before toolsMenu is made visible.
function toolsMenu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toolsMenu (see VARARGIN)

% Choose default command line output for toolsMenu
%handles.output = hObject;
handles.dir = varargin{1};
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes toolsMenu wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toolsMenu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
delete(handles.figure1);


% --- Executes on button press in longAdjust_button.
function longAdjust_button_Callback(hObject, eventdata, handles)
% hObject    handle to longAdjust_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pre_longAdjustTool(handles.dir);
guidata(hObject, handles);



% --- Executes on button press in vertDim_button.
function vertDim_button_Callback(hObject, eventdata, handles)
% hObject    handle to vertDim_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dimensionado_vertical(handles.dir);
guidata(hObject, handles);


% --- Executes on button press in aleironDim_button.
function aleironDim_button_Callback(hObject, eventdata, handles)
% hObject    handle to aleironDim_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dimensionado_aleron(handles.dir);
guidata(hObject, handles);


% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1,'waitstatus'), 'waiting')
     % The GUI is still in UIWAIT, us UIRESUME
     uiresume(handles.figure1);
     else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1,'waitstatus'), 'waiting')
     % The GUI is still in UIWAIT, us UIRESUME
     uiresume(handles.figure1);
     else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end
