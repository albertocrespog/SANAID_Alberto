function varargout = resultsMenu(varargin)
% RESULTSMENU MATLAB code for resultsMenu.fig
%      RESULTSMENU, by itself, creates a new RESULTSMENU or raises the existing
%      singleton*.
%
%      H = RESULTSMENU returns the handle to a new RESULTSMENU or the handle to
%      the existing singleton*.
%
%      RESULTSMENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTSMENU.M with the given input arguments.
%
%      RESULTSMENU('Property','Value',...) creates a new RESULTSMENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resultsMenu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resultsMenu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resultsMenu

% Last Modified by GUIDE v2.5 25-Aug-2015 11:21:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resultsMenu_OpeningFcn, ...
                   'gui_OutputFcn',  @resultsMenu_OutputFcn, ...
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


% --- Executes just before resultsMenu is made visible.
function resultsMenu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to resultsMenu (see VARARGIN)

% Choose default command line output for resultsMenu
handles.modelo  = varargin{1};
handles.dir     = varargin{2};
handles.output = handles.modelo;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes resultsMenu wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = resultsMenu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.modelo;
delete(handles.figure1);


% --- Executes on button press in trimLong_button.
function trimLong_button_Callback(hObject, eventdata, handles)
% hObject    handle to trimLong_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

trimado_longitudinal(handles.modelo,handles.dir);
guidata(hObject, handles);


% --- Executes on button press in OEIvientoSol_button.
function OEIvientoSol_button_Callback(hObject, eventdata, handles)
% hObject    handle to OEIvientoSol_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

OEI_vientoResults(handles.modelo,handles.dir);
guidata(hObject, handles);


% --- Executes on button press in derEst_button.
function derEst_button_Callback(hObject, eventdata, handles)
% hObject    handle to derEst_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelo = stabilityDerivatives_calc(handles.modelo, handles.dir);
guidata(hObject, handles);

% --- Executes on button press in estDin_button.
function estDin_button_Callback(hObject, eventdata, handles)
% hObject    handle to estDin_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resultDinamica(handles.modelo, handles.dir);
guidata(hObject, handles);

% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1,'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
    guidata(hObject, handles);
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
    guidata(hObject, handles);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end
