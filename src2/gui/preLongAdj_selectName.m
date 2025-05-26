function varargout = preLongAdj_selectName(varargin)
% PRELONGADJ_SELECTNAME MATLAB code for preLongAdj_selectName.fig
%      PRELONGADJ_SELECTNAME, by itself, creates a new PRELONGADJ_SELECTNAME or raises the existing
%      singleton*.
%
%      H = PRELONGADJ_SELECTNAME returns the handle to a new PRELONGADJ_SELECTNAME or the handle to
%      the existing singleton*.
%
%      PRELONGADJ_SELECTNAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRELONGADJ_SELECTNAME.M with the given input arguments.
%
%      PRELONGADJ_SELECTNAME('Property','Value',...) creates a new PRELONGADJ_SELECTNAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preLongAdj_selectName_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preLongAdj_selectName_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preLongAdj_selectName

% Last Modified by GUIDE v2.5 18-Jan-2016 17:03:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preLongAdj_selectName_OpeningFcn, ...
                   'gui_OutputFcn',  @preLongAdj_selectName_OutputFcn, ...
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


% --- Executes just before preLongAdj_selectName is made visible.
function preLongAdj_selectName_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preLongAdj_selectName (see VARARGIN)

% Choose default command line output for preLongAdj_selectName

handles.output = varargin{1};

set(hObject, 'Name', 'Select a name for the new model');
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

set(handles.nameText,'String',handles.output);
handles.save = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preLongAdj_selectName wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = preLongAdj_selectName_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.save;
delete(handles.figure1);



function nameText_Callback(hObject, eventdata, handles)
% hObject    handle to nameText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nameText as text
%        str2double(get(hObject,'String')) returns contents of nameText as a double
handles.output = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nameText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nameText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.save = 1;

if isequal(get(handles.figure1,'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    guidata(hObject, handles);
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
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
