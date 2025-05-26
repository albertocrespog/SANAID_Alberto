function varargout = weightProperties(varargin)
% WEIGHTPROPERTIES MATLAB code for weightProperties.fig
%      WEIGHTPROPERTIES, by itself, creates a new WEIGHTPROPERTIES or raises the existing
%      singleton*.
%
%      H = WEIGHTPROPERTIES returns the handle to a new WEIGHTPROPERTIES or the handle to
%      the existing singleton*.
%
%      WEIGHTPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WEIGHTPROPERTIES.M with the given input arguments.
%
%      WEIGHTPROPERTIES('Property','Value',...) creates a new WEIGHTPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before weightProperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to weightProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help weightProperties

% Last Modified by GUIDE v2.5 16-Oct-2015 16:55:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @weightProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @weightProperties_OutputFcn, ...
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


% --- Executes just before weightProperties is made visible.
function weightProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to weightProperties (see VARARGIN)

% Choose default command line output for weightProperties
handles.output = varargin{1};
handles.general = handles.output;
%     handles.general.Ixx = [];
%     handles.general.Iyy = [];
%     handles.general.Izz = [];
%     handles.general.Iyz = [];



set(handles.mtow_val,'String',num2str(handles.general.mtow));
set(handles.Xcg_val,'String',num2str(handles.general.Xcg));
set(handles.Ixx_val,'String',num2str(handles.general.Ixx));
set(handles.Iyy_val,'String',num2str(handles.general.Iyy));
set(handles.Izz_val,'String',num2str(handles.general.Izz));
set(handles.Ixz_val,'String',num2str(handles.general.Ixz));
handles.changes = varargin{3};
handles.dir     = varargin{2};
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes weightProperties wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = weightProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.changes;
delete(handles.figure1);


function mtow_val_Callback(hObject, eventdata, handles)
% hObject    handle to mtow_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mtow_val as text
%        str2double(get(hObject,'String')) returns contents of mtow_val as a double
handles.general.mtow = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function mtow_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mtow_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xcg_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xcg_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xcg_val as text
%        str2double(get(hObject,'String')) returns contents of Xcg_val as a double
handles.general.Xcg = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xcg_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xcg_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ixx_val_Callback(hObject, eventdata, handles)
% hObject    handle to Ixx_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ixx_val as text
%        str2double(get(hObject,'String')) returns contents of Ixx_val as a double
handles.general.Ixx = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function Ixx_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ixx_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Iyy_val_Callback(hObject, eventdata, handles)
% hObject    handle to Iyy_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Iyy_val as text
%        str2double(get(hObject,'String')) returns contents of Iyy_val as a double
handles.general.Iyy = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Iyy_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Iyy_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Izz_val_Callback(hObject, eventdata, handles)
% hObject    handle to Izz_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Izz_val as text
%        str2double(get(hObject,'String')) returns contents of Izz_val as a double
handles.general.Izz = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Izz_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Izz_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ixz_val_Callback(hObject, eventdata, handles)
% hObject    handle to Ixz_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ixz_val as text
%        str2double(get(hObject,'String')) returns contents of Ixz_val as a double
handles.general.Ixz = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Ixz_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ixz_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isequaln(handles.output,handles.general)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.general;
            handles.changes = 1;
            guidata(hObject, handles);
            uiresume(handles.figure1);
        case 'DISCARD'
            guidata(hObject, handles);
            uiresume(handles.figure1);
        case 'CANCEL'
            % Do nothing
    end
else
    uiresume(handles.figure1);
end


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isequaln(handles.output, handles.general)
    handles.changes = 1;
end

handles.output = handles.general;
guidata(hObject, handles);
