function varargout = propulsiveProperties(varargin)
% PROPULSIVEPROPERTIES MATLAB code for propulsiveProperties.fig
%      PROPULSIVEPROPERTIES, by itself, creates a new PROPULSIVEPROPERTIES or raises the existing
%      singleton*.
%
%      H = PROPULSIVEPROPERTIES returns the handle to a new PROPULSIVEPROPERTIES or the handle to
%      the existing singleton*.
%
%      PROPULSIVEPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROPULSIVEPROPERTIES.M with the given input arguments.
%
%      PROPULSIVEPROPERTIES('Property','Value',...) creates a new PROPULSIVEPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before propulsiveProperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to propulsiveProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help propulsiveProperties

% Last Modified by GUIDE v2.5 17-Mar-2016 10:31:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @propulsiveProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @propulsiveProperties_OutputFcn, ...
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


% --- Executes just before propulsiveProperties is made visible.
function propulsiveProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to propulsiveProperties (see VARARGIN)

% Choose default command line output for propulsiveProperties
handles.shp2W = 745.69987;
handles.kW2W = 1000;
handles.modelo = varargin{1}; 
handles.propulsion = handles.modelo.propulsion;
handles.output = handles.propulsion;
if isempty(handles.propulsion.F_OEI)
    handles.propulsion.F_OEI = 1.25;
end
handles.dir     = varargin{2};
handles.changes = varargin{3};


handles.topView = imread('engTopView.png');
handles.sideView = imread('engSideView.png');
handles.prop = imread('propeller.png');
axes(handles.axes2)
imshow(handles.sideView);
axes(handles.axes3)
imshow(handles.topView);
axes(handles.axes4)
imshow(handles.prop);
% handles.propulsion.nBlades      = [];
% handles.propulsion.betaBlade    = [];

if handles.propulsion.F_OEI == 1.25
    set(handles.propType_val,'Value',1);
else
    set(handles.propType_val,'Value',2);
end
set(handles.Pmax_val,'String',num2str(handles.propulsion.Pmax/handles.kW2W));
set(handles.deltaT_val,'String',num2str(handles.propulsion.deltaT));
set(handles.deltaT_val,'visible','off');
set(handles.text4,'visible','off');
set(handles.nEng_val,'String',num2str(handles.propulsion.n));
set(handles.rendProp_val,'String',num2str(handles.propulsion.rendProp));
set(handles.Dprop_val,'String',num2str(handles.propulsion.Dprop));
set(handles.wR_03,'String',num2str(handles.propulsion.w_R_30));
set(handles.wR_06,'String',num2str(handles.propulsion.w_R_60));
set(handles.wR_09,'String',num2str(handles.propulsion.w_R_90));
set(handles.X_val,'String',num2str(handles.propulsion.X));
set(handles.Y_val,'String',num2str(handles.propulsion.Y));
set(handles.Z_val,'String',num2str(handles.propulsion.Z));
set(handles.A_val,'String',num2str(handles.propulsion.Acoef/handles.kW2W));
set(handles.B_val,'String',num2str(handles.propulsion.Bcoef/handles.kW2W));
set(handles.C_val,'String',num2str(handles.propulsion.Ccoef/handles.kW2W));
set(handles.nBlades_val,'String',num2str(handles.propulsion.nBlades));
set(handles.beta_val,'String',num2str(handles.propulsion.betaBlade));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes propulsiveProperties wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = propulsiveProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.changes;
delete(handles.figure1);



function Pmax_val_Callback(hObject, eventdata, handles)
% hObject    handle to Pmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pmax_val as text
%        str2double(get(hObject,'String')) returns contents of Pmax_val as a double
handles.propulsion.Pmax = str2double(get(hObject,'String'))*handles.kW2W;
if ~isempty(handles.modelo.general.h) && ~isempty(handles.modelo.general.Minf) && ~isempty(handles.propulsion.Pmax)
    cd(handles.dir.code);
    [handles.propulsion.Acoef,handles.propulsion.Bcoef,handles.propulsion.Ccoef] = getPropModel(handles.propulsion.Pmax, 1, handles.modelo.general.Minf, handles.modelo.general.h);
    cd(handles.dir.gui);
    set(handles.A_val,'String',num2str(handles.propulsion.Acoef/handles.kW2W));
    set(handles.B_val,'String',num2str(handles.propulsion.Bcoef/handles.kW2W));
    set(handles.C_val,'String',num2str(handles.propulsion.Ccoef/handles.kW2W));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Pmax_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function deltaT_val_Callback(hObject, eventdata, handles)
% hObject    handle to deltaT_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltaT_val as text
%        str2double(get(hObject,'String')) returns contents of deltaT_val as a double
handles.propulsion.deltaT = str2double(get(hObject,'String'));
if ~isempty(handles.modelo.general.h) && ~isempty(handles.modelo.general.Minf) && ~isempty(handles.propulsion.Pmax)
    cd(handles.dir.code);
    [handles.propulsion.Acoef,handles.propulsion.Bcoef,handles.propulsion.Ccoef] = getPropModel(handles.propulsion.Pmax, 1, handles.modelo.general.Minf, handles.modelo.general.h);
    cd(handles.dir.gui);
    set(handles.A_val,'String',num2str(handles.propulsion.Acoef/handles.kW2W));
    set(handles.B_val,'String',num2str(handles.propulsion.Bcoef/handles.kW2W));
    set(handles.C_val,'String',num2str(handles.propulsion.Ccoef/handles.kW2W));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function deltaT_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaT_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in propType_val.
function propType_val_Callback(hObject, eventdata, handles)
% hObject    handle to propType_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns propType_val contents as cell array
%        contents{get(hObject,'Value')} returns selected item from propType_val
if get(hObject,'Value') == 1
    handles.propulsion.F_OEI = 1.25;
else
    handles.propulsion.F_OEI = 1.10;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function propType_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propType_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nEng_val_Callback(hObject, eventdata, handles)
% hObject    handle to nEng_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nEng_val as text
%        str2double(get(hObject,'String')) returns contents of nEng_val as a double
handles.propulsion.n = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function nEng_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nEng_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dprop_val_Callback(hObject, eventdata, handles)
% hObject    handle to Dprop_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dprop_val as text
%        str2double(get(hObject,'String')) returns contents of Dprop_val as a double
handles.propulsion.Dprop = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Dprop_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dprop_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_val_Callback(hObject, eventdata, handles)
% hObject    handle to B_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_val as text
%        str2double(get(hObject,'String')) returns contents of B_val as a double
handles.propulsion.Bcoef = str2double(get(hObject,'String'))*handles.kW2W;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function B_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C_val_Callback(hObject, eventdata, handles)
% hObject    handle to C_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_val as text
%        str2double(get(hObject,'String')) returns contents of C_val as a double
handles.propulsion.Ccoef = str2double(get(hObject,'String'))*handles.kW2W;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function C_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A_val_Callback(hObject, eventdata, handles)
% hObject    handle to A_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_val as text
%        str2double(get(hObject,'String')) returns contents of A_val as a double
handles.propulsion.Acoef = str2double(get(hObject,'String'))*handles.kW2W;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function A_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_val (see GCBO)
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
if ~isequaln(handles.output,handles.propulsion)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.propulsion;
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
if ~isequaln(handles.output,handles.propulsion)
    handles.changes = 1;
end
handles.output = handles.propulsion;
 
guidata(hObject, handles);



function X_val_Callback(hObject, eventdata, handles)
% hObject    handle to X_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of X_val as text
%        str2double(get(hObject,'String')) returns contents of X_val as a double
handles.propulsion.X = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function X_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Z_val_Callback(hObject, eventdata, handles)
% hObject    handle to Z_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z_val as text
%        str2double(get(hObject,'String')) returns contents of Z_val as a double
handles.propulsion.Z = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Z_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_val_Callback(hObject, eventdata, handles)
% hObject    handle to Y_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y_val as text
%        str2double(get(hObject,'String')) returns contents of Y_val as a double
handles.propulsion.Y = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Y_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wR_03_Callback(hObject, eventdata, handles)
% hObject    handle to wR_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wR_03 as text
%        str2double(get(hObject,'String')) returns contents of wR_03 as a double
handles.propulsion.w_R_30 = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function wR_03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wR_03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wR_06_Callback(hObject, eventdata, handles)
% hObject    handle to wR_06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wR_06 as text
%        str2double(get(hObject,'String')) returns contents of wR_06 as a double
handles.propulsion.w_R_60 = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function wR_06_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wR_06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wR_09_Callback(hObject, eventdata, handles)
% hObject    handle to wR_09 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wR_09 as text
%        str2double(get(hObject,'String')) returns contents of wR_09 as a double
handles.propulsion.w_R_90 = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function wR_09_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wR_09 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rendProp_val_Callback(hObject, eventdata, handles)
% hObject    handle to rendProp_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rendProp_val as text
%        str2double(get(hObject,'String')) returns contents of rendProp_val as a double
handles.propulsion.rendProp = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function rendProp_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rendProp_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nBlades_val_Callback(hObject, eventdata, handles)
% hObject    handle to nBlades_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nBlades_val as text
%        str2double(get(hObject,'String')) returns contents of nBlades_val as a double
handles.propulsion.nBlades = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nBlades_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nBlades_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_val_Callback(hObject, eventdata, handles)
% hObject    handle to beta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_val as text
%        str2double(get(hObject,'String')) returns contents of beta_val as a double
handles.propulsion.betaBlade = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function beta_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if ~isequaln(handles.output,handles.propulsion)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.propulsion;
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


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen([handles.dir.help,'propulsion_guide.pdf']);
