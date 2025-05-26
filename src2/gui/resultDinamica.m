function varargout = resultDinamica(varargin)
% RESULTDINAMICA MATLAB code for resultDinamica.fig
%      RESULTDINAMICA, by itself, creates a new RESULTDINAMICA or raises the existing
%      singleton*.
%
%      H = RESULTDINAMICA returns the handle to a new RESULTDINAMICA or the handle to
%      the existing singleton*.
%
%      RESULTDINAMICA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTDINAMICA.M with the given input arguments.
%
%      RESULTDINAMICA('Property','Value',...) creates a new RESULTDINAMICA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resultDinamica_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resultDinamica_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resultDinamica

% Last Modified by GUIDE v2.5 26-Oct-2015 16:48:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resultDinamica_OpeningFcn, ...
                   'gui_OutputFcn',  @resultDinamica_OutputFcn, ...
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


% --- Executes just before resultDinamica is made visible.
function resultDinamica_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to resultDinamica (see VARARGIN)

% Choose default command line output for resultDinamica
handles.modelo  = varargin{1};
handles.dir         = varargin{2};


handles.results     = dynamic_analysis(handles.modelo);



set(handles.eigLon_1,'String',num2str(handles.results.long.pole1));
set(handles.eigLon_2,'String',num2str(handles.results.long.pole2));
set(handles.eigLon_3,'String',num2str(handles.results.long.pole3));
set(handles.eigLon_4,'String',num2str(handles.results.long.pole4));

set(handles.wnSP_val,'String',num2str(handles.results.long.SP_wn));
set(handles.dampSP_val,'String',num2str(handles.results.long.SP_damp));
set(handles.tSP_val,'String',num2str(handles.results.long.SP_T));
set(handles.t2SP_val,'String',num2str(handles.results.long.SP_T2));

set(handles.wnPH_val,'String',num2str(handles.results.long.PH_wn));
set(handles.dampPH_val,'String',num2str(handles.results.long.PH_damp));
set(handles.tPH_val,'String',num2str(handles.results.long.PH_T));
set(handles.t2PH_val,'String',num2str(handles.results.long.PH_T2));


set(handles.eigLat_1,'String',num2str(handles.results.lat.pole1));
set(handles.eigLat_2,'String',num2str(handles.results.lat.pole2));
set(handles.eigLat_3,'String',num2str(handles.results.lat.pole3));
set(handles.eigLat_4,'String',num2str(handles.results.lat.pole4));
set(handles.eigLat_5,'String',num2str(handles.results.lat.pole5));

set(handles.wnDR_val,'String',num2str(handles.results.lat.DR_wn));
set(handles.dampDR_val,'String',num2str(handles.results.lat.DR_damp));
set(handles.tDR_val,'String',num2str(handles.results.lat.DR_T));
set(handles.t2DR_val,'String',num2str(handles.results.lat.DR_T2));

if handles.results.lat.Spiral_pole > 0
    set(handles.t2ESP_txt,'String', 'Double Time (s)');
else
    set(handles.t2ESP_txt,'String', 'Half Time (s)');
end
set(handles.t2ESP_val,'String',num2str(handles.results.lat.ESP_T2));

if handles.results.lat.Roll_pole > 0
    set(handles.t2ROL_txt,'String', 'Double Time (s)');
else
    set(handles.t2ROL_txt,'String', 'Half Time (s)');
end
set(handles.t2ROL_val,'String',num2str(handles.results.lat.ROL_T2));



handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes resultDinamica wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = resultDinamica_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function eigLon_1_Callback(hObject, eventdata, handles)
% hObject    handle to eigLon_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLon_1 as text
%        str2double(get(hObject,'String')) returns contents of eigLon_1 as a double


% --- Executes during object creation, after setting all properties.
function eigLon_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLon_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLon_2_Callback(hObject, eventdata, handles)
% hObject    handle to eigLon_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLon_2 as text
%        str2double(get(hObject,'String')) returns contents of eigLon_2 as a double


% --- Executes during object creation, after setting all properties.
function eigLon_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLon_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
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



function eigLon_3_Callback(hObject, eventdata, handles)
% hObject    handle to eigLon_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLon_3 as text
%        str2double(get(hObject,'String')) returns contents of eigLon_3 as a double


% --- Executes during object creation, after setting all properties.
function eigLon_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLon_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLon_4_Callback(hObject, eventdata, handles)
% hObject    handle to eigLon_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLon_4 as text
%        str2double(get(hObject,'String')) returns contents of eigLon_4 as a double


% --- Executes during object creation, after setting all properties.
function eigLon_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLon_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wnPH_val_Callback(hObject, eventdata, handles)
% hObject    handle to wnPH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wnPH_val as text
%        str2double(get(hObject,'String')) returns contents of wnPH_val as a double


% --- Executes during object creation, after setting all properties.
function wnPH_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wnPH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dampPH_val_Callback(hObject, eventdata, handles)
% hObject    handle to dampPH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dampPH_val as text
%        str2double(get(hObject,'String')) returns contents of dampPH_val as a double


% --- Executes during object creation, after setting all properties.
function dampPH_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dampPH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2PH_val_Callback(hObject, eventdata, handles)
% hObject    handle to t2PH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2PH_val as text
%        str2double(get(hObject,'String')) returns contents of t2PH_val as a double


% --- Executes during object creation, after setting all properties.
function t2PH_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2PH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tPH_val_Callback(hObject, eventdata, handles)
% hObject    handle to tPH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tPH_val as text
%        str2double(get(hObject,'String')) returns contents of tPH_val as a double


% --- Executes during object creation, after setting all properties.
function tPH_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tPH_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wnSP_val_Callback(hObject, eventdata, handles)
% hObject    handle to wnSP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wnSP_val as text
%        str2double(get(hObject,'String')) returns contents of wnSP_val as a double


% --- Executes during object creation, after setting all properties.
function wnSP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wnSP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dampSP_val_Callback(hObject, eventdata, handles)
% hObject    handle to dampSP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dampSP_val as text
%        str2double(get(hObject,'String')) returns contents of dampSP_val as a double


% --- Executes during object creation, after setting all properties.
function dampSP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dampSP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2SP_val_Callback(hObject, eventdata, handles)
% hObject    handle to t2SP_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2SP_txt as text
%        str2double(get(hObject,'String')) returns contents of t2SP_txt as a double


% --- Executes during object creation, after setting all properties.
function t2SP_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2SP_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tSP_val_Callback(hObject, eventdata, handles)
% hObject    handle to tSP_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tSP_txt as text
%        str2double(get(hObject,'String')) returns contents of tSP_txt as a double


% --- Executes during object creation, after setting all properties.
function tSP_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tSP_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wnDR_val_Callback(hObject, eventdata, handles)
% hObject    handle to wnDR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wnDR_val as text
%        str2double(get(hObject,'String')) returns contents of wnDR_val as a double


% --- Executes during object creation, after setting all properties.
function wnDR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wnDR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dampDR_val_Callback(hObject, eventdata, handles)
% hObject    handle to dampDR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dampDR_val as text
%        str2double(get(hObject,'String')) returns contents of dampDR_val as a double


% --- Executes during object creation, after setting all properties.
function dampDR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dampDR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2DR_val_Callback(hObject, eventdata, handles)
% hObject    handle to t2DR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2DR_val as text
%        str2double(get(hObject,'String')) returns contents of t2DR_val as a double


% --- Executes during object creation, after setting all properties.
function t2DR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2DR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tDR_val_Callback(hObject, eventdata, handles)
% hObject    handle to tDR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tDR_val as text
%        str2double(get(hObject,'String')) returns contents of tDR_val as a double


% --- Executes during object creation, after setting all properties.
function tDR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tDR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wnESP_val_Callback(hObject, eventdata, handles)
% hObject    handle to wnESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wnESP_val as text
%        str2double(get(hObject,'String')) returns contents of wnESP_val as a double


% --- Executes during object creation, after setting all properties.
function wnESP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wnESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dampESP_val_Callback(hObject, eventdata, handles)
% hObject    handle to dampESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dampESP_val as text
%        str2double(get(hObject,'String')) returns contents of dampESP_val as a double


% --- Executes during object creation, after setting all properties.
function dampESP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dampESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2ESP_val_Callback(hObject, eventdata, handles)
% hObject    handle to t2ESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2ESP_val as text
%        str2double(get(hObject,'String')) returns contents of t2ESP_val as a double


% --- Executes during object creation, after setting all properties.
function t2ESP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2ESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tESP_val_Callback(hObject, eventdata, handles)
% hObject    handle to tESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tESP_val as text
%        str2double(get(hObject,'String')) returns contents of tESP_val as a double


% --- Executes during object creation, after setting all properties.
function tESP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tESP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLat_1_Callback(hObject, eventdata, handles)
% hObject    handle to eigLat_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLat_1 as text
%        str2double(get(hObject,'String')) returns contents of eigLat_1 as a double


% --- Executes during object creation, after setting all properties.
function eigLat_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLat_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLat_2_Callback(hObject, eventdata, handles)
% hObject    handle to eigLat_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLat_2 as text
%        str2double(get(hObject,'String')) returns contents of eigLat_2 as a double


% --- Executes during object creation, after setting all properties.
function eigLat_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLat_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLat_3_Callback(hObject, eventdata, handles)
% hObject    handle to eigLat_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLat_3 as text
%        str2double(get(hObject,'String')) returns contents of eigLat_3 as a double


% --- Executes during object creation, after setting all properties.
function eigLat_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLat_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLat_4_Callback(hObject, eventdata, handles)
% hObject    handle to eigLat_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLat_4 as text
%        str2double(get(hObject,'String')) returns contents of eigLat_4 as a double


% --- Executes during object creation, after setting all properties.
function eigLat_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLat_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigLat_5_Callback(hObject, eventdata, handles)
% hObject    handle to eigLat_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eigLat_5 as text
%        str2double(get(hObject,'String')) returns contents of eigLat_5 as a double


% --- Executes during object creation, after setting all properties.
function eigLat_5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigLat_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function t2SP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2SP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function tSP_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tSP_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2ROL_val_Callback(hObject, eventdata, handles)
% hObject    handle to t2ROL_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2ROL_val as text
%        str2double(get(hObject,'String')) returns contents of t2ROL_val as a double


% --- Executes during object creation, after setting all properties.
function t2ROL_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2ROL_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
