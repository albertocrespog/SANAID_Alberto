function varargout = dimensionado_aleron(varargin)
% DIMENSIONADO_ALERON MATLAB code for dimensionado_aleron.fig
%      DIMENSIONADO_ALERON, by itself, creates a new DIMENSIONADO_ALERON or raises the existing
%      singleton*.
%
%      H = DIMENSIONADO_ALERON returns the handle to a new DIMENSIONADO_ALERON or the handle to
%      the existing singleton*.
%
%      DIMENSIONADO_ALERON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIMENSIONADO_ALERON.M with the given input arguments.
%
%      DIMENSIONADO_ALERON('Property','Value',...) creates a new DIMENSIONADO_ALERON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dimensionado_aleron_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dimensionado_aleron_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dimensionado_aleron

% Last Modified by GUIDE v2.5 13-Aug-2015 17:18:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dimensionado_aleron_OpeningFcn, ...
                   'gui_OutputFcn',  @dimensionado_aleron_OutputFcn, ...
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


% --- Executes just before dimensionado_aleron is made visible.
function dimensionado_aleron_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dimensionado_aleron (see VARARGIN)

% Choose default command line output for dimensionado_aleron

set(hObject, 'Name', 'Ailerons Design');
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

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dimensionado_aleron wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dimensionado_aleron_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
delete(hObject);



function b_val_Callback(hObject, eventdata, handles)
% hObject    handle to b_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_val as text
%        str2double(get(hObject,'String')) returns contents of b_val as a double
handles.b = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function b_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S_val_Callback(hObject, eventdata, handles)
% hObject    handle to S_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of S_val as text
%        str2double(get(hObject,'String')) returns contents of S_val as a double
handles.S = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function S_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cr_val_Callback(hObject, eventdata, handles)
% hObject    handle to yf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yf_val as text
%        str2double(get(hObject,'String')) returns contents of yf_val as a double
handles.cr = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ct_val_Callback(hObject, eventdata, handles)
% hObject    handle to y0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y0_val as text
%        str2double(get(hObject,'String')) returns contents of y0_val as a double
handles.ct = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ct_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAM_val_Callback(hObject, eventdata, handles)
% hObject    handle to yf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yf_val as text
%        str2double(get(hObject,'String')) returns contents of yf_val as a double
handles.LAM = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function LAM_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cac_val_Callback(hObject, eventdata, handles)
% hObject    handle to cr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr_val as text
%        str2double(get(hObject,'String')) returns contents of cr_val as a double
handles.ca_c = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cac_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cla_val_Callback(hObject, eventdata, handles)
% hObject    handle to ct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ct_val as text
%        str2double(get(hObject,'String')) returns contents of ct_val as a double
handles.Cla = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cla_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tc_val_Callback(hObject, eventdata, handles)
% hObject    handle to cac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cac_val as text
%        str2double(get(hObject,'String')) returns contents of cac_val as a double
handles.t_c = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tc_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y0_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cla_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cla_val as text
%        str2double(get(hObject,'String')) returns contents of Cla_val as a double
handles.y0 = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function y0_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cla_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yf_val_Callback(hObject, eventdata, handles)
% hObject    handle to tc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tc_val as text
%        str2double(get(hObject,'String')) returns contents of tc_val as a double
handles.yf = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yf_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vcor_val_Callback(hObject, eventdata, handles)
% hObject    handle to y0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y0_val as text
%        str2double(get(hObject,'String')) returns contents of y0_val as a double
handles.Vcor = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Vcor_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cldareq_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cldareq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cldareq_val as text
%        str2double(get(hObject,'String')) returns contents of Cldareq_val as a double


% --- Executes during object creation, after setting all properties.
function Cldareq_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cldareq_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h_val_Callback(hObject, eventdata, handles)
% hObject    handle to yf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yf_val as text
%        str2double(get(hObject,'String')) returns contents of yf_val as a double
handles.h = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function h_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Clda_calc_button.
function Clda_calc_button_Callback(hObject, eventdata, handles)
% hObject    handle to Clda_calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'b') || ~isfield(handles,'S') || ~isfield(handles,'cr') || ~isfield(handles,'ct') || ~isfield(handles,'LAM') || ~isfield(handles,'ca_c') ||...
   ~isfield(handles,'Cla') || ~isfield(handles,'t_c') || ~isfield(handles,'y0') || ~isfield(handles,'yf') || ~isfield(handles,'Vcor') || ~isfield(handles,'h') ||...
   isempty(handles.b) || isempty(handles.S) || isempty(handles.cr) || isempty(handles.ct) || isempty(handles.LAM) || isempty(handles.ca_c) ||...
   isempty(handles.Cla) || isempty(handles.t_c) || isempty(handles.y0) || isempty(handles.yf) || isempty(handles.Vcor) || isempty(handles.h)
else
    [~,~,~,handles.ainf] = atmos_inter(handles.h);
    handles.Minf = handles.Vcor/handles.ainf;
    handles.AR = handles.b^2/handles.S;
    handles.LAMc4 = atan((handles.b*tan(handles.LAM)/2 - 0.25*(handles.cr-handles.ct))/handles.b/2);
    handles.TR = handles.ct/handles.cr;
    
    % Calculo de de Cl_da
    beta = sqrt(1-handles.Minf^2);
    betaCldelta_kf  = betaCldelta_k_calc(2*handles.yf/handles.b,handles.TR,handles.AR,handles.LAMc4,handles.Minf);
    betaCldelta_k0  = betaCldelta_k_calc(2*handles.y0/handles.b,handles.TR,handles.AR,handles.LAMc4,handles.Minf);
    Cldelta_prima   = 1/beta*(betaCldelta_kf - betaCldelta_k0);

    Cldelta_int    = Cldelta_Cldeltatheory_calc(handles.ca_c,handles.Cla)*Cldeltatheory_calc(handles.ca_c,handles.t_c);
    alpha_delta = Cldelta_int/handles.Cla;

    handles.Cl_da = Cldelta_prima*alpha_delta;
    set(handles.Clda_val,'String',num2str(handles.Cl_da));
    
end
guidata(hObject, handles);


function P_val_Callback(hObject, eventdata, handles)
% hObject    handle to P_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_val as text
%        str2double(get(hObject,'String')) returns contents of P_val as a double
handles.P = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function P_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function damax_val_Callback(hObject, eventdata, handles)
% hObject    handle to damax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of damax_val as text
%        str2double(get(hObject,'String')) returns contents of damax_val as a double
handles.da_max = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function damax_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to damax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Cldareq_calc_button.
function Cldareq_calc_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cldareq_calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'b') || ~isfield(handles,'S') || ~isfield(handles,'cr') || ~isfield(handles,'ct') ||...
   ~isfield(handles,'LAM') || ~isfield(handles,'Vcor') || ~isfield(handles,'h') ||...
   isempty(handles.b) || isempty(handles.S) || isempty(handles.cr) || isempty(handles.ct) ||...
   isempty(handles.LAM) || isempty(handles.Vcor) || isempty(handles.h)
else
    [~,~,~,handles.ainf] = atmos_inter(handles.h);
    handles.Minf = handles.Vcor/handles.ainf;
    handles.AR = handles.b^2/handles.S;
    handles.LAMc4 = atan((handles.b*tan(handles.LAM)/2 - 0.25*(handles.cr-handles.ct))/handles.b/2);
    handles.TR = handles.ct/handles.cr;
    handles.Clp = C_lp_calc(handles.TR,handles.AR,handles.LAMc4,handles.Minf);
    handles.Clda_req = handles.P*handles.b*handles.Clp/2/handles.Vcor/handles.da_max;
    set(handles.Cldareq_val,'String',num2str(handles.Clda_req));
end
guidata(hObject, handles);



function Clda_val_Callback(hObject, eventdata, handles)
% hObject    handle to Clda_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Clda_val as text
%        str2double(get(hObject,'String')) returns contents of Clda_val as a double


% --- Executes during object creation, after setting all properties.
function Clda_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Clda_val (see GCBO)
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
uiresume(handles.figure1);



% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
