function varargout = dimensionado_vertical(varargin)
% DIMENSIONADO_VERTICAL MATLAB code for dimensionado_vertical.fig
%      DIMENSIONADO_VERTICAL, by itself, creates a new DIMENSIONADO_VERTICAL or raises the existing
%      singleton*.
%
%      H = DIMENSIONADO_VERTICAL returns the handle to a new DIMENSIONADO_VERTICAL or the handle to
%      the existing singleton*.
%
%      DIMENSIONADO_VERTICAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIMENSIONADO_VERTICAL.M with the given input arguments.
%
%      DIMENSIONADO_VERTICAL('Property','Value',...) creates a new DIMENSIONADO_VERTICAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dimensionado_vertical_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dimensionado_vertical_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dimensionado_vertical

% Last Modified by GUIDE v2.5 28-Oct-2015 15:33:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dimensionado_vertical_OpeningFcn, ...
                   'gui_OutputFcn',  @dimensionado_vertical_OutputFcn, ...
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


% --- Executes just before dimensionado_vertical is made visible.
function dimensionado_vertical_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dimensionado_vertical (see VARARGIN)

% Choose default command line output for dimensionado_vertical

set(hObject, 'Name', 'Vertical Stabilizer Design');
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

handles.type = 'convencional';
handles.F_OEI = 1.10;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dimensionado_vertical wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dimensionado_vertical_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
delete(hObject);

function Pmax_val_Callback(hObject, eventdata, handles)
% hObject    handle to Pmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pmax_val as text
%        str2double(get(hObject,'String')) returns contents of Pmax_val as a double
kW2W = 1000;
handles.Pmax = 0.5*kW2W*str2double(get(hObject,'String'));
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


function Vs_val_Callback(hObject, eventdata, handles)
% hObject    handle to Vs_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vs_val as text
%        str2double(get(hObject,'String')) returns contents of Vs_val as a double
handles.Vstall = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Vs_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vs_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function h_val_Callback(hObject, eventdata, handles)
% hObject    handle to h_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h_val as text
%        str2double(get(hObject,'String')) returns contents of h_val as a double

handles.h = str2double(get(hObject, 'String'));

[~,handles.rho,~,handles.ainf] = atmos_inter(handles.h);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function h_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drmax_val_Callback(hObject, eventdata, handles)
% hObject    handle to drmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drmax_val as text
%        str2double(get(hObject,'String')) returns contents of drmax_val as a double

handles.drmax = pi/180*str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function drmax_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CndrOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to CndrOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CndrOEI_val as text
%        str2double(get(hObject,'String')) returns contents of CndrOEI_val as a double


% --- Executes during object creation, after setting all properties.
function CndrOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CndrOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcularOEI_button.
function calcularOEI_button_Callback(hObject, eventdata, handles)
% hObject    handle to calcularOEI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles
% ~isfield(handles,'Vstall')
% ~isfield(handles,'Pmax')
% ~isfield(handles,'F_OEI')
% ~isfield(handles,'dmotor')
% ~isfield(handles,'rho')
% ~isfield(handles,'Sref')
% ~isfield(handles,'bw')
% ~isfield(handles,'rendprop')
% ~isfield(handles,'drmax')
% isempty(handles.Vstall)
% isempty(handles.Pmax)
% isempty(handles.F_OEI)
% isempty(handles.dmotor)
% isempty(handles.rho)
% isempty(handles.Sref)
% isempty(handles.bw)
% isempty(handles.rendprop)
% isempty(handles.drmax)


if  ~isfield(handles,'Vstall') || ~isfield(handles,'Pmax') || ~isfield(handles,'F_OEI') || ~isfield(handles,'dmotor') || ~isfield(handles,'rho') || ~isfield(handles,'Sref') ||...
    ~isfield(handles,'bw') || ~isfield(handles,'rendprop') || ~isfield(handles,'drmax') || isempty(handles.Vstall) || isempty(handles.Pmax) || isempty(handles.F_OEI) ||...
    isempty(handles.dmotor) || isempty(handles.rho) || isempty(handles.Sref) || isempty(handles.bw) || isempty(handles.rendprop) || isempty(handles.drmax)
    
else
    V = 1.2*handles.Vstall;
    F = handles.rendprop*handles.Pmax/V;
    N = handles.F_OEI*F*handles.dmotor;
    handles.Cndr_OEI = -2*N/handles.Sref/V^2/handles.rho/handles.bw/handles.drmax;
    set(handles.CndrOEI_val,'String',num2str(handles.Cndr_OEI));
end

guidata(hObject, handles);



function Cndr_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cndr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cndr_val as text
%        str2double(get(hObject,'String')) returns contents of Cndr_val as a double


% --- Executes during object creation, after setting all properties.
function Cndr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cndr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcularCndr_button.
function calcularCndr_button_Callback(hObject, eventdata, handles)
% hObject    handle to calcularCndr_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if      ~isfield(handles,'CLav') || ~isfield(handles,'Sref') || ~isfield(handles,'CLav') || ~isfield(handles,'lv') || ~isfield(handles,'bw')...
        || ~isfield(handles,'Sr_Sv') || ~isfield(handles,'etav') || isempty(handles.CLav) || isempty(handles.Sref) || isempty(handles.lv)...
        || isempty(handles.bw) || isempty(handles.Sr_Sv) || isempty(handles.etav)
else
    
    handles.tau =   controlEffec_calc(handles.Sr_Sv);
    handles.Cndr = -handles.CLav*handles.tau*handles.Sv/handles.Sref*handles.lv/handles.bw*handles.etav;
    if(strcmp(handles.type, 'twin_vertical'))
        % TODO: Improve twin vertical control derivative estimation
        handles.Cndr = 2*handles.Cndr;
    end
    set(handles.Cndr_val,'String',num2str(handles.Cndr));
end
guidata(hObject, handles);

% --- Executes on button press in CalcularCLav_button.
function CalcularCLav_button_Callback(hObject, eventdata, handles)
% hObject    handle to CalcularCLav_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if      ~isfield(handles,'bv') || ~isfield(handles,'crv') || ~isfield(handles,'TRv') || ~isfield(handles,'LAMv') || ~isfield(handles,'Shor')...
        || ~isfield(handles,'zhor') || ~isfield(handles,'Dfusv') || ~isfield(handles,'h') || ~isfield(handles,'Vstall') || isempty(handles.bv)...
        || isempty(handles.crv) || isempty(handles.TRv) || isempty(handles.LAMv) || isempty(handles.Shor) || isempty(handles.zhor)...
        || isempty(handles.Dfusv) || isempty(handles.h) || isempty(handles.Vstall)
else
    if ~isfield(handles,'Cla')
        handles.Cla = [];
    end
    handles.AR = handles.bv^2/handles.Sv;
    handles.ctv = handles.crv*handles.TRv;
    handles.LAMc2 = atan((handles.bv*tan(handles.LAMv) - 0.5*(handles.crv-handles.ctv))/handles.bv);
    handles.Minf = 1.2*handles.Vstall/handles.ainf;
    
    handles.AReffec = getARv_eff(handles.Sv, handles.bv, handles.TRv, handles.Shor, 0.25*handles.crv, handles.zhor, handles.Dfusv, handles.type);
    handles.CLav = getCLa_geometryBased(handles.AReffec, handles.LAMc2, handles.Minf, handles.Cla);
    
    set(handles.CLav_val,'String', num2str(handles.CLav));
end

guidata(hObject, handles);

% --- Executes on button press in pasovariable_button.
function pasovariable_button_Callback(hObject, eventdata, handles)
% hObject    handle to pasovariable_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pasovariable_button

handles.F_OEI = 1.10;
guidata(hObject, handles);


% --- Executes on button press in pasofijo_button.
function pasofijo_button_Callback(hObject, eventdata, handles)
% hObject    handle to pasofijo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pasofijo_button

handles.F_OEI = 1.25;
guidata(hObject, handles);


function Sv_val_Callback(hObject, eventdata, handles)
% hObject    handle to etav_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etav_val as text
%        str2double(get(hObject,'String')) returns contents of etav_val as a double
handles.Sv = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Sv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etav_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLav_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLav_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etav_val as text
%        str2double(get(hObject,'String')) returns contents of etav_val as a double
handles.CLav = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CLav_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLav_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Sr_Sv_val_Callback(hObject, eventdata, handles)
% hObject    handle to etav_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etav_val as text
%        str2double(get(hObject,'String')) returns contents of etav_val as a double
handles.Sr_Sv = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Sr_Sv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etav_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lv_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sr_Sv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sr_Sv_val as text
%        str2double(get(hObject,'String')) returns contents of Sr_Sv_val as a double
handles.lv = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function lv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sr_Sv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bw_val_Callback(hObject, eventdata, handles)
% hObject    handle to bw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bw_val as text
%        str2double(get(hObject,'String')) returns contents of bw_val as a double
handles.bw = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function bw_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etav_val_Callback(hObject, eventdata, handles)
% hObject    handle to lv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lv_val as text
%        str2double(get(hObject,'String')) returns contents of lv_val as a double

handles.etav = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function etav_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bv_val_Callback(hObject, eventdata, handles)
% hObject    handle to bv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bv_val as text
%        str2double(get(hObject,'String')) returns contents of bv_val as a double

handles.bv = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function bv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function crv_val_Callback(hObject, eventdata, handles)
% hObject    handle to crv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of crv_val as text
%        str2double(get(hObject,'String')) returns contents of crv_val as a double
handles.crv = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function crv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TRv_val_Callback(hObject, eventdata, handles)
% hObject    handle to TRv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TRv_val as text
%        str2double(get(hObject,'String')) returns contents of TRv_val as a double
handles.TRv = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TRv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TRv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMv_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMv_val as text
%        str2double(get(hObject,'String')) returns contents of LAMv_val as a double
handles.LAMv = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LAMv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Shor_val_Callback(hObject, eventdata, handles)
% hObject    handle to Shor_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shor_val as text
%        str2double(get(hObject,'String')) returns contents of Shor_val as a double
handles.Shor = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Shor_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shor_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zhor_val_Callback(hObject, eventdata, handles)
% hObject    handle to zhor_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zhor_val as text
%        str2double(get(hObject,'String')) returns contents of zhor_val as a double
handles.zhor = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function zhor_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zhor_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dfusv_val_Callback(hObject, eventdata, handles)
% hObject    handle to Dfusv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dfusv_val as text
%        str2double(get(hObject,'String')) returns contents of Dfusv_val as a double

handles.Dfusv = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Dfusv_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dfusv_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Minf_val_Callback(hObject, eventdata, handles)
% hObject    handle to Minf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Minf_val as text
%        str2double(get(hObject,'String')) returns contents of Minf_val as a double
handles.Minf = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Minf_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Minf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cla_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cla_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cla_val as text
%        str2double(get(hObject,'String')) returns contents of Cla_val as a double
handles.Cla = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cla_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cla_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sref_val_Callback(hObject, eventdata, handles)
% hObject    handle to bw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bw_val as text
%        str2double(get(hObject,'String')) returns contents of bw_val as a double
handles.Sref = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Sref_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in convencional_button.
function convencional_button_Callback(hObject, eventdata, handles)
% hObject    handle to convencional_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convencional_button

handles.type = 'convencional';
guidata(hObject, handles);


% --- Executes on button press in twin_button.
function twin_button_Callback(hObject, eventdata, handles)
% hObject    handle to twin_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twin_button

handles.type = 'twin_vertical';
guidata(hObject, handles);


function dmotor_val_Callback(hObject, eventdata, handles)
% hObject    handle to dmotor_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dmotor_val as text
%        str2double(get(hObject,'String')) returns contents of dmotor_val as a double
handles.dmotor = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dmotor_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dmotor_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rendprop_val_Callback(hObject, eventdata, handles)
% hObject    handle to rendprop_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rendprop_val as text
%        str2double(get(hObject,'String')) returns contents of rendprop_val as a double
handles.rendprop = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rendprop_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rendprop_val (see GCBO)
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
uiresume(handles.figure1)


% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
