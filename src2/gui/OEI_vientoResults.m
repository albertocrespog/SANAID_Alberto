function varargout = OEI_vientoResults(varargin)
% OEI_VIENTORESULTS MATLAB code for OEI_vientoResults.fig
%      OEI_VIENTORESULTS, by itself, creates a new OEI_VIENTORESULTS or raises the existing
%      singleton*.
%
%      H = OEI_VIENTORESULTS returns the handle to a new OEI_VIENTORESULTS or the handle to
%      the existing singleton*.
%
%      OEI_VIENTORESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OEI_VIENTORESULTS.M with the given input arguments.
%
%      OEI_VIENTORESULTS('Property','Value',...) creates a new OEI_VIENTORESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OEI_vientoResults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OEI_vientoResults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OEI_vientoResults

% Last Modified by GUIDE v2.5 25-Jan-2016 18:45:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OEI_vientoResults_OpeningFcn, ...
                   'gui_OutputFcn',  @OEI_vientoResults_OutputFcn, ...
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


% --- Executes just before OEI_vientoResults is made visible.
function OEI_vientoResults_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OEI_vientoResults (see VARARGIN)

% Choose default command line output for OEI_vientoResults

set(hObject, 'Name', 'OEI & Sideslip conditions analysis');
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


handles.output  = hObject;

handles.modelo = varargin{1};
handles.dir     = varargin{2};

axes(handles.axes1);
handles.systemImage = imread('imagenSistema2.png');
imshow(handles.systemImage);

% handles.der.Cybeta  = -1.1874;
% handles.der.Cyda    = 0;
% handles.der.Cydr    = 0.5384;
% handles.der.Clbeta  = -0.1566;
% handles.der.Clda    = 0.163;
% handles.der.Cldr    = 0.0199;
% handles.der.Cnbeta  = 0.0515;
% handles.der.Cnda    = -0.012;
% handles.der.Cndr    = -0.1445;

switch(handles.modelo.confVert)
    case 'no_vert'
        vertField = fieldnames(handles.modelo.vertical);
        nField = length(vertField);
        for k = 1:nField
            handles.modelo.vertical.(vertField{k}) = 0;
        end
    case 'convencional' 
        handles.vertFact = 1;
        
    case 'twin_vertical'
        handles.vertFact = 2;
        handles.modelo.vertical.S = handles.modelo.vertical.S*handles.vertFact;
end

handles = getDerivatives(handles);
% handles.origColor   = [0.83 0.82 0.78];
% handles.modColor    = [1, 0.6, 0.6];

% cd(handles.dir.code);
% handles.der.Cybeta  = getCybeta(handles.modelo);
% handles.der.Cyda    = getCyda(handles.modelo);
% handles.der.Cydr    = getCydr(handles.modelo);
% handles.der.Clbeta  = getClbeta(handles.modelo);
% handles.der.Clda    = getClda(handles.modelo);
% handles.der.Cldr    = getCldr(handles.modelo,0);
% handles.der.Cnbeta  = getCnbeta(handles.modelo);
% handles.der.Cnda    = getCnda(handles.modelo);
% handles.der.Cndr    = getCndr(handles.modelo, 0);
% cd(handles.dir.gui);
% 
% handles.der.matrix  =  [handles.der.Cybeta, handles.der.Cyda, handles.der.Cydr; 
%                         handles.der.Clbeta, handles.der.Clda, handles.der.Cldr; 
%                         handles.der.Cnbeta, handles.der.Cnda, handles.der.Cndr];
%                     
% 
% set(handles.Cybeta_val, 'String', num2str(handles.der.Cybeta));
% set(handles.Cyda_val, 'String', num2str(handles.der.Cyda));
% set(handles.Cydr_val, 'String', num2str(handles.der.Cydr));
% set(handles.Clbeta_val, 'String', num2str(handles.der.Clbeta));
% set(handles.Clda_val, 'String', num2str(handles.der.Clda));
% set(handles.Cldr_val, 'String', num2str(handles.der.Cldr));
% set(handles.Cnbeta_val, 'String', num2str(handles.der.Cnbeta));
% set(handles.Cnda_val, 'String', num2str(handles.der.Cnda));
% set(handles.Cndr_val, 'String', num2str(handles.der.Cndr));

set(handles.betaMaxOEI_val,'Visible', 'off');
set(handles.VmaxOEI_val,'Visible', 'off');
set(handles.PmaxOEI_val,'Visible', 'off');
set(handles.ejeX_popup,'Visible', 'off');
set(handles.OEIbarrido_button,'Visible', 'off');
set(handles.MinOEI_text,'Visible', 'off');
set(handles.MaxOEI_text,'Visible', 'off');
set(handles.ejeX_text,'Visible', 'off');


set(handles.VMaxViento_val,'Visible', 'off');
set(handles.betaMaxViento_val,'Visible', 'off');
set(handles.MinViento_text,'Visible', 'off');
set(handles.MaxViento_text,'Visible', 'off');
set(handles.Vientobarrido_button,'Visible', 'off');



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OEI_vientoResults wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OEI_vientoResults_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
delete(hObject);

function handles = getDerivatives(handles)
handles.origColor   = [0.83 0.82 0.78];
handles.modColor    = [1, 0.6, 0.6];

%cd(handles.dir.code);
handles.der.Cybeta  = getCybeta(handles.modelo);
handles.der.Cyda    = getCyda(handles.modelo);
handles.der.Cydr    = getCydr(handles.modelo);
handles.der.Clbeta  = getClbeta(handles.modelo);
handles.der.Clda    = getClda(handles.modelo);
handles.der.Cldr    = getCldr(handles.modelo,0);
handles.der.Cnbeta  = getCnbeta(handles.modelo);
handles.der.Cnda    = getCnda(handles.modelo);
handles.der.Cndr    = getCndr(handles.modelo, 0);
%cd(handles.dir.gui);

handles.der.matrix  =  refreshMatrix(handles.der);
%[handles.der.Cybeta, handles.der.Cyda, handles.der.Cydr; 
%                        handles.der.Clbeta, handles.der.Clda, handles.der.Cldr; 
%                        handles.der.Cnbeta, handles.der.Cnda, handles.der.Cndr];
                    
set(handles.Cybeta_val, 'String', num2str(handles.der.Cybeta), 'BackgroundColor', handles.origColor);
set(handles.Cyda_val, 'String', num2str(handles.der.Cyda), 'BackgroundColor', handles.origColor);
set(handles.Cydr_val, 'String', num2str(handles.der.Cydr), 'BackgroundColor', handles.origColor);
set(handles.Clbeta_val, 'String', num2str(handles.der.Clbeta), 'BackgroundColor', handles.origColor);
set(handles.Clda_val, 'String', num2str(handles.der.Clda), 'BackgroundColor', handles.origColor);
set(handles.Cldr_val, 'String', num2str(handles.der.Cldr), 'BackgroundColor', handles.origColor);
set(handles.Cnbeta_val, 'String', num2str(handles.der.Cnbeta), 'BackgroundColor', handles.origColor);
set(handles.Cnda_val, 'String', num2str(handles.der.Cnda), 'BackgroundColor', handles.origColor);
set(handles.Cndr_val, 'String', num2str(handles.der.Cndr), 'BackgroundColor', handles.origColor);
                    
function matrix = refreshMatrix(der)
matrix  =  [der.Cybeta, der.Cyda, der.Cydr;
            der.Clbeta, der.Clda, der.Cldr; 
            der.Cnbeta, der.Cnda, der.Cndr];
%set(handles.Cybeta_val, 'String', num2str(handles.der.Cybeta));
%set(handles.Cyda_val, 'String', num2str(handles.der.Cyda));
%set(handles.Cydr_val, 'String', num2str(handles.der.Cydr));
%set(handles.Clbeta_val, 'String', num2str(handles.der.Clbeta));
%set(handles.Clda_val, 'String', num2str(handles.der.Clda));
%set(handles.Cldr_val, 'String', num2str(handles.der.Cldr));
%set(handles.Cnbeta_val, 'String', num2str(handles.der.Cnbeta));
%set(handles.Cnda_val, 'String', num2str(handles.der.Cnda));
%set(handles.Cndr_val, 'String', num2str(handles.der.Cndr));



function Cybeta_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cybeta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cybeta_val as text
%        str2double(get(hObject,'String')) returns contents of Cybeta_val as a double

set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cybeta = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cybeta_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cybeta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cyda_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cyda_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cyda_val as text
%        str2double(get(hObject,'String')) returns contents of Cyda_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cyda = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cyda_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cyda_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cydr_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cydr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cydr_val as text
%        str2double(get(hObject,'String')) returns contents of Cydr_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cydr = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cydr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cydr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Clbeta_val_Callback(hObject, eventdata, handles)
% hObject    handle to Clbeta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Clbeta_val as text
%        str2double(get(hObject,'String')) returns contents of Clbeta_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Clbeta = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Clbeta_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Clbeta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Clda_val_Callback(hObject, eventdata, handles)
% hObject    handle to Clda_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Clda_val as text
%        str2double(get(hObject,'String')) returns contents of Clda_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Clda = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


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



function Cldr_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cldr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cldr_val as text
%        str2double(get(hObject,'String')) returns contents of Cldr_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cldr = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cldr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cldr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cnbeta_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cnbeta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cnbeta_val as text
%        str2double(get(hObject,'String')) returns contents of Cnbeta_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cnbeta = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cnbeta_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cnbeta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cnda_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cnda_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cnda_val as text
%        str2double(get(hObject,'String')) returns contents of Cnda_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cnda = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Cnda_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cnda_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cndr_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cndr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cndr_val as text
%        str2double(get(hObject,'String')) returns contents of Cndr_val as a double
set(hObject, 'BackgroundColor', handles.modColor);
handles.der.Cndr = str2double(get(hObject,'String'));
handles.der.matrix  =  refreshMatrix(handles.der);
guidata(hObject, handles);


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


% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = getDerivatives(handles);
% cd(handles.dir.code);
% handles.der.Cybeta  = getCybeta(handles.modelo);
% handles.der.Cyda    = getCyda(handles.modelo);
% handles.der.Cydr    = getCydr(handles.modelo);
% handles.der.Clbeta  = getClbeta(handles.modelo);
% handles.der.Clda    = getClda(handles.modelo);
% handles.der.Cldr    = getCldr(handles.modelo,0);
% handles.der.Cnbeta  = getCnbeta(handles.modelo);
% handles.der.Cnda    = getCnda(handles.modelo);
% handles.der.Cndr    = getCndr(handles.modelo, 0);
% cd(handles.dir.gui);
% 
% 
% handles.der.matrix  =  [handles.der.Cybeta, handles.der.Cyda, handles.der.Cydr; 
%                         handles.der.Clbeta, handles.der.Clda, handles.der.Cldr; 
%                         handles.der.Cnbeta, handles.der.Cnda, handles.der.Cndr];
%  
% set(handles.Cybeta_val, 'String', num2str(handles.der.Cybeta), 'BackgroundColor', handles.origColor);
% set(handles.Cyda_val, 'String', num2str(handles.der.Cyda), 'BackgroundColor', handles.origColor);
% set(handles.Cydr_val, 'String', num2str(handles.der.Cydr), 'BackgroundColor', handles.origColor);
% set(handles.Clbeta_val, 'String', num2str(handles.der.Clbeta), 'BackgroundColor', handles.origColor);
% set(handles.Clda_val, 'String', num2str(handles.der.Clda), 'BackgroundColor', handles.origColor);
% set(handles.Cldr_val, 'String', num2str(handles.der.Cldr), 'BackgroundColor', handles.origColor);
% set(handles.Cnbeta_val, 'String', num2str(handles.der.Cnbeta), 'BackgroundColor', handles.origColor);
% set(handles.Cnda_val, 'String', num2str(handles.der.Cnda), 'BackgroundColor', handles.origColor);
% set(handles.Cndr_val, 'String', num2str(handles.der.Cndr), 'BackgroundColor', handles.origColor);

guidata(hObject, handles);

function beta_val_Callback(hObject, eventdata, handles)
% hObject    handle to beta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_val as text
%        str2double(get(hObject,'String')) returns contents of beta_val as a double


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



% --- Executes on selection change in ejeX_popup.
function ejeX_popup_Callback(hObject, eventdata, handles)
% hObject    handle to ejeX_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ejeX_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ejeX_popup

switch get(hObject, 'Value')
    case 1
        handles.ejeX = 'V';
        set(handles.VmaxOEI_val, 'Visible', 'on');
        set(handles.PmaxOEI_val, 'Visible', 'off');
    case 2
        handles.ejeX = 'P';
        set(handles.VmaxOEI_val, 'Visible', 'off');
        set(handles.PmaxOEI_val, 'Visible', 'on');
end


% --- Executes during object creation, after setting all properties.
function ejeX_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ejeX_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in barridoOEI_checkbox.
function barridoOEI_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to barridoOEI_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of barridoOEI_checkbox

set(hObject,'Value',0);
switch get(handles.barridoOEI_checkbox,'Value')    
    case 1
        set(handles.ejeX_popup,'Value', 1);
        set(handles.betaMaxOEI_val,'Visible', 'on');
        set(handles.VmaxOEI_val,'Visible', 'on');
        set(handles.PmaxOEI_val, 'Visible', 'off');
        set(handles.ejeX_popup,'Visible', 'on');
        set(handles.OEIbarrido_button,'Visible', 'on');
        set(handles.MinOEI_text,'Visible', 'on');
        set(handles.MaxOEI_text,'Visible', 'on');
        set(handles.ejeX_text,'Visible', 'on');
    case 0
        set(handles.betaMaxOEI_val,'Visible', 'off');
        set(handles.VmaxOEI_val,'Visible', 'off');
        set(handles.PmaxOEI_val, 'Visible', 'off');
        set(handles.ejeX_popup,'Visible', 'off');
        set(handles.OEIbarrido_button,'Visible', 'off');
        set(handles.MinOEI_text,'Visible', 'off');
        set(handles.MaxOEI_text,'Visible', 'off');
        set(handles.ejeX_text,'Visible', 'off');
end

guidata(hObject, handles);


% --- Executes on button press in OEIbarrido_button.
function OEIbarrido_button_Callback(hObject, eventdata, handles)
% hObject    handle to OEIbarrido_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nBeta   = 5;
handles.graphColor      = hsv(handles.nBeta);
handles.beta    = linspace(handles.betaMinOEI, handles.betaMaxOEI, handles.nBeta);

switch handles.ejeX
    case 'V'
        handles.V       = handles.VminOEI:0.01:handles.VmaxOEI;
        handles.P       = handles.PminOEI;
    case 'P'
        handles.P       = handles.PminOEI:0.01:handles.PmaxOEI;
        handles.V       = handles.VminOEI;
end

[handles.phi, handles.da, handles.dr] = OEIsystemResol('OEI', handles.modelo, handles.der.matrix, handles.rhoOEI, handles.beta, handles.V, handles.P);

[handles.fileName, handles.pathName] = uiputfile('*.png', 'Guardar resultados como', '../');
phiFigure   = figure();
hold on;
grid on;
switch handles.ejeX
    case 'V'
        for jBeta = 1:handles.nBeta
            plot(V, handles.phi(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
            phiLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi), 'º'];
        end
        xlabel('V');
        ylabel('\phi');
        legend(phiLegend);
    case 'P'
        for jBeta = 1:handles.nBeta
            plot(P, handles.phi(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
            phiLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi), 'º'];
        end
        xlabel('P');
        ylabel('\phi');
        legend(phiLegend);
end

set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print('-dpng','-r300',[handles.pathName,'phi_', handles.fileName]);
close(phiFigure);

daFigure    = figure();
hold on;
grid on;
switch handles.ejeX
    case 'V'
        for jBeta = 1:handles.nBeta
            plot(V, handles.da(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
            daLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi), 'º'];
        end
        xlabel('V');
        ylabel('\delta_a');
        legend(daLegend);
    case 'P'
        for jBeta = 1:handles.nBeta
            plot(P, handles.da(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
            daLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi), 'º'];
        end
        xlabel('P');
        ylabel('\delta_a');
        legend(daLegend);
end
legend(daLegend)
set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print('-dpng','-r300',[handles.pathName,'da_', handles.fileName]);
close(daFigure);

drFigure    = figure();
hold on;
grid on;
switch handles.ejeX
    case 'V'
        for jBeta = 1:handles.nBeta
            plot(V, handles.dr(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
            drLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi), 'º'];
        end
        xlabel('V');
        ylabel('\delta_r');
        legend(drLegend);
    case 'P'
        for jBeta = 1:handles.nBeta
            plot(P, handles.dr(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
            drLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi), 'º'];
        end
        xlabel('P');
        ylabel('\delta_r');
        legend(drLegend);
end

set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print('-dpng','-r300',[handles.pathName,'dr_', handles.fileName]);
close(drFigure);

winopen(handles.pathName);

guidata(hObject, handles);


% --- Executes on button press in OEISystem_button.
function OEISystem_button_Callback(hObject, eventdata, handles)
% hObject    handle to OEISystem_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.phi, handles.da, handles.dr] = OEIsystemResol('OEI', handles.modelo, handles.der.matrix, handles.rhoOEI, handles.betaMinOEI, handles.VminOEI, handles.PminOEI);


set(handles.phiOEI_val,'String', num2str(handles.phi*180/pi));
set(handles.daOEI_val,'String', num2str(handles.da*180/pi));
set(handles.drOEI_val,'String', num2str(handles.dr*180/pi));

guidata(hObject, handles);


function hOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to hOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hOEI_val as text
%        str2double(get(hObject,'String')) returns contents of hOEI_val as a double
handles.hOEI = str2double(get(hObject,'String'));

cd(handles.dir.code);
[~, handles.rhoOEI, ~, ~] = atmos_inter(handles.hOEI);  
cd(handles.dir.gui);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function hOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function betaMinOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to betaMinOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of betaMinOEI_val as text
%        str2double(get(hObject,'String')) returns contents of betaMinOEI_val as a double
handles.betaMinOEI = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function betaMinOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betaMinOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function betaMaxOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to betaMaxOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of betaMaxOEI_val as text
%        str2double(get(hObject,'String')) returns contents of betaMaxOEI_val as a double
handles.betaMaxOEI = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function betaMaxOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betaMaxOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function VminOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to VminOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VminOEI_val as text
%        str2double(get(hObject,'String')) returns contents of VminOEI_val as a double
handles.modelo.general.Vstall
handles.VminOEI = str2double(get(hObject,'String'))*handles.modelo.general.Vstall;
handles.VminOEI
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function VminOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VminOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VmaxOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to VmaxOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VmaxOEI_val as text
%        str2double(get(hObject,'String')) returns contents of VmaxOEI_val as a double

handles.VMaxOEI = str2double(get(hObject,'String'))*handles.modelo.general.Vstall;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function VmaxOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VmaxOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PminOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to PminOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PminOEI_val as text
%        str2double(get(hObject,'String')) returns contents of PminOEI_val as a double
handles.PminOEI = str2double(get(hObject,'String'))*handles.modelo.propulsion.Pmax;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PminOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PminOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PmaxOEI_val_Callback(hObject, eventdata, handles)
% hObject    handle to PmaxOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PmaxOEI_val as text
%        str2double(get(hObject,'String')) returns contents of PmaxOEI_val as a double
handles.PMaxOEI = str2double(get(hObject,'String'))*handles.modelo.propulsion.Pmax;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PmaxOEI_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PmaxOEI_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in barridoViento_checkbox.
function barridoViento_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to barridoViento_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of barridoViento_checkbox
set(hObject,'Value',0);
switch get(hObject,'Value')    
    case 1
        set(handles.VMaxViento_val,'Visible', 'on');
        set(handles.betaMaxViento_val,'Visible', 'on');
        set(handles.MinViento_text,'Visible', 'on');
        set(handles.MaxViento_text,'Visible', 'on');
        set(handles.Vientobarrido_button,'Visible', 'on');
    case 0
        set(handles.VMaxViento_val,'Visible', 'off');
        set(handles.betaMaxViento_val,'Visible', 'off');
        set(handles.MinViento_text,'Visible', 'off');
        set(handles.MaxViento_text,'Visible', 'off');
        set(handles.Vientobarrido_button,'Visible', 'off');
end



% --- Executes on button press in VientoSystem_button.
function VientoSystem_button_Callback(hObject, eventdata, handles)
% hObject    handle to VientoSystem_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.phi, handles.da, handles.dr] = OEIsystemResol('Viento', handles.modelo, handles.der.matrix, handles.rhoViento, handles.betaMinViento, handles.VminViento);

set(handles.phiViento_val,'String', num2str(handles.phi*180/pi));
set(handles.daViento_val,'String', num2str(handles.da*180/pi));
set(handles.drViento_val,'String', num2str(handles.dr*180/pi));

guidata(hObject, handles);


function hViento_val_Callback(hObject, eventdata, handles)
% hObject    handle to hViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hViento_val as text
%        str2double(get(hObject,'String')) returns contents of hViento_val as a double

handles.hViento = str2double(get(hObject,'String'));

cd(handles.dir.code);
[~, handles.rhoViento, ~, ~] = atmos_inter(handles.hViento);  
cd(handles.dir.gui);

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function hViento_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Vientobarrido_button.
function Vientobarrido_button_Callback(hObject, eventdata, handles)
% hObject    handle to Vientobarrido_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nBeta   = 5;
handles.graphColor      = hsv(handles.nBeta);
handles.beta    = linspace(handles.betaMinOEI, handles.betaMaxOEI, handles.nBeta);

switch handles.ejeX
    case 'V'
        handles.V       = handles.VminOEI:0.01:handles.VmaxOEI;
        handles.P       = handles.PminOEI;
    case 'P'
        handles.P       = handles.PinOEI:0.01:handles.PmaxOEI;
        handles.V       = handles.VminOEI;
end

[handles.phi, handles.da, handles.dr] = OEIsystemResol('OEI', handles.modelo, handles.der.matrix, handles.rhoOEI, handles.beta, handles.V, handles.P);

[handles.fileName, handles.pathName] = uiputfile('*.png', 'Guardar resultados como', '../');
phiFigure   = figure();
hold on;
grid on;

for jBeta = 1:handles.nBeta
    plot(V, handles.phi(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
    daLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi,3), 'º'];
end

xlabel('V');
ylabel('\phi');
legend(phiLegend);

set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print('-dpng','-r300',[handles.pathName,'phi_', handles.fileName]);
close(phiFigure);

daFigure    = figure();
hold on;
grid on;

for jBeta = 1:handles.nBeta
    plot(V, handles.da(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
    daLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi,3), 'º'];
end

xlabel('V');
ylabel('\delta_a');
legend(daLegend);

legend(daLegend)
set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print('-dpng','-r300',[handles.pathName,'da_', handles.fileName]);
close(daFigure);

drFigure    = figure();
hold on;
grid on;

for jBeta = 1:handles.nBeta
    plot(V, handles.dr(:,jBeta)*180/pi, handles.graphColor(jBeta,:));
    drLegend{jBeta} = ['\beta = ', num2str(beta(jBeta)*180/pi,2), 'º'];
end

xlabel('V');
ylabel('\delta_r');
legend(drLegend);

set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print('-dpng','-r300',[handles.pathName,'dr_', handles.fileName]);
close(drFigure);

winopen(handles.pathName);
guidata(hObject, handles);


function betaMinViento_val_Callback(hObject, eventdata, handles)
% hObject    handle to betaMinViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of betaMinViento_val as text
%        str2double(get(hObject,'String')) returns contents of betaMinViento_val as a double
handles.betaMinViento = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function betaMinViento_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betaMinViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function betaMaxViento_val_Callback(hObject, eventdata, handles)
% hObject    handle to betaMaxViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of betaMaxViento_val as text
%        str2double(get(hObject,'String')) returns contents of betaMaxViento_val as a double
handles.betaMaxViento = str2double(get(hObject,'String'))*pi/180;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function betaMaxViento_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to betaMaxViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VMinViento_val_Callback(hObject, eventdata, handles)
% hObject    handle to VminViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VminViento_val as text
%        str2double(get(hObject,'String')) returns contents of VminViento_val as a double
handles.VminViento = str2double(get(hObject,'String'))*handles.modelo.general.Vstall;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function VMinViento_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VMinViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VMaxViento_val_Callback(hObject, eventdata, handles)
% hObject    handle to VMaxViento_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VMaxViento_val as text
%        str2double(get(hObject,'String')) returns contents of VMaxViento_val as a double
handles.VmaxViento = str2double(get(hObject,'String'))*handles.modelo.general.Vstall;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function VMaxViento_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VMaxViento_val (see GCBO)
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
if isequal(get(handles.figure1,'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    guidata(hObject, handles);
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end
guidata(hObject,handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1,'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    guidata(hObject, handles);
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end
guidata(hObject,handles);


% --- Executes on button press in editButton.
function editButton_Callback(hObject, eventdata, handles)
% hObject    handle to editButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.modelo,handles.changes] = editMainMenu(handles.modelo, handles.dir);
handles = getDerivatives(handles);
guidata(hObject, handles);


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
modelo = handles.modelo;
[modelName,saveBol] = preLongAdj_selectName(modelo.name);
if saveBol
    fileName = [handles.dir.model,modelName,'.mat'];
    save(fileName, 'modelo');
end
guidata(hObject,handles)
