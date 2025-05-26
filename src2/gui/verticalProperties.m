function varargout = verticalProperties(varargin)
%VERTICALPROPERTIES M-file for verticalProperties.fig
%      VERTICALPROPERTIES, by itself, creates a new VERTICALPROPERTIES or raises the existing
%      singleton*.
%
%      H = VERTICALPROPERTIES returns the handle to a new VERTICALPROPERTIES or the handle to
%      the existing singleton*.
%
%      VERTICALPROPERTIES('Property','Value',...) creates a new VERTICALPROPERTIES using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to verticalProperties_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      VERTICALPROPERTIES('CALLBACK') and VERTICALPROPERTIES('CALLBACK',hObject,...) call the
%      local function named CALLBACK in VERTICALPROPERTIES.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help verticalProperties

% Last Modified by GUIDE v2.5 17-Mar-2016 10:52:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @verticalProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @verticalProperties_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before verticalProperties is made visible.
function verticalProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for verticalProperties
modelo                      = varargin{1};
handles.output              = modelo.vertical;
handles.vertical            = handles.output;
handles.dir                 = varargin{2};
handles.changes             = varargin{3};

rotate3d(handles.axes1);
axis(handles.axes1,'equal');


set(handles.S_val,'String', num2str(handles.vertical.S,3));
set(handles.b_val,'String', num2str(handles.vertical.b,3));
set(handles.AR_val,'String', num2str(handles.vertical.AR,3));

set(handles.cr_val,'String', num2str(handles.vertical.cr,3));
set(handles.ct_val,'String', num2str(handles.vertical.ct,3));
set(handles.TR_val,'String', num2str(handles.vertical.TR,3));
set(handles.cMAC_val,'String', num2str(handles.vertical.MAC,3));
set(handles.tc_val,'String', num2str(handles.vertical.t_c,3));
set(handles.LAMle_val,'String', num2str(handles.vertical.LAM*180/pi,3));




set(handles.Xca_val,'String', num2str(handles.vertical.Xca,3));
set(handles.Zca_val,'String', num2str(handles.vertical.Zca,3));
set(handles.eta_val,'String', num2str(handles.vertical.eta,3));

set(handles.CLa_val,'String', num2str(handles.vertical.CLa,3));
set(handles.Cla2_val,'String', num2str(handles.vertical.Cla,3));


set(handles.cac_val,'String', num2str(handles.vertical.cm_c,3));
set(handles.y0b2_val,'String', num2str(handles.vertical.y0_b2,3));
set(handles.y1b2_val,'String', num2str(handles.vertical.y1_b2,3));

handles.chordDistr = [handles.vertical.z/handles.vertical.b,... 
                     handles.vertical.le_z/handles.vertical.cr,...
                     handles.vertical.c_z/handles.vertical.cr]; 
set(handles.MATfile_val,'String', num2str(handles.chordDistr, 3));

if ~isempty(handles.chordDistr)                   
    handles = refreshVertdata(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wingProperties wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = verticalProperties_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.changes;
delete(hObject);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isequaln(handles.output, handles.vertical)
    handles.changes = 1;
end

handles.output = handles.vertical;

guidata(hObject, handles);


% --- Executes on button press in loadMAT_button.
function loadMAT_button_Callback(hObject, eventdata, handles)
% hObject    handle to loadMAT_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[matFileName,matPathName,~] = uigetfile('*.mat','Seleccione los datos de la superficie aerodinámica');

if ~isequal(matPathName, 0) && ~isequal(matFileName, 0)
    
    handles.chordDistr = load([matPathName, matFileName]);
    obj_name = fieldnames(handles.chordDistr);
    handles.chordDistr = handles.chordDistr.(obj_name{1});

    if size(handles.chordDistr,2) ~= 3 || size(handles.chordDistr,1) < 2
        valid = 0;
        handles = []; % rmfield(handles,'chordDistr');
    else
        valid = 1;
        set(handles.MATfile_val,'String',num2str(handles.chordDistr));
        
    end

    if valid && ~isempty(handles.vertical.b) && ~isempty(handles.vertical.S)
        handles = refreshVertdata(handles);
    end
end
guidata(hObject, handles);


function MATfile_val_Callback(hObject, eventdata, handles)
% hObject    handle to MATfile_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MATfile_val as text
%        str2double(get(hObject,'String')) returns contents of MATfile_val as a double
handles.chordDistr = str2num(get(handles.MATfile_val,'String'));

if isempty(get(handles.MATfile_val,'String')) || (size(handles.chordDistr,2) ~= 3) || (size(handles.chordDistr,1) < 2)
    valid = 0;
    handles.chordDistr = []; %rmfield(handles,'chordDistr');
else
    valid = 1;
    
end

if valid && ~isempty(handles.vertical.b) && ~isempty(handles.vertical.S)   
    handles = refreshVertdata(handles);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MATfile_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MATfile_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function cac_val_Callback(hObject, eventdata, handles)
% hObject    handle to cac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cac_val as text
%        str2double(get(hObject,'String')) returns contents of cac_val as a double
handles.vertical.cm_c = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cac_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y0b2_val_Callback(hObject, eventdata, handles)
% hObject    handle to y0b2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y0b2_val as text
%        str2double(get(hObject,'String')) returns contents of y0b2_val as a double
handles.vertical.y0_b2 = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function y0b2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0b2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y1b2_val_Callback(hObject, eventdata, handles)
% hObject    handle to y1b2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y1b2_val as text
%        str2double(get(hObject,'String')) returns contents of y1b2_val as a double
handles.vertical.y1_b2 = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function y1b2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y1b2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CL0_val_Callback(hObject, eventdata, handles)
% hObject    handle to CL0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CL0_val as text
%        str2double(get(hObject,'String')) returns contents of CL0_val as a double
handles.vertical.CL0 = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CL0_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CL0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CLa_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLa_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLa_val as text
%        str2double(get(hObject,'String')) returns contents of CLa_val as a double
handles.vertical.CLa = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CLa_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLa_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cla2_val_Callback(hObject, eventdata, handles)
% hObject    handle to Cla2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cla2_val as text
%        str2double(get(hObject,'String')) returns contents of Cla2_val as a double
handles.vertical.Cla = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Cla2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cla2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CM0_val_Callback(hObject, eventdata, handles)
% hObject    handle to CM0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CM0_val as text
%        str2double(get(hObject,'String')) returns contents of CM0_val as a double
handles.vertical.CM0 = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CM0_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CM0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eta_val_Callback(hObject, eventdata, handles)
% hObject    handle to eta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eta_val as text
%        str2double(get(hObject,'String')) returns contents of eta_val as a double
handles.vertical.eta = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function eta_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eta_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diedro_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xca_val as text
%        str2double(get(hObject,'String')) returns contents of Xca_val as a double
handles.vertical.diedro = str2double(get(hObject,'String'))*pi/180;
handles.vertical.zca     = handles.vertical.yca*sin(handles.vertical.diedro);
set(handles.zca2_val, 'String', handles.vertical.zca);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function diedro_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMc4_val_Callback(hObject, eventdata, handles)
% hObject    handle to yca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yca2_val as text
%        str2double(get(hObject,'String')) returns contents of yca2_val as a double
handles.vertical.LAMc4 = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LAMc4_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMc2_val_Callback(hObject, eventdata, handles)
% hObject    handle to zca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zca2_val as text
%        str2double(get(hObject,'String')) returns contents of zca2_val as a double
handles.vertical.LAMc2 = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LAMc2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zca_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMle_val as text
%        str2double(get(hObject,'String')) returns contents of LAMle_val as a double
handles.vertical.Zca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Zca_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cMAC_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xca_val as text
%        str2double(get(hObject,'String')) returns contents of Xca_val as a double
handles.vertical.MAC = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cMAC_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xca2_val_Callback(hObject, eventdata, handles)
% hObject    handle to xca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xca2_val as text
%        str2double(get(hObject,'String')) returns contents of xca2_val as a double
handles.vertical.xca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xca2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yca2_val_Callback(hObject, eventdata, handles)
% hObject    handle to yca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yca2_val as text
%        str2double(get(hObject,'String')) returns contents of yca2_val as a double
handles.vertical.yca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function yca2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zca2_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMc2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMc2_val as text
%        str2double(get(hObject,'String')) returns contents of LAMc2_val as a double
handles.vertical.zca = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function zca2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMc2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_val_Callback(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zca_val as text
%        str2double(get(hObject,'String')) returns contents of Zca_val as a double
handles.vertical.i = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function i_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tc_val_Callback(hObject, eventdata, handles)
% hObject    handle to zca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zca2_val as text
%        str2double(get(hObject,'String')) returns contents of zca2_val as a double
handles.vertical.t_c = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tc_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sexp_val_Callback(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zca_val as text
%        str2double(get(hObject,'String')) returns contents of Zca_val as a double
handles.vertical.Se = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Sexp_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
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
handles.vertical.S = str2double(get(handles.S_val,'String'));
if ~isempty(handles.vertical.b)
    handles.vertical.AR = handles.vertical.b^2/handles.vertical.S;
    set(handles.AR_val,'String',num2str(handles.vertical.AR));
elseif ~isempty(handles.vertical.AR)
    handles.vertical.b = sqrt(handles.vertical.AR*handles.vertical.S);
    set(handles.b_val,'String',num2str(handles.vertical.b));
end

if ~isempty(handles.chordDistr) && ~isempty(handles.vertical.b) 
    handles = refreshVertdata(handles);
end

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



function b_val_Callback(hObject, eventdata, handles)
% hObject    handle to b_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_val as text
%        str2double(get(hObject,'String')) returns contents of b_val as a double
handles.vertical.b = str2double(get(handles.b_val,'String'));

if ~isempty(handles.vertical.S) && isempty(handles.vertical.AR)
    handles.vertical.AR = handles.vertical.b^2/handles.vertical.S;
    set(handles.AR_val,'String',num2str(handles.vertical.AR));
elseif ~isempty(handles.vertical.AR) && isempty(handles.vertical.S)
    handles.vertical.S = handles.vertical.b^2*handles.vertical.AR;
    set(handles.S_val,'String',num2str(handles.vertical.S));
end

if ~isempty(handles.chordDistr) && ~isempty(handles.vertical.S)
    handles = refreshVertdata(handles);    
end

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



function AR_val_Callback(hObject, eventdata, handles)
% hObject    handle to AR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AR_val as text
%        str2double(get(hObject,'String')) returns contents of AR_val as a double
handles.vertical.AR = str2double(get(hObject,'String'));

if ~isempty(handles.vertical.S) && isempty(handles.vertical.b)
    handles.vertical.b = sqrt(handles.vertical.S*handles.vertical.AR); 
    set(handles.bval,'String',num2str(handles.vertical.b));
elseif isempty(handles.vertical.S) && ~isempty(handles.vertical.b)
    handles.vertical.S = handles.vertical.b^2/handles.vertical.AR; 
    set(handles.Sval,'String',num2str(handles.vertical.S));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function AR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = refreshVertdata(handles)
    handles.vertical.MAC = handles.vertical.S/handles.vertical.b;
    cd(handles.dir.code);
    [handles.vertical.cr, handles.vertical.ct, handles.vertical.xca, handles.vertical.zca,...
    ~, ~, handles.vertical.LAM] =aeroSurf_geom(2*handles.vertical.S, 2*handles.vertical.b,...
    handles.chordDistr(:,1), handles.chordDistr(:,3), handles.chordDistr(:,2));
    cd(handles.dir.gui);
    
    handles.vertical.TR      = handles.vertical.ct/handles.vertical.cr;
    
    handles.vertical.le_z    = handles.chordDistr(:,2)*handles.vertical.cr;    % Coordenadas x del borde de ataque de vertical
    handles.vertical.c_z     = handles.chordDistr(:,3)*handles.vertical.cr;    % Cuerdas de las secciones del vertical.
    handles.vertical.z       = handles.chordDistr(:,1)*handles.vertical.b;   % Coordenadas y de las secciones del vertical
       
    cla(handles.axes1);
    hold on; grid on;
    if isempty(handles.vertical.t_c)
        t_c = 0.15;
    else
        t_c = handles.vertical.t_c;
    end
    cd(handles.dir.code);
    [x,y]   = getNACA(t_c);
    cd(handles.dir.gui);
    n       = length(x);
    x_mesh  = handles.vertical.c_z*x + handles.vertical.le_z*ones(1,n);
    y_mesh  = handles.vertical.c_z*y;
    z_mesh  = handles.vertical.z*ones(1,n);

    mesh(x_mesh,y_mesh,z_mesh);
    grid on;
%     y_plot  = [-handles.vertical.z(end:-1:1); handles.vertical.z];
%     le_plot = [handles.vertical.le_z(end:-1:1); handles.vertical.le_z];
%     te_plot = le_plot + [handles.vertical.c_z(end:-1:1); handles.vertical.c_z];
%     
%     plot(y_plot, le_plot, 'b');
%     plot(y_plot, te_plot, 'b');
%     plot(handles.vertical.yca, handles.vertical.xca, 'rx');
    
    set(handles.cr_val,'String', num2str(handles.vertical.cr,3));
    set(handles.ct_val,'String', num2str(handles.vertical.ct,3));
    set(handles.TR_val,'String', num2str(handles.vertical.TR,3));
    set(handles.cMAC_val,'String', num2str(handles.vertical.MAC,3));
    set(handles.LAMle_val,'String', num2str(handles.vertical.LAM*180/pi,2));


function TR_val_Callback(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zca_val as text
%        str2double(get(hObject,'String')) returns contents of Zca_val as a double
if ~isempty(handles.vertical.cr) && isempty(handles.vertical.ct)
    handles.vertical.ct = handles.vertical.cr*handles.vertical.TR; 
    set(handles.ct_val,'String',num2str(handles.vertical.ct));
elseif ~isempty(handles.vertical.ct) && isempty(handles.vertical.cr)
    handles.vertical.cr = handles.vertical.ct/handles.vertical.TR; 
    set(handles.cr_val,'String',num2str(handles.vertical.cr));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cr_val_Callback(hObject, eventdata, handles)
% hObject    handle to cr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr_val as text
%        str2double(get(hObject,'String')) returns contents of cr_val as a double
handles.vertical.cr = str2double(get(handles.cr_val,'String'));

if ~isempty(handles.vertical.ct) && isempty(handles.vertical.TR)
    handles.vertical.TR = handles.vertical.ct/handles.vertical.cr;
    set(handles.TR_val,'String',num2str(handles.vertical.TR));
elseif ~isempty(handles.vertical.TR) && isempty(handles.vertical.ct)
    handles.vertical.ct = handles.vertical.cr*handles.vertical.TR;
    set(handles.ct_val,'String',num2str(handles.vertical.ct));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ct_val_Callback(hObject, eventdata, handles)
% hObject    handle to ct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ct_val as text
%        str2double(get(hObject,'String')) returns contents of ct_val as a double
handles.vertical.ct = str2double(get(handles.ct_val,'String'));

if ~isempty(handles.vertical.cr) && isempty(handles.vertical.TR)
    handles.vertical.TR = handles.vertical.ct/handles.vertical.cr;
    set(handles.TR_val,'String',num2str(handles.vertical.TR));
elseif ~isempty(handles.vertical.TR) && isempty(handles.vertical.cr)
    handles.vertical.cr = handles.vertical.ct*handles.vertical.TR;
    set(handles.cr_val,'String',num2str(handles.vertical.cr));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ct_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMle_val_Callback(hObject, eventdata, handles)
% hObject    handle to TR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR_val as text
%        str2double(get(hObject,'String')) returns contents of TR_val as a double
handles.vertical.LAM = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LAMle_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Xca_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xca_val as text
%        str2double(get(hObject,'String')) returns contents of Xca_val as a double
handles.vertical.Xca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xca_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xca_val (see GCBO)
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
if ~isequaln(handles.output,handles.vertical)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.vertical;
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

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if ~isequaln(handles.output,handles.vertical)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.vertical;
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



% --- Executes during object creation, after setting all properties.
function verticalType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verticalType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen([handles.dir.help,'vtp_guide.pdf']);
