function varargout = wingProperties(varargin)
% WINGPROPERTIES MATLAB code for wingProperties.fig
%      WINGPROPERTIES, by itself, creates a new WINGPROPERTIES or raises the existing
%      singleton*.
%
%      H = WINGPROPERTIES returns the handle to a new WINGPROPERTIES or the handle to
%      the existing singleton*.
%
%      WINGPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINGPROPERTIES.M with the given input arguments.
%
%      WINGPROPERTIES('Property','Value',...) creates a new WINGPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wingProperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wingProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wingProperties

% Last Modified by GUIDE v2.5 02-Jan-2016 13:32:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wingProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @wingProperties_OutputFcn, ...
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


% --- Executes just before wingProperties is made visible.
function wingProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wingProperties (see VARARGIN)

% Choose default command line output for wingProperties
modelo          = varargin{1};
handles.output  = modelo.ala;
handles.ala     = handles.output;

handles.dir     = varargin{2};

axis(handles.axes1,'equal');
rotate3d(handles.axes1);

handles.changes = varargin{3};


set(handles.S_val,'String', num2str(handles.ala.S,3));
set(handles.b_val,'String', num2str(handles.ala.b,3));
set(handles.AR_val,'String', num2str(handles.ala.AR,3));
set(handles.TR_val,'String', num2str(handles.ala.TR,3));
set(handles.cr_val,'String', num2str(handles.ala.cr,3));
set(handles.ct_val,'String', num2str(handles.ala.ct,3));
set(handles.LAMle_val,'String', num2str(handles.ala.LAM*180/pi,3));
set(handles.Xca_val,'String', num2str(handles.ala.Xca,3));
set(handles.Zca_val,'String', num2str(handles.ala.Zca,3));
set(handles.diedro_val,'String', num2str(handles.ala.diedro*180/pi,3));
set(handles.LAMc4_val,'String', num2str(handles.ala.LAMc4*180/pi,3));
set(handles.LAMc2_val,'String', num2str(handles.ala.LAMc2*180/pi,3));
set(handles.cMAC_val,'String', num2str(handles.ala.MAC,3));
set(handles.xca2_val,'String', num2str(handles.ala.xca,3));
set(handles.yca2_val,'String', num2str(handles.ala.yca,3));
set(handles.zca2_val,'String', num2str(handles.ala.zca,3));
set(handles.i_val,'String', num2str(handles.ala.i*180/pi,3));
set(handles.tc_val,'String', num2str(handles.ala.t_c,3));
set(handles.Sexp_val,'String', num2str(handles.ala.Se,3));
set(handles.CL0_val,'String', num2str(handles.ala.CL0,3));
set(handles.CM0_val,'String', num2str(handles.ala.CM0,3));
set(handles.CLa_val,'String', num2str(handles.ala.CLa,3));
set(handles.Cla2_val,'String', num2str(handles.ala.Cla,3));
set(handles.eta_val,'String', num2str(handles.ala.eta,3));
set(handles.cac_val,'String', num2str(handles.ala.cm_c,3));
set(handles.y0b2_val,'String', num2str(handles.ala.y0_b2,3));
set(handles.y1b2_val,'String', num2str(handles.ala.y1_b2,3));
handles.chordDistr =    [handles.ala.y/handles.ala.b*2,... 
                        handles.ala.le_y/handles.ala.cr,...
                        handles.ala.c_y/handles.ala.cr];
                    
set(handles.MATfile_val,'String', num2str(handles.chordDistr,3));
                              
if isempty(get(handles.MATfile_val,'String')) || isempty(get(handles.S_val,'String')) || isempty(get(handles.b_val,'String')) 
else
    handles = refreshSurfdata(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wingProperties wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wingProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.changes;
delete(hObject);


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
        handles.chordDistr = []; % rmfield(handles,'chordDistr');
    else
        valid = 1;
        set(handles.MATfile_val,'String',num2str(handles.chordDistr));
    end

    if valid && ~isempty(handles.ala.b) && ~isempty(handles.ala.S)
        handles = refreshSurfdata(handles);
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

if valid && ~isempty(handles.ala.b) && ~isempty(handles.ala.S)   
    handles = refreshSurfdata(handles);
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

function handles = refreshSurfdata(handles)
    handles.ala.MAC = handles.ala.S/handles.ala.b;
    cd(handles.dir.code);
    [handles.ala.cr, handles.ala.ct, handles.ala.xca, handles.ala.yca,...
    handles.ala.LAMc4, handles.ala.LAMc2, handles.ala.LAM] =...
    aeroSurf_geom(handles.ala.S, handles.ala.b, handles.chordDistr(:,1), handles.chordDistr(:,3), handles.chordDistr(:,2));
    cd(handles.dir.gui);
    
    handles.ala.TR      = handles.ala.ct/handles.ala.cr;
    handles.ala.zca     = handles.ala.yca*sin(handles.ala.diedro);
    
    handles.ala.le_y    = handles.chordDistr(:,2)*handles.ala.cr;    % Coordenadas x del borde de ataque de ala
    handles.ala.c_y     = handles.chordDistr(:,3)*handles.ala.cr;    % Cuerdas de las secciones del ala.
    handles.ala.y       = handles.chordDistr(:,1)*handles.ala.b/2;   % Coordenadas y de las secciones del ala
       
%     cla(handles.axes1);
%     axis(handles.axes1,'equal');
%     hold on; grid on;
    
    cla(handles.axes1);
    hold on; grid on;
    if isempty(handles.ala.t_c)
        t_c = 0.15;
    else
        t_c = handles.ala.t_c;
    end
    cd(handles.dir.code);
    [x,z]   = getNACA(t_c);
    cd(handles.dir.gui);
    n       = length(x);
    x_mesh  = handles.ala.c_y*x + handles.ala.le_y*ones(1,n);
    z_mesh  = handles.ala.c_y*z;
    y_mesh  = handles.ala.y*ones(1,n);
    
    x_mesh  = [x_mesh(end:(-1):1,:); x_mesh];
    y_mesh  = [-y_mesh(end:(-1):1,:); y_mesh];
    z_mesh  = [z_mesh(end:(-1):1,:); z_mesh];

    mesh(y_mesh,x_mesh,z_mesh);
    
    
    size = 30;
    color = [1 0 0];
    scatter3(0, handles.ala.xca, handles.ala.cr*t_c*0.6,size,color,'fill')
    grid on;
    
   
%     y_plot  = [-handles.ala.y(end:-1:1); handles.ala.y];
%     le_plot = [handles.ala.le_y(end:-1:1); handles.ala.le_y];
%     te_plot = le_plot + [handles.ala.c_y(end:-1:1); handles.ala.c_y];
%     
%     plot(y_plot, le_plot, 'b');
%     plot(y_plot, te_plot, 'b');
%     plot(handles.ala.yca, handles.ala.xca, 'rx');
    
    set(handles.cr_val,'String', num2str(handles.ala.cr,3));
    set(handles.ct_val,'String', num2str(handles.ala.ct,3));
    set(handles.TR_val,'String', num2str(handles.ala.TR,3));
    set(handles.cMAC_val,'String', num2str(handles.ala.MAC,3));
    set(handles.xca2_val,'String', num2str(handles.ala.xca,3));
    set(handles.yca2_val,'String', num2str(handles.ala.yca,3));
    set(handles.zca2_val,'String', num2str(handles.ala.zca,3));
    set(handles.LAMle_val,'String', num2str(handles.ala.LAM*180/pi,2));
    set(handles.LAMc4_val,'String', num2str(handles.ala.LAMc4*180/pi,2));
    set(handles.LAMc2_val,'String', num2str(handles.ala.LAMc2*180/pi,2));




function CL0_val_Callback(hObject, eventdata, handles)
% hObject    handle to CL0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CL0_val as text
%        str2double(get(hObject,'String')) returns contents of CL0_val as a double
handles.ala.CL0 = str2double(get(hObject,'String'));

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
handles.ala.CLa = str2double(get(hObject,'String'));
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
handles.ala.Cla = str2double(get(hObject,'String'));

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
handles.ala.CM0 = str2double(get(hObject,'String'));

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
handles.ala.eta = str2double(get(hObject,'String'));

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


function Xca_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMle_val as text
%        str2double(get(hObject,'String')) returns contents of LAMle_val as a double
handles.ala.Xca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xca_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AR_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMle_val as text
%        str2double(get(hObject,'String')) returns contents of LAMle_val as a double
handles.ala.AR = str2double(get(hObject,'String'));

if ~isempty(handles.ala.S) && isempty(handles.ala.b)
    handles.ala.b = sqrt(handles.ala.S*handles.ala.AR); 
    set(handles.bval,'String',num2str(handles.ala.b));
elseif isempty(handles.ala.S) && ~isempty(handles.ala.b)
    handles.ala.S = handles.ala.b^2/handles.ala.AR; 
    set(handles.Sval,'String',num2str(handles.ala.S));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function AR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMle_val as text
%        str2double(get(hObject,'String')) returns contents of LAMle_val as a double
handles.ala.S = str2double(get(handles.S_val,'String'));
if ~isempty(handles.ala.b)
    handles.ala.AR = handles.ala.b^2/handles.ala.S;
    set(handles.AR_val,'String',num2str(handles.ala.AR));
elseif ~isempty(handles.ala.AR)
    handles.ala.b = sqrt(handles.ala.AR*handles.ala.S);
    set(handles.b_val,'String',num2str(handles.ala.b));
end

if ~isempty(handles.chordDistr) && ~isempty(handles.ala.b) 
    handles = refreshSurfdata(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function S_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function b_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMle_val as text
%        str2double(get(hObject,'String')) returns contents of LAMle_val as a double
handles.ala.b = str2double(get(handles.b_val,'String'));

if ~isempty(handles.ala.S) && isempty(handles.ala.AR)
    handles.ala.AR = handles.ala.b^2/handles.ala.S;
    set(handles.AR_val,'String',num2str(handles.ala.AR));
elseif ~isempty(handles.ala.AR) && isempty(handles.ala.S)
    handles.ala.S = handles.ala.b^2*handles.ala.AR;
    set(handles.S_val,'String',num2str(handles.ala.S));
end

if ~isempty(handles.chordDistr) && ~isempty(handles.ala.S)
    handles = refreshSurfdata(handles);    
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function b_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cr_val_Callback(hObject, eventdata, handles)
% hObject    handle to xca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xca2_val as text
%        str2double(get(hObject,'String')) returns contents of xca2_val as a double
handles.ala.cr = str2double(get(handles.cr_val,'String'));

if ~isempty(handles.ala.ct) && isempty(handles.ala.TR)
    handles.ala.TR = handles.ala.ct/handles.ala.cr;
    set(handles.TR_val,'String',num2str(handles.ala.TR));
elseif ~isempty(handles.ala.TR) && isempty(handles.ala.ct)
    handles.ala.ct = handles.ala.cr*handles.ala.TR;
    set(handles.ct_val,'String',num2str(handles.ala.ct));
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ct_val_Callback(hObject, eventdata, handles)
% hObject    handle to cr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr_val as text
%        str2double(get(hObject,'String')) returns contents of cr_val as a double
handles.ala.ct = str2double(get(handles.ct_val,'String'));

if ~isempty(handles.ala.cr) && isempty(handles.ala.TR)
    handles.ala.TR = handles.ala.ct/handles.ala.cr;
    set(handles.TR_val,'String',num2str(handles.ala.TR));
elseif ~isempty(handles.ala.TR) && isempty(handles.ala.cr)
    handles.ala.cr = handles.ala.ct*handles.ala.TR;
    set(handles.cr_val,'String',num2str(handles.ala.cr));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ct_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TR_val_Callback(hObject, eventdata, handles)
% hObject    handle to xca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xca2_val as text
%        str2double(get(hObject,'String')) returns contents of xca2_val as a double
if ~isempty(handles.ala.cr) && isempty(handles.ala.ct)
    handles.ala.ct = handles.ala.cr*handles.ala.TR; 
    set(handles.ct_val,'String',num2str(handles.ala.ct));
elseif ~isempty(handles.ala.ct) && isempty(handles.ala.cr)
    handles.ala.cr = handles.ala.ct/handles.ala.TR; 
    set(handles.cr_val,'String',num2str(handles.ala.cr));
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function TR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMle_val_Callback(hObject, eventdata, handles)
% hObject    handle to ct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ct_val as text
%        str2double(get(hObject,'String')) returns contents of ct_val as a double
handles.ala.LAM = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function LAMle_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ct_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diedro_val_Callback(hObject, eventdata, handles)
% hObject    handle to diedro_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diedro_val as text
%        str2double(get(hObject,'String')) returns contents of diedro_val as a double
handles.ala.diedro  = str2double(get(hObject,'String'))*pi/180;
handles.ala.zca     = handles.ala.yca*sin(handles.ala.diedro);
set(handles.zca2_val, 'String', handles.ala.zca);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function diedro_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diedro_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMc4_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMc4_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMc4_val as text
%        str2double(get(hObject,'String')) returns contents of LAMc4_val as a double
handles.ala.LAMc4 = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LAMc4_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMc4_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LAMc2_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMc2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMc2_val as text
%        str2double(get(hObject,'String')) returns contents of LAMc2_val as a double
handles.ala.LAMc2 = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LAMc2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMc2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zca_val_Callback(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zca_val as text
%        str2double(get(hObject,'String')) returns contents of Zca_val as a double
handles.ala.Zca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Zca_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cMAC_val_Callback(hObject, eventdata, handles)
% hObject    handle to cMAC_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cMAC_val as text
%        str2double(get(hObject,'String')) returns contents of cMAC_val as a double
handles.ala.MAC = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cMAC_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cMAC_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xca2_val_Callback(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LAMle_val as text
%        str2double(get(hObject,'String')) returns contents of LAMle_val as a double
handles.ala.xca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xca2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LAMle_val (see GCBO)
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
handles.ala.yca = str2double(get(hObject,'String'));

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



function cac_val_Callback(hObject, eventdata, handles)
% hObject    handle to cac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cac_val as text
%        str2double(get(hObject,'String')) returns contents of cac_val as a double
handles.ala.cm_c = str2double(get(hObject,'String'));

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
handles.ala.y0_b2 = str2double(get(hObject,'String'));

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
handles.ala.y1_b2 = str2double(get(hObject,'String'));

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



function zca2_val_Callback(hObject, eventdata, handles)
% hObject    handle to zca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zca2_val as text
%        str2double(get(hObject,'String')) returns contents of zca2_val as a double
handles.ala.zca = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function zca2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zca2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_val_Callback(hObject, eventdata, handles)
% hObject    handle to i_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_val as text
%        str2double(get(hObject,'String')) returns contents of i_val as a double
handles.ala.i = str2double(get(hObject,'String'))*pi/180;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function i_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tc_val_Callback(hObject, eventdata, handles)
% hObject    handle to tc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tc_val as text
%        str2double(get(hObject,'String')) returns contents of tc_val as a double
handles.ala.t_c = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tc_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sexp_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sexp_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sexp_val as text
%        str2double(get(hObject,'String')) returns contents of Sexp_val as a double
handles.ala.Se = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Sexp_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sexp_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isequaln(handles.output, handles.ala)
    handles.changes = 1;
end

handles.output = handles.ala;
guidata(hObject, handles);

% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1,'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    if ~isequaln(handles.output,handles.ala)
        saveAnswer = unsavedData_dialog(handles.dir);
        switch saveAnswer
            case 'SAVE'
                handles.output = handles.ala;
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
    if ~isequaln(handles.output,handles.ala)
        saveAnswer = unsavedData_dialog(handles.dir);
        switch saveAnswer
            case 'SAVE'
                handles.output = handles.ala;
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
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen([handles.dir.help,'wing_guide.pdf']);
