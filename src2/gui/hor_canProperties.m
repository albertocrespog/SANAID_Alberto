function varargout = hor_canProperties(varargin)
%HOR_CANPROPERTIES M-file for hor_canProperties.fig
%      HOR_CANPROPERTIES, by itself, creates a new HOR_CANPROPERTIES or raises the existing
%      singleton*.
%
%      H = HOR_CANPROPERTIES returns the handle to a new HOR_CANPROPERTIES or the handle to
%      the existing singleton*.
%
%      HOR_CANPROPERTIES('Property','Value',...) creates a new HOR_CANPROPERTIES using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to hor_canProperties_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      HOR_CANPROPERTIES('CALLBACK') and HOR_CANPROPERTIES('CALLBACK',hObject,...) call the
%      local function named CALLBACK in HOR_CANPROPERTIES.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hor_canProperties

% Last Modified by GUIDE v2.5 03-Sep-2015 17:43:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hor_canProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @hor_canProperties_OutputFcn, ...
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


% --- Executes just before hor_canProperties is made visible.
function hor_canProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for hor_canProperties


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



modelo                      = varargin{1};
handles.surfName            = varargin{2};
handles.output              = modelo.(handles.surfName);
handles.(handles.surfName)  = handles.output;
handles.dir                 = varargin{3};
handles.changes             = varargin{4};
axis(handles.axes1,'equal');
rotate3d(handles.axes1);

switch handles.surfName
    case 'canard'
        set(handles.windowsTitle, 'String', 'CANARD EDITION');
        set(hObject, 'Name', 'Canard Data Edition');
    case 'horizontal'
        set(handles.windowsTitle, 'String', 'HORIZONTAL STABILIZER EDITION');
        set(hObject, 'Name', 'Horizontal Data Edition');
end

handles.chordDistr = [];

set(handles.S_val,'String', num2str(handles.(handles.surfName).S,3));
set(handles.b_val,'String', num2str(handles.(handles.surfName).b,3));
set(handles.AR_val,'String', num2str(handles.(handles.surfName).AR,3));
set(handles.Sexp_val,'String', num2str(handles.(handles.surfName).Se,3));
set(handles.cr_val,'String', num2str(handles.(handles.surfName).cr,3));
set(handles.ct_val,'String', num2str(handles.(handles.surfName).ct,3));
set(handles.TR_val,'String', num2str(handles.(handles.surfName).TR,3));
set(handles.cMAC_val,'String', num2str(handles.(handles.surfName).MAC,3));
set(handles.tc_val,'String', num2str(handles.(handles.surfName).t_c,3));
set(handles.LAMle_val,'String', num2str(handles.(handles.surfName).LAM*180/pi,3));
set(handles.LAMc2_val,'String', num2str(handles.(handles.surfName).LAMc2*180/pi,3));
set(handles.LAMc4_val,'String', num2str(handles.(handles.surfName).LAMc4*180/pi,3));
set(handles.i_val,'String', num2str(handles.(handles.surfName).i*180/pi,3));
set(handles.diedro_val,'String', num2str(handles.(handles.surfName).diedro*180/pi,3));
set(handles.Xca_val,'String', num2str(handles.(handles.surfName).Xca,3));
set(handles.Zca_val,'String', num2str(handles.(handles.surfName).Zca,3));
set(handles.xca2_val,'String', num2str(handles.(handles.surfName).xca,3));
set(handles.yca2_val,'String', num2str(handles.(handles.surfName).yca,3));
set(handles.zca2_val,'String', num2str(handles.(handles.surfName).zca,3));

set(handles.CL0_val,'String', num2str(handles.(handles.surfName).CL0,3));
set(handles.CM0_val,'String', num2str(handles.(handles.surfName).CM0,3));
set(handles.CLa_val,'String', num2str(handles.(handles.surfName).CLa,3));
set(handles.Cla2_val,'String', num2str(handles.(handles.surfName).Cla,3));
set(handles.eta_val,'String', num2str(handles.(handles.surfName).eta,3));

set(handles.cac_val,'String', num2str(handles.(handles.surfName).cm_c,3));
set(handles.y0b2_val,'String', num2str(handles.(handles.surfName).y0_b2,3));
set(handles.y1b2_val,'String', num2str(handles.(handles.surfName).y1_b2,3));

handles.chordDistr =    [handles.(handles.surfName).y/handles.(handles.surfName).b*2,... 
                        handles.(handles.surfName).le_y/handles.(handles.surfName).cr,...
                        handles.(handles.surfName).c_y/handles.(handles.surfName).cr];
                    
set(handles.MATfile_val,'String', num2str(handles.chordDistr, 3));

if isempty(get(handles.MATfile_val,'String')) || isempty(get(handles.S_val,'String')) || isempty(get(handles.b_val,'String')) 
else
    handles = refreshSurfdata(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wingProperties wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = hor_canProperties_OutputFcn(hObject, eventdata, handles)
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
if ~isequaln(handles.output, handles.(handles.surfName))
    handles.changes = 1;
end
handles.output = handles.(handles.surfName);
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

    if valid && ~isempty(handles.(handles.surfName).b) && ~isempty(handles.(handles.surfName).S)
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

if valid && ~isempty(handles.(handles.surfName).b) && ~isempty(handles.(handles.surfName).S)   
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
    handles.(handles.surfName).MAC = handles.(handles.surfName).S/handles.(handles.surfName).b;
    cd(handles.dir.code);
    [handles.(handles.surfName).cr, handles.(handles.surfName).ct, handles.(handles.surfName).xca, handles.(handles.surfName).yca,...
    handles.(handles.surfName).LAMc4, handles.(handles.surfName).LAMc2, handles.(handles.surfName).LAM] =...
    aeroSurf_geom(handles.(handles.surfName).S, handles.(handles.surfName).b, handles.chordDistr(:,1), handles.chordDistr(:,3), handles.chordDistr(:,2));
    cd(handles.dir.gui);
    
    handles.(handles.surfName).TR       = handles.(handles.surfName).ct/handles.(handles.surfName).cr;
    handles.(handles.surfName).zca      = handles.(handles.surfName).yca*sin(handles.(handles.surfName).diedro);
    
    handles.(handles.surfName).le_y     = handles.chordDistr(:,2)*handles.(handles.surfName).cr;    % Coordenadas x del borde de ataque de ala
    handles.(handles.surfName).c_y      = handles.chordDistr(:,3)*handles.(handles.surfName).cr;    % Cuerdas de las secciones del ala.
    handles.(handles.surfName).y        = handles.chordDistr(:,1)*handles.(handles.surfName).b/2;   % Coordenadas y de las secciones del ala
       
    cla(handles.axes1);
    axis(handles.axes1,'equal');
    hold on; grid on;
    
     cla(handles.axes1);
    hold on; grid on;
    if isempty(handles.(handles.surfName).t_c)
        t_c = 0.15;
    else
        t_c = handles.(handles.surfName).t_c;
    end
    cd(handles.dir.code);
    [x,z]   = getNACA(t_c);
    cd(handles.dir.gui);
    n       = length(x);
    x_mesh  = handles.(handles.surfName).c_y*x + handles.(handles.surfName).le_y*ones(1,n);
    z_mesh  = handles.(handles.surfName).c_y*z;
    y_mesh  = handles.(handles.surfName).y*ones(1,n);
    
    x_mesh  = [x_mesh(end:(-1):1,:); x_mesh];
    y_mesh  = [-y_mesh(end:(-1):1,:); y_mesh];
    z_mesh  = [z_mesh(end:(-1):1,:); z_mesh];

    mesh(y_mesh,x_mesh,z_mesh);
    
    size = 30;
    color = [1 0 0];
       
    scatter3(0, handles.(handles.surfName).xca, handles.(handles.surfName).cr*t_c*0.6,size,color,'fill')
    grid on;
    
%     y_plot  = [-handles.(handles.surfName).y(end:-1:1); handles.(handles.surfName).y];
%     le_plot = [handles.(handles.surfName).le_y(end:-1:1); handles.(handles.surfName).le_y];
%     te_plot = le_plot + [handles.(handles.surfName).c_y(end:-1:1); handles.(handles.surfName).c_y];
%     
%     plot(y_plot, le_plot, 'b');
%     plot(y_plot, te_plot, 'b');
%     plot(handles.(handles.surfName).yca, handles.(handles.surfName).xca, 'rx');
    
    set(handles.cr_val,'String', num2str(handles.(handles.surfName).cr,3));
    set(handles.ct_val,'String', num2str(handles.(handles.surfName).ct,3));
    set(handles.TR_val,'String', num2str(handles.(handles.surfName).TR,3));
    set(handles.cMAC_val,'String', num2str(handles.(handles.surfName).MAC,3));
    set(handles.xca2_val,'String', num2str(handles.(handles.surfName).xca,3));
    set(handles.yca2_val,'String', num2str(handles.(handles.surfName).yca,3));
    set(handles.zca2_val,'String', num2str(handles.(handles.surfName).zca,3));
    set(handles.LAMle_val,'String', num2str(handles.(handles.surfName).LAM*180/pi,2));
    set(handles.LAMc4_val,'String', num2str(handles.(handles.surfName).LAMc4*180/pi,2));
    set(handles.LAMc2_val,'String', num2str(handles.(handles.surfName).LAMc2*180/pi,2));
    

function cac_val_Callback(hObject, eventdata, handles)
% hObject    handle to cac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cac_val as text
%        str2double(get(hObject,'String')) returns contents of cac_val as a double
handles.(handles.surfName).cm_c = str2double(get(hObject,'String'));

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
handles.(handles.surfName).y0_b2 = str2double(get(hObject,'String'));

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
handles.(handles.surfName).y1_b2 = str2double(get(hObject,'String'));

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
handles.(handles.surfName).CL0 = str2double(get(hObject,'String'));

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
handles.(handles.surfName).CLa = str2double(get(hObject,'String'));

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
handles.(handles.surfName).Cla = str2double(get(hObject,'String'));

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
handles.(handles.surfName).CM0 = str2double(get(hObject,'String'));

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
handles.(handles.surfName).eta = str2double(get(hObject,'String'));

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
handles.(handles.surfName).diedro = str2double(get(hObject,'String'))*pi/180;
handles.(handles.surfName).zca     = handles.(handles.surfName).yca*sin(handles.(handles.surfName).diedro);
set(handles.zca2_val, 'String', handles.(handles.surfName).zca);

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
handles.(handles.surfName).LAMc4 = str2double(get(hObject,'String'))*pi/180;

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
handles.(handles.surfName).LAMc2 = str2double(get(hObject,'String'))*pi/180;

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
handles.(handles.surfName).Zca = str2double(get(hObject,'String'));

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
handles.(handles.surfName).MAC = str2double(get(hObject,'String'));

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
handles.(handles.surfName).xca = str2double(get(hObject,'String'));

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
handles.(handles.surfName).yca = str2double(get(hObject,'String'));

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
handles.(handles.surfName).zca = str2double(get(hObject,'String'));

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
handles.(handles.surfName).i = str2double(get(hObject,'String'))*pi/180;

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
handles.(handles.surfName).t_c = str2double(get(hObject,'String'));

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
handles.(handles.surfName).Se = str2double(get(hObject,'String'));

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
handles.(handles.surfName).S = str2double(get(handles.S_val,'String'));
if ~isempty(handles.(handles.surfName).b)
    handles.(handles.surfName).AR = handles.(handles.surfName).b^2/handles.(handles.surfName).S;
    set(handles.AR_val,'String',num2str(handles.(handles.surfName).AR));
elseif ~isempty(handles.(handles.surfName).AR)
    handles.(handles.surfName).b = sqrt(handles.(handles.surfName).AR*handles.(handles.surfName).S);
    set(handles.b_val,'String',num2str(handles.(handles.surfName).b));
end

if ~isempty(handles.chordDistr) && ~isempty(handles.(handles.surfName).b) 
    handles = refreshSurfdata(handles);
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
handles.(handles.surfName).b = str2double(get(handles.b_val,'String'));

if ~isempty(handles.(handles.surfName).S) && isempty(handles.(handles.surfName).AR)
    handles.(handles.surfName).AR = handles.(handles.surfName).b^2/handles.(handles.surfName).S;
    set(handles.AR_val,'String',num2str(handles.(handles.surfName).AR));
elseif ~isempty(handles.(handles.surfName).AR) && isempty(handles.(handles.surfName).S)
    handles.(handles.surfName).S = handles.(handles.surfName).b^2*handles.(handles.surfName).AR;
    set(handles.S_val,'String',num2str(handles.(handles.surfName).S));
end

if ~isempty(handles.chordDistr) && ~isempty(handles.(handles.surfName).S)
    handles = refreshSurfdata(handles);    
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
handles.(handles.surfName).AR = str2double(get(hObject,'String'));

if ~isempty(handles.(handles.surfName).S) && isempty(handles.(handles.surfName).b)
    handles.(handles.surfName).b = sqrt(handles.(handles.surfName).S*handles.(handles.surfName).AR); 
    set(handles.bval,'String',num2str(handles.(handles.surfName).b));
elseif isempty(handles.(handles.surfName).S) && ~isempty(handles.(handles.surfName).b)
    handles.(handles.surfName).S = handles.(handles.surfName).b^2/handles.(handles.surfName).AR; 
    set(handles.Sval,'String',num2str(handles.(handles.surfName).S));
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



function TR_val_Callback(hObject, eventdata, handles)
% hObject    handle to Zca_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zca_val as text
%        str2double(get(hObject,'String')) returns contents of Zca_val as a double
if ~isempty(handles.(handles.surfName).cr) && isempty(handles.(handles.surfName).ct)
    handles.(handles.surfName).ct = handles.(handles.surfName).cr*handles.(handles.surfName).TR; 
    set(handles.ct_val,'String',num2str(handles.(handles.surfName).ct));
elseif ~isempty(handles.(handles.surfName).ct) && isempty(handles.(handles.surfName).cr)
    handles.(handles.surfName).cr = handles.(handles.surfName).ct/handles.(handles.surfName).TR; 
    set(handles.cr_val,'String',num2str(handles.(handles.surfName).cr));
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
handles.(handles.surfName).cr = str2double(get(handles.cr_val,'String'));

if ~isempty(handles.(handles.surfName).ct) && isempty(handles.(handles.surfName).TR)
    handles.(handles.surfName).TR = handles.(handles.surfName).ct/handles.(handles.surfName).cr;
    set(handles.TR_val,'String',num2str(handles.(handles.surfName).TR));
elseif ~isempty(handles.(handles.surfName).TR) && isempty(handles.(handles.surfName).ct)
    handles.(handles.surfName).ct = handles.(handles.surfName).cr*handles.(handles.surfName).TR;
    set(handles.ct_val,'String',num2str(handles.(handles.surfName).ct));
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
handles.(handles.surfName).ct = str2double(get(handles.ct_val,'String'));

if ~isempty(handles.(handles.surfName).cr) && isempty(handles.(handles.surfName).TR)
    handles.(handles.surfName).TR = handles.(handles.surfName).ct/handles.(handles.surfName).cr;
    set(handles.TR_val,'String',num2str(handles.(handles.surfName).TR));
elseif ~isempty(handles.(handles.surfName).TR) && isempty(handles.(handles.surfName).cr)
    handles.(handles.surfName).cr = handles.(handles.surfName).ct*handles.(handles.surfName).TR;
    set(handles.cr_val,'String',num2str(handles.(handles.surfName).cr));
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
handles.(handles.surfName).LAM = str2double(get(hObject,'String'))*pi/180;

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
handles.(handles.surfName).Xca = str2double(get(hObject,'String'));

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
if ~isequaln(handles.output,handles.(handles.surfName))
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.(handles.surfName);
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
if ~isequaln(handles.output,handles.(handles.surfName))
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.(handles.surfName);
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
