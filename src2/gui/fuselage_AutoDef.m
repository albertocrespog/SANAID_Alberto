function varargout = fuselage_AutoDef(varargin)
%FUSELAGE_AUTODEF M-file for fuselage_AutoDef.fig
%      FUSELAGE_AUTODEF, by itself, creates a new FUSELAGE_AUTODEF or raises the existing
%      singleton*.
%
%      H = FUSELAGE_AUTODEF returns the handle to a new FUSELAGE_AUTODEF or the handle to
%      the existing singleton*.
%
%      FUSELAGE_AUTODEF('Property','Value',...) creates a new FUSELAGE_AUTODEF using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to fuselage_AutoDef_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      FUSELAGE_AUTODEF('CALLBACK') and FUSELAGE_AUTODEF('CALLBACK',hObject,...) call the
%      local function named CALLBACK in FUSELAGE_AUTODEF.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fuselage_AutoDef

% Last Modified by GUIDE v2.5 22-Mar-2016 20:39:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fuselage_AutoDef_OpeningFcn, ...
                   'gui_OutputFcn',  @fuselage_AutoDef_OutputFcn, ...
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


% --- Executes just before fuselage_AutoDef is made visible.
function fuselage_AutoDef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for fuselage_AutoDef
set(hObject, 'Name', 'Fuselage Edition');
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

if isempty(varargin)
    handles.fuselaje.l          = [];
    handles.fuselaje.D          = [];
    handles.fuselaje.W          = [];
    handles.fuselaje.Sfront     = [];
    handles.fuselaje.Sside      = [];
    handles.fuselaje.Stop       = [];
    handles.fuselaje.vol        = [];
    handles.fuselaje.S0         = [];
    handles.fuselaje.CLa        = [];
    handles.output              = handles.fuselaje;
else
    handles.modelo      = varargin{1};
    handles.fuselaje    = handles.modelo.fuselaje;
    handles.output      = handles.fuselaje;
   
    set(handles.lfus_val, 'String', num2str(handles.fuselaje.l));
    set(handles.Dmax_val, 'String', num2str(handles.fuselaje.D));
    set(handles.Wmax_val, 'String', num2str(handles.fuselaje.W));
    set(handles.Sfront_val, 'String', num2str(handles.fuselaje.Sfront));
    set(handles.Sside_val, 'String', num2str(handles.fuselaje.Sside));
    set(handles.Stop_val, 'String', num2str(handles.fuselaje.Stop));
    set(handles.vol_val, 'String', num2str(handles.fuselaje.vol));
    set(handles.CLa_val, 'String', num2str(handles.fuselaje.CLa));
    
    if isempty(handles.fuselaje.meshData)
        rotate3d(handles.axes1);
    else    
        refresh_fusGraph(handles.fuselaje.meshData, handles.axes1);
        rotate3d(handles.axes1);
    end
%     if ~isempty(handles.fuselaje.meshData)
%         meshData = handles.fuselaje.meshData;
%         hold on;
%         C(:,:,1) = 0.2*ones(size(meshData{1}));
%         C(:,:,2) = 0.3*ones(size(meshData{1}));
%         C(:,:,3) = 0.7*ones(size(meshData{1}));
%         mesh(meshData{1},meshData{2},meshData{3},C);
%         axis equal;
%         rotate3d(handles.axes1);
%         grid on;
%     end
end

handles.dir     = varargin{2};
handles.changes = varargin{3};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fuselage_AutoDef wait for user response (see UIRESUME)
uiwait(handles.figure1);

function refresh_fusGraph(meshData,fus_axes)

C(:,:,1) = 0.2*ones(size(meshData{1}));
C(:,:,2) = 0.3*ones(size(meshData{1}));
C(:,:,3) = 0.7*ones(size(meshData{1}));
axes(fus_axes);
cla(fus_axes);
hold on;
mesh(meshData{1},meshData{2},meshData{3},C);
axis equal;
grid on;



% --- Outputs from this function are returned to the command line.
function varargout = fuselage_AutoDef_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.changes;
delete(handles.figure1);



% --- Executes on button press in openFuselageFile_button.
function openFuselageFile_button_Callback(hObject, eventdata, handles)
% hObject    handle to openFuselageFile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentDir = pwd;
cd(handles.dir.main);
[handles.FileName, handles.PathName, ~] = uigetfile('*.txt','Select fuselage .txt file');
cd(currentDir);
if ~handles.PathName
    handles.fuselageFilePath = [handles.PathName, handles.FileName];
    set(handles.filePath,'String',handles.fuselageFilePath);
end
guidata(hObject, handles);



% --- Executes on button press in loadFuselageFile_button.
function loadFuselageFile_button_Callback(hObject, eventdata, handles)
% hObject    handle to loadFuselageFile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fuselageFilePath = [handles.PathName, handles.FileName];
% Chech if there is a preexistent fuselage in the model
emptyFus = 0;
if isempty(handles.fuselaje.l)
    emptyFus = 1;
elseif isnan(handles.fuselaje.l)
    emptyFus = 1;
else
    lnew = handles.fuselaje.l;
end
cd(handles.dir.code);
[handles.fuselaje] = getAdimFus(handles.fuselageFilePath, handles.nPoint, handles.nSection, handles.skip);
if emptyFus
    lnew = handles.fuselaje.l;
end
handles.fuselaje = fusRefresh(handles.fuselaje,lnew,handles.modelo.general.Sref);
cd(handles.dir.gui);

refresh_fusGraph(handles.fuselaje.meshData, handles.axes1);

% handles.fuselaje.meshData = meshData;
% hold on;
% C(:,:,1) = 0.2*ones(size(meshData{1}));
% C(:,:,2) = 0.3*ones(size(meshData{1}));
% C(:,:,3) = 0.7*ones(size(meshData{1}));
% mesh(meshData{1},meshData{2},meshData{3},C);
% axis equal;
% rotate3d(handles.axes1);
% grid on;

set(handles.lfus_val,'String', num2str(handles.fuselaje.l));
set(handles.Dmax_val,'String', num2str(handles.fuselaje.D));
set(handles.Wmax_val,'String', num2str(handles.fuselaje.W));
set(handles.Sfront_val,'String', num2str(handles.fuselaje.Sfront));
set(handles.Sside_val,'String', num2str(handles.fuselaje.Sside));
set(handles.Stop_val, 'String', num2str(handles.fuselaje.Stop));
set(handles.vol_val,'String', num2str(handles.fuselaje.vol));
set(handles.CLa_val,'String', num2str(handles.fuselaje.CLa));

guidata(hObject, handles);



function filePath_Callback(hObject, eventdata, handles)
% hObject    handle to filePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filePath as text
%        str2double(get(hObject,'String')) returns contents of filePath as a double
handles.fuselageFilePath = get(handles.filePath,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function filePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filePath (see GCBO)
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
if ~isequaln(handles.output,handles.fuselaje)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.fuselaje;
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
if ~isequaln(handles.output, handles.fuselaje)
    handles.changes = 1;
end

handles.output = handles.fuselaje;
guidata(hObject, handles);



function lfus_val_Callback(hObject, eventdata, handles)
% hObject    handle to lfus_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lfus_val as text
%        str2double(get(hObject,'String')) returns contents of lfus_val as a double
lnew = str2double(get(hObject,'String'));

handles.fuselaje = fusRefresh(handles.fuselaje, lnew, handles.modelo.general.Sref);

refresh_fusGraph(handles.fuselaje.meshData, handles.axes1);


set(handles.lfus_val,'String', num2str(handles.fuselaje.l));
set(handles.Dmax_val,'String', num2str(handles.fuselaje.D));
set(handles.Wmax_val,'String', num2str(handles.fuselaje.W));
set(handles.Sfront_val,'String', num2str(handles.fuselaje.Sfront));
set(handles.Sside_val,'String', num2str(handles.fuselaje.Sside));
set(handles.Stop_val, 'String', num2str(handles.fuselaje.Stop));
set(handles.vol_val,'String', num2str(handles.fuselaje.vol));
set(handles.CLa_val,'String', num2str(handles.fuselaje.CLa));


guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function lfus_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lfus_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dmax_val_Callback(hObject, eventdata, handles)
% hObject    handle to Dmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dmax_val as text
%        str2double(get(hObject,'String')) returns contents of Dmax_val as a double
handles.fuselaje.D = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Dmax_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wmax_val_Callback(hObject, eventdata, handles)
% hObject    handle to Wmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Wmax_val as text
%        str2double(get(hObject,'String')) returns contents of Wmax_val as a double
handles.fuselaje.W = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Wmax_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wmax_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sfront_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sfront_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sfront_val as text
%        str2double(get(hObject,'String')) returns contents of Sfront_val as a double

handles.fuselaje.Sfront = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Sfront_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sfront_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sside_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sside_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sside_val as text
%        str2double(get(hObject,'String')) returns contents of Sside_val as a double
handles.fuselaje.Sside = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Sside_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sside_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Stop_val_Callback(hObject, eventdata, handles)
% hObject    handle to Stop_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Stop_val as text
%        str2double(get(hObject,'String')) returns contents of Stop_val as a double
handles.fuselaje.Stop = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Stop_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stop_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vol_val_Callback(hObject, eventdata, handles)
% hObject    handle to vol_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vol_val as text
%        str2double(get(hObject,'String')) returns contents of vol_val as a double
handles.fuselaje.vol = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function vol_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vol_val (see GCBO)
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
handles.fuselaje.CLa = str2double(get(hObject,'String'));

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

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isequaln(handles.output,handles.fuselaje)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.fuselaje;
            handles.changes = 1;
            uiresume(handles.figure1);
        case 'DISCARD'
            uiresume(handles.figure1);
        case 'CANCEL'
            % Do nothing
    end
else
    uiresume(handles.figure1);
end

% Hint: delete(hObject) closes the figure
% guidata(hObject, handles);
% uiresume(handles.figure1);



function nSection_val_Callback(hObject, eventdata, handles)
% hObject    handle to nSection_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nSection_val as text
%        str2double(get(hObject,'String')) returns contents of nSection_val as a double
handles.nSection = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nSection_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nSection_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nPoint_val_Callback(hObject, eventdata, handles)
% hObject    handle to nPoint_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nPoint_val as text
%        str2double(get(hObject,'String')) returns contents of nPoint_val as a double
handles.nPoint = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nPoint_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nPoint_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function skip_val_Callback(hObject, eventdata, handles)
% hObject    handle to skip_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of skip_val as text
%        str2double(get(hObject,'String')) returns contents of skip_val as a double
handles.skip = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function skip_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to skip_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen([handles.dir.help,'fuselage_guide.pdf']);
