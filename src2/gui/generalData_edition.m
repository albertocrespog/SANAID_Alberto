function varargout = generalData_edition(varargin)
% GENERALDATA_EDITION MATLAB code for generalData_edition.fig
%      GENERALDATA_EDITION, by itself, creates a new GENERALDATA_EDITION or raises the existing
%      singleton*.
%
%      H = GENERALDATA_EDITION returns the handle to a new GENERALDATA_EDITION or the handle to
%      the existing singleton*.
%
%      GENERALDATA_EDITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENERALDATA_EDITION.M with the given input arguments.
%
%      GENERALDATA_EDITION('Property','Value',...) creates a new GENERALDATA_EDITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before generalData_edition_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to generalData_edition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help generalData_edition

% Last Modified by GUIDE v2.5 10-Jan-2016 21:38:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @generalData_edition_OpeningFcn, ...
                   'gui_OutputFcn',  @generalData_edition_OutputFcn, ...
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


% --- Executes just before generalData_edition is made visible.
function generalData_edition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to generalData_edition (see VARARGIN)

% Choose default command line output for generalData_edition

set(hObject, 'Name', 'General Data Edition');
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


handles.dir     = varargin{2};

handles.conv_im         = imread([handles.dir.images,'conv.png']);
handles.canard_im       = imread([handles.dir.images,'canard.png']);
handles.conv_canard_im  = imread([handles.dir.images,'conv_canard.png']);
handles.ala_vol_im      = imread([handles.dir.images,'ala_vol.png']);


handles.changes = varargin{3};

axes(handles.axes1);
axis off
if isempty(varargin{1})
    handles.output = [];
    handles.model = handles.output;
    imshow(handles.conv_im);
    set(handles.conv_button, 'Value',1);
else
    handles.output = varargin{1};
    handles.model = handles.output;
    handles = showModelData(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes generalData_edition wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = generalData_edition_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;
varargout{2} = handles.changes;
delete(handles.figure1);

function handles = showModelData(handles)

set(handles.modelName_val, 'String',handles.model.name);
if isfield(handles.model,'general')
    if isfield(handles.model.general,'mtow'), set(handles.mtow_val,'String',num2str(handles.model.general.mtow)); end
    if isfield(handles.model.general,'Sref'), set(handles.Sref_val,'String',num2str(handles.model.general.Sref)); end
    if isfield(handles.model.general,'h'), set(handles.h_val,'String',num2str(handles.model.general.h)); end
    if isfield(handles.model.general,'rhoinf'), set(handles.rho_val,'String',num2str(handles.model.general.rhoinf)); end
    if isfield(handles.model.general,'Tinf'), set(handles.T_val,'String',num2str(handles.model.general.Tinf)); end
    if isfield(handles.model.general,'pinf'), set(handles.p_val,'String',num2str(handles.model.general.pinf)); end
    if isfield(handles.model.general,'Minf'), set(handles.Minf_val,'String',num2str(handles.model.general.Minf)); end
    if isfield(handles.model.general,'Vinf'), set(handles.Vinf_val,'String',num2str(handles.model.general.Vinf)); end
    if isfield(handles.model.general,'qinf'), set(handles.qinf_val,'String',num2str(handles.model.general.qinf)); end
    if isfield(handles.model.general,'w_w0'), set(handles.w_w0_val,'String',num2str(handles.model.general.w_w0)); end
    if isfield(handles.model.general,'L'), set(handles.Lplane_val,'String',num2str(handles.model.general.L)); end
    if isfield(handles.model.general,'Vstall'), set(handles.Vstall_val,'String',num2str(handles.model.general.Vstall)); end
end
if ~isempty(handles.model.general.polar)
    switch length(handles.model.general.polar)
        case 3
            set(handles.k1_val,'String',num2str(handles.model.general.polar(2)));
            set(handles.k2_val,'String',num2str(handles.model.general.polar(3)));
            set(handles.cd0_val,'String',num2str(handles.model.general.polar(1)));
        case 2
            set(handles.cd0_val,'String',num2str(handles.model.general.polar(1)));
            set(handles.k1_val,'String',num2str(handles.model.general.polar(2)));
        case 1
            set(handles.cd0_val,'String',num2str(handles.model.general.polar(1)));
    end
            
    
end

switch handles.model.conf
    case 'convencional'
        imshow(handles.conv_im);
        set(handles.conv_button, 'Value',1);
    case 'canard'
        imshow(handles.canard_im);
        set(handles.ala_canard_button, 'Value',1);
    case 'convencional_canard'
        imshow(handles.conv_canard_im);
        set(handles.conv_canard_button, 'Value',1);
    case 'flyWing'
        imshow(handles.ala_vol_im);
        set(handles.flyWing_button, 'Value',1);
end

switch handles.model.confVert
    case 'convencional'
        set(handles.convVert_button, 'Value',1);
    case 'twin_vertical'
        set(handles.twinVert_button, 'Value',1);
    case 'no_vert'
        set(handles.noVert_button, 'Value',1);
end


function modelName_val_Callback(hObject, eventdata, handles)
% hObject    handle to modelName_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modelName_val as text
%        str2double(get(hObject,'String')) returns contents of modelName_val as a double

handles.model.name = get(handles.modelName_val, 'String');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function modelName_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelName_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mtow_val_Callback(hObject, eventdata, handles)
% hObject    handle to mtow_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mtow_val as text
%        str2double(get(hObject,'String')) returns contents of mtow_val as a double

handles.model.general.mtow = str2double(get(handles.mtow_val, 'String'));

% Update handles structure
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



function Sref_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sref_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sref_val as text
%        str2double(get(hObject,'String')) returns contents of Sref_val as a double

handles.model.general.Sref = str2double(get(handles.Sref_val, 'String'));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Sref_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sref_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in conv_button.
function conv_button_Callback(hObject, eventdata, handles)
% hObject    handle to conv_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conv_button
handles.model.conf = 'convencional';
% handles.model.ala = [];
% handles.model.horizontal = [];
% if isfield(handles.model,'canard'), handles.model = rmfield(handles.model,'canard'); end

axes(handles.axes1);
cla(handles.axes1);
axis off
imshow(handles.conv_im);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ala_canard_button.
function ala_canard_button_Callback(hObject, eventdata, handles)
% hObject    handle to ala_canard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ala_canard_button
handles.model.conf = 'ala_canard';

% handles.model.ala = [];
% handles.model.canard = [];
% if isfield(handles.model,'horizontal'), handles.model = rmfield(handles.model,'horizontal'); end

axes(handles.axes1);
cla(handles.axes1);
axis off
imshow(handles.canard_im);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in conv_canard_button.
function conv_canard_button_Callback(hObject, eventdata, handles)
% hObject    handle to conv_canard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conv_canard_button
handles.model.conf = 'convencional_canard';

% handles.model.ala = [];
% handles.model.horizontal = [];
% handles.model.canard = [];

axes(handles.axes1);
cla(handles.axes1);
axis off
imshow(handles.conv_canard_im);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in flyWing_button.
function flyWing_button_Callback(hObject, eventdata, handles)
% hObject    handle to flyWing_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flyWing_button
handles.model.conf = 'flyWing';
axes(handles.axes1);
cla(handles.axes1);
axis off
imshow(handles.ala_vol_im);


% Update handles structure
guidata(hObject, handles);


function Minf_val_Callback(hObject, eventdata, handles)
% hObject    handle to Minf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Minf_val as text
%        str2double(get(hObject,'String')) returns contents of Minf_val as a double

handles.model.general.Minf = str2double(get(handles.Minf_val, 'String'));

if isfield(handles.model.general,'ainf')
    handles.model.general.Vinf = handles.model.general.Minf*handles.model.general.ainf;
    set(handles.Vinf_val,'String',num2str(handles.model.general.Vinf));
    handles.model.general.qinf = 0.5*handles.model.general.rhoinf*(handles.model.general.Vinf)^2;
    set(handles.qinf_val,'String',num2str(handles.model.general.qinf));
end

if ~isempty(handles.model.general.h) && ~isempty(handles.model.general.Minf) && ~isempty(handles.model.propulsion.Pmax)
    cd(handles.dir.code);
    [handles.model.propulsion.Acoef,handles.model.propulsion.Bcoef,handles.model.propulsion.Ccoef] = getPropModel(handles.model.propulsion.Pmax, 1, handles.model.general.Minf, handles.model.general.h);
    cd(handles.dir.gui);
end

% Update handles structure
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


function h_val_Callback(hObject, eventdata, handles)
% hObject    handle to h_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h_val as text
%        str2double(get(hObject,'String')) returns contents of h_val as a double

handles.model.general.h = str2double(get(handles.h_val, 'String'));
cd(handles.dir.code);
[handles.model.general.Tinf, handles.model.general.rhoinf, handles.model.general.pinf, handles.model.general.ainf] = atmos_inter(handles.model.general.h);
cd(handles.dir.gui);

set(handles.rho_val,'String',num2str(handles.model.general.rhoinf));
set(handles.T_val,'String',num2str(handles.model.general.Tinf));
set(handles.p_val,'String',num2str(handles.model.general.pinf));

if isfield(handles.model.general,'Minf')
    handles.model.general.Vinf = handles.model.general.Minf*handles.model.general.ainf;
    set(handles.Vinf_val,'String',num2str(handles.model.general.Vinf));
    handles.model.general.qinf = 0.5*handles.model.general.rhoinf*(handles.model.general.Vinf)^2;
    set(handles.qinf_val,'String',num2str(handles.model.general.qinf));
end

if ~isempty(handles.model.general.h) && ~isempty(handles.model.general.Minf) && ~isempty(handles.model.propulsion.Pmax)
    cd(handles.dir.code);
    [handles.model.propulsion.Acoef,handles.model.propulsion.Bcoef,handles.model.propulsion.Ccoef] = getPropModel(handles.model.propulsion.Pmax, 1, handles.model.general.Minf, handles.model.general.h);
    cd(handles.dir.gui);
end

% Update handles structure
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



function rho_val_Callback(hObject, eventdata, handles)
% hObject    handle to rho_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rho_val as text
%        str2double(get(hObject,'String')) returns contents of rho_val as a double


% --- Executes during object creation, after setting all properties.
function rho_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rho_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_val_Callback(hObject, eventdata, handles)
% hObject    handle to T_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_val as text
%        str2double(get(hObject,'String')) returns contents of T_val as a double


% --- Executes during object creation, after setting all properties.
function T_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_val_Callback(hObject, eventdata, handles)
% hObject    handle to p_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_val as text
%        str2double(get(hObject,'String')) returns contents of p_val as a double



% --- Executes during object creation, after setting all properties.
function p_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vinf_val_Callback(hObject, eventdata, handles)
% hObject    handle to Vinf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vinf_val as text
%        str2double(get(hObject,'String')) returns contents of Vinf_val as a double



% --- Executes during object creation, after setting all properties.
function Vinf_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vinf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qinf_val_Callback(hObject, eventdata, handles)
% hObject    handle to qinf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qinf_val as text
%        str2double(get(hObject,'String')) returns contents of qinf_val as a double



% --- Executes during object creation, after setting all properties.
function qinf_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qinf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_w0_val_Callback(hObject, eventdata, handles)
% hObject    handle to w_w0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_w0_val as text
%        str2double(get(hObject,'String')) returns contents of w_w0_val as a double

handles.model.general.w_w0 = str2double(get(handles.w_w0_val, 'String'));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function w_w0_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_w0_val (see GCBO)
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
if ~isequaln(handles.output, handles.model)
    handles.changes = 1;
end
handles.output = handles.model;
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isequaln(handles.output,handles.model)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.model;
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


% --- Executes on button press in convVert_button.
function convVert_button_Callback(hObject, eventdata, handles)
% hObject    handle to convVert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convVert_button

handles.model.confVert = 'convencional';
guidata(hObject, handles);



% --- Executes on button press in twinVert_button.
function twinVert_button_Callback(hObject, eventdata, handles)
% hObject    handle to twinVert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twinVert_button

handles.model.confVert = 'twin_vertical';

guidata(hObject, handles);

% --- Executes on button press in noVert_button.
function noVert_button_Callback(hObject, eventdata, handles)
% hObject    handle to noVert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noVert_button
handles.model.confVert = 'no_vert';

guidata(hObject, handles);



function k1_val_Callback(hObject, eventdata, handles)
% hObject    handle to k1_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1_val as text
%        str2double(get(hObject,'String')) returns contents of k1_val as a double

handles.model.general.polar(2) = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function k1_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2_val_Callback(hObject, eventdata, handles)
% hObject    handle to k2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2_val as text
%        str2double(get(hObject,'String')) returns contents of k2_val as a double

handles.model.general.polar(3) = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function k2_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cd0_val_Callback(hObject, eventdata, handles)
% hObject    handle to cd0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cd0_val as text
%        str2double(get(hObject,'String')) returns contents of cd0_val as a double
handles.model.general.polar(1) = str2double(get(hObject,'String'));

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function cd0_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cd0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lplane_val_Callback(hObject, eventdata, handles)
% hObject    handle to Lplane_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lplane_val as text
%        str2double(get(hObject,'String')) returns contents of Lplane_val as a double
handles.model.general.L = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Lplane_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lplane_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Vstall_val_Callback(hObject, eventdata, handles)
% hObject    handle to Vstall_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vstall_val as text
%        str2double(get(hObject,'String')) returns contents of Vstall_val as a double
handles.model.general.Vstall = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Vstall_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vstall_val (see GCBO)
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
winopen([handles.dir.help,'GeneralData_guide.pdf']);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if ~isequaln(handles.output,handles.model)
    saveAnswer = unsavedData_dialog(handles.dir);
    switch saveAnswer
        case 'SAVE'
            handles.output = handles.model;
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
