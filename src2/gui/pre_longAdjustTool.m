function varargout = pre_longAdjustTool(varargin)
% PRE_LONGADJUSTTOOL MATLAB code for pre_longAdjustTool.fig
%      PRE_LONGADJUSTTOOL, by itself, creates a new PRE_LONGADJUSTTOOL or raises the existing
%      singleton*.
%
%      H = PRE_LONGADJUSTTOOL returns the handle to a new PRE_LONGADJUSTTOOL or the handle to
%      the existing singleton*.
%
%      PRE_LONGADJUSTTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRE_LONGADJUSTTOOL.M with the given input arguments.
%
%      PRE_LONGADJUSTTOOL('Property','Value',...) creates a new PRE_LONGADJUSTTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pre_longAdjustTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pre_longAdjustTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pre_longAdjustTool

% Last Modified by GUIDE v2.5 27-Mar-2016 03:32:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pre_longAdjustTool_OpeningFcn, ...
                   'gui_OutputFcn',  @pre_longAdjustTool_OutputFcn, ...
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


% --- Executes just before pre_longAdjustTool is made visible.
function pre_longAdjustTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pre_longAdjustTool (see VARARGIN)

% Choose default command line output for pre_longAdjustTool

%handles.output = hObject;

% Se guardan los directorios del codigo y las gui
set(hObject, 'Name', 'Preliminary Longitudinal Design');
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

handles.dir = varargin{1};

% Se adquieren los datos del modelo
handles.model   = getDefaultModel();
cd(handles.dir.code);
[handles.auxData.eq_F,handles.auxData.eq_M] = getEquilEquations();
cd(handles.dir.gui);

% Se inicia auxData que contiene datos de los Sliders
handles = getAuxData(handles);
handles = getResults(handles);
handles = refreshWindDat(handles);
% Se inicia la interfaz una vez se tienen los datos del modelo
%handles = setWindowsConf(handles);

% Se inicia la grafica
axes(handles.axes1);
[handles.surfObj, handles.centers] = initGraph(handles.model,handles.dir);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject, handles);

uiwait(handles.figure1);

function handles = getAuxData(handles)
    model   = handles.model;
    auxData = handles.auxData;
    % Center of Gravity
    auxData.Xcg.min = 0;
    auxData.Xcg.max = model.general.L;
    if isempty(model.general.Xcg) || ...
        model.general.Xcg > auxData.Xcg.max
        model.general.Xcg = 0.5*(auxData.Xcg.min+auxData.Xcg.max);
    end
    auxData.Xcg.val     = model.general.Xcg; 
    
    % Wing
    auxData.wing.Xmin   = 0;
    auxData.wing.Xmax   = model.general.L;
    if      isempty(model.ala.Xca) || ...
            model.ala.Xca < auxData.wing.Xmin || ...
            model.ala.Xca > auxData.wing.Xmax
        model.ala.Xca = 0.5*(auxData.wing.Xmin+auxData.wing.Xmax);
    end
    auxData.wing.X      = model.ala.Xca;
    
    % Horizontal
    if strcmp(model.conf,'convencional') || strcmp(model.conf,'convencional_canard')
        auxData.hor.Xmin   = model.ala.Xca;
        auxData.hor.Xmax   = model.general.L;
        if  isempty(model.horizontal.Xca) || ... 
            (model.horizontal.Xca < auxData.hor.Xmin) || ...  
            (model.horizontal.Xca > auxData.hor.Xmax)
            model.horizontal.Xca = 0.5*(auxData.hor.Xmin+auxData.hor.Xmax);
        end
        auxData.hor.X      = model.horizontal.Xca;
    end
    
    % Canard
    if strcmp(model.conf,'canard') || strcmp(model.conf,'convencional_canard')
        auxData.can.Xmin   = 0;
        auxData.can.Xmax   = model.ala.Xca;
        if  isempty(model.canard.Xca) || ...
            (model.canard.Xca < auxData.can.Xmin) || ...
            (model.canard.Xca > auxData.can.Xmax)
            model.canard.Xca = 0.5*(auxData.can.Xmin+auxData.can.Xmax);
        end
        auxData.can.X      = model.canard.Xca;
    end
    % Incidences
    if ~isfield(auxData,'iFix')
        auxData.iFix = 'w';
        auxData.iFix_val = model.ala.i;
    end
    
    if isempty(model.ala.i)
        model.ala.i = 0;
    end
    
    handles.auxData = auxData;
    handles.model = model;

    
function handles = refreshWindDat(handles)
    auxData = handles.auxData;
    model   = handles.model;
    
    set(handles.L_val,'String',model.general.L);
    set(handles.Sref_val,'String',model.general.Sref);
    set(handles.c_val,'String',model.ala.MAC);
    
    % Center of Gravity Windows Data
    set(handles.Xcg_min, 'String', num2str(auxData.Xcg.min));
    set(handles.Xcg_max, 'String', num2str(auxData.Xcg.max));
    set(handles.Xcg_val, 'String', num2str(auxData.Xcg.val));
    set(handles.Xcg_slide, 'Value', auxData.Xcg.val, 'Min', auxData.Xcg.min, 'Max', auxData.Xcg.max);    
    
    % Wing Windows Data
    set(handles.Sw_val,'String',model.ala.S);
    set(handles.CLaw_val,'String',model.ala.CLa);
    set(handles.Xw_slide,'Visible','on','Value',auxData.wing.X,'Max',auxData.wing.Xmax,'Min',auxData.wing.Xmin);
    set(handles.Xw_val, 'Visible', 'on','String',num2str(auxData.wing.X,3));
    set(handles.Xw_min_text,'Visible','on','String', num2str(auxData.wing.Xmin));
    set(handles.Xw_max_text,'Visible','on','String', num2str(auxData.wing.Xmax));
    
    set(handles.iw_text,'Visible','on');
    set(handles.iw_val,'Visible','on','String',num2str(handles.model.ala.i*180/pi));
    
    % Stabilizing Surfaces Windows Data
    if strcmp(model.conf,'convencional_canard')
        set(handles.convCanard_button,'Value',1);
        
        set(handles.Sh_text,'Visible','on')
        set(handles.Sh_val,'Visible','on','String',model.horizontal.S);
        set(handles.CLah_text,'Visible','on');
        set(handles.CLah_val,'Visible','on','String',model.horizontal.CLa);
        set(handles.Xh_slide,'Visible','on','Value',auxData.hor.X,'Max',auxData.hor.Xmax,'Min',auxData.hor.Xmin);
        set(handles.Xh_text, 'Visible','on');
        set(handles.Xh_val, 'Visible', 'on','String',num2str(auxData.hor.X,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(auxData.hor.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(auxData.hor.Xmax));
        
        
        set(handles.Sc_text,'Visible','on')
        set(handles.Sc_val,'Visible','on','String',model.canard.S);
        set(handles.CLac_text,'Visible','on');
        set(handles.CLac_val,'Visible','on','String',model.canard.CLa);
        set(handles.Xc_slide,'Visible','on','Value',auxData.can.X,'Max',auxData.can.Xmax,'Min',auxData.can.Xmin);
        set(handles.Xc_text, 'Visible','on');
        set(handles.Xc_val, 'Visible', 'on','String',num2str(auxData.can.X,3));
        set(handles.Xc_min_text,'Visible','on','String', num2str(auxData.can.Xmin));
        set(handles.Xc_max_text,'Visible','on','String', num2str(auxData.can.Xmax));
        
        set(handles.i_buttonGroup, 'Visible', 'on')
        set(handles.ih_text,'Visible','on');
        set(handles.ih_val,'Visible','on','String',num2str(handles.model.horizontal.i*180/pi));
        set(handles.ic_text,'Visible','on');
        set(handles.ic_val,'Visible','on','String',num2str(handles.model.canard.i*180/pi));
        %set(handles.iw_button,'Value', 1);
        set(handles.fixInc_val,'String', num2str(auxData.iFix_val*180/pi));
        
    
    elseif strcmp(model.conf,'convencional')
        set(handles.conv_button,'Value',1);
        
        set(handles.Sh_text,'Visible','on');
        set(handles.Sh_val,'Visible','on','String',model.horizontal.S);
        set(handles.CLah_text,'Visible','on');
        set(handles.CLah_val,'Visible','on','String',model.horizontal.CLa);
        set(handles.Xh_slide,'Visible','on','Value',auxData.hor.X,'Max',auxData.hor.Xmax,'Min',auxData.hor.Xmin);
        set(handles.Xh_text, 'Visible','on');
        set(handles.Xh_val, 'Visible', 'on','String',num2str(auxData.hor.X,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(auxData.hor.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(auxData.hor.Xmax));
        set(handles.ih_text,'Visible','on');
        set(handles.ih_val,'Visible','on','String',num2str(handles.model.horizontal.i*180/pi));
        
        set(handles.Sc_text,'Visible','off');
        set(handles.Sc_val,'Visible','off');
        set(handles.CLac_text,'Visible','off');
        set(handles.CLac_val,'Visible','off');
        set(handles.Xc_slide,'Visible','off');
        set(handles.Xc_text, 'Visible','off');
        set(handles.Xc_val, 'Visible','off');
        set(handles.Xc_min_text,'Visible','off');
        set(handles.Xc_max_text,'Visible','off');
        set(handles.i_buttonGroup, 'Visible', 'off');
        set(handles.ic_text,'Visible','off');
        set(handles.ic_val,'Visible','off');
        
        set(handles.i_buttonGroup, 'Visible', 'off');
        
    elseif strcmp(model.conf,'canard')
        set(handles.canard_button,'Value',1);
        
        set(handles.Sc_text,'Visible','on');
        set(handles.Sc_val,'Visible','on','String',model.canard.S);
        set(handles.CLac_text,'Visible','on');
        set(handles.CLac_val,'Visible','on','String',model.canard.CLa);
        set(handles.Xc_slide,'Visible','on','Value',auxData.can.X,'Max',auxData.can.Xmax,'Min',auxData.can.Xmin);
        set(handles.Xc_text, 'Visible','on');
        set(handles.Xc_val, 'Visible', 'on','String',num2str(auxData.can.X,3));
        set(handles.Xc_min_text,'Visible','on','String', num2str(auxData.can.Xmin));
        set(handles.Xc_max_text,'Visible','on','String', num2str(auxData.can.Xmax));
        set(handles.i_buttonGroup, 'Visible', 'off');
        set(handles.ic_text,'Visible','on');
        set(handles.ic_val,'Visible','on','String',num2str(handles.model.canard.i*180/pi));
            
        set(handles.Sh_text,'Visible','off');
        set(handles.Sh_val,'Visible','off');
        set(handles.CLah_text,'Visible','off');
        set(handles.CLah_val,'Visible','off');
        set(handles.Xh_slide,'Visible','off');
        set(handles.Xh_text, 'Visible','off');
        set(handles.Xh_val, 'Visible','off');
        set(handles.Xh_min_text,'Visible','off');
        set(handles.Xh_max_text,'Visible','off');
        set(handles.ih_text,'Visible','off');
        set(handles.ih_val,'Visible','off');
        
        set(handles.i_buttonGroup, 'Visible', 'off');
    
    elseif strcmp(model.conf,'flyWing')
        set(handles.flyWing_button,'Value',1);
        
        set(handles.Sc_text,'Visible','off');
        set(handles.Sc_val,'Visible','off');
        set(handles.CLac_text,'Visible','off');
        set(handles.CLac_val,'Visible','off');
        set(handles.Xc_slide,'Visible','off');
        set(handles.Xc_text, 'Visible','off');
        set(handles.Xc_val, 'Visible','off');
        set(handles.Xc_min_text,'Visible','off');
        set(handles.Xc_max_text,'Visible','off');
        set(handles.ic_text,'Visible','off');
        set(handles.ic_val,'Visible','off');
        
        set(handles.Sh_text,'Visible','off');
        set(handles.Sh_val,'Visible','off');
        set(handles.CLah_text,'Visible','off');
        set(handles.CLah_val,'Visible','off');
        set(handles.Xh_slide,'Visible','off');
        set(handles.Xh_text, 'Visible','off');
        set(handles.Xh_val, 'Visible','off');
        set(handles.Xh_min_text,'Visible','off');
        set(handles.Xh_max_text,'Visible','off');
        set(handles.ih_text,'Visible','off');
        set(handles.ih_val,'Visible','off');
        
        set(handles.i_buttonGroup, 'Visible', 'off');
    end


    
function handles = getResults(handles)
    cd(handles.dir.code);
    handles.model.general.Xna = pre_aeroCenter(handles.model);
    cd(handles.dir.gui);

    handles.auxData.ME = 1/handles.model.ala.MAC*(handles.model.general.Xna - handles.model.general.Xcg);
    set(handles.Xna_val,'String', num2str(handles.model.general.Xna,3));
    set(handles.ME_val,'String', num2str(handles.auxData.ME,3));


    
% --- Outputs from this function are returned to the command line.
function varargout = pre_longAdjustTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
delete(handles.figure1);


% --- Executes on slider movement.
function Xh_slide_Callback(hObject, eventdata, handles)
% hObject    handle to Xh_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Actualizacion de los datos internos de la aplicacion

handles.model.horizontal.Xca = get(handles.Xh_slide,'Value');
set(handles.Xh_val,'String', num2str(handles.model.horizontal.Xca,3));

% Ejecucion del calculo del punto neutro y almacenamiento de dicha variable
% en la estructura correspondiente. Actualizacion de los campos de texto de
% la aplicacion.
handles = getAuxData(handles); 
handles = getResults(handles);

% Actualizacion de la grafica
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Xh_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xh_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Xw_slide_Callback(hObject, eventdata, handles)
% hObject    handle to Xw_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.model.ala.Xca = get(handles.Xw_slide,'Value');

handles = getAuxData(handles);
handles = refreshWindDat(handles);

handles = getResults(handles);

refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Xw_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xw_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Xc_slide_Callback(hObject, eventdata, handles)
% hObject    handle to Xc_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.model.canard.Xca = get(handles.Xc_slide,'Value');
set(handles.Xc_val,'String', num2str(handles.model.canard.Xca,3));

handles = getResults(handles);

% Actualizacion de la grafica
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Xc_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xc_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Xw_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xw_val as text
%        str2double(get(hObject,'String')) returns contents of Xw_val as a double
%str2double(get(hObject,'String'))

% Actualizacion de la informacion interna de la aplicacion
handles.model.ala.Xca =  get(handles.Xw_val,'Value');
handles = getAuxData(handles);
handles = refreshWindDat(handles);

handles = getResults(handles);

% Actualizacion de las graficas
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Xw_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Xc_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xc_val as text
%        str2double(get(hObject,'String')) returns contents of Xc_val as a double
handles.auxData.can.X = get(handles.Xc_val,'Value'); 
if  (handles.auxData.can.X <= handles.auxData.can.Xmax) && ...
    (handles.auxData.can.X >= handles.auxData.can.Xmin)
    handles.model.canard.Xca =  handles.auxData.can.X;
    set(handles.Xc_slider,'Value',handles.auxData.can.X);
    handles = getResults(handles);
    refreshGraph(handles.surfObj, handles.centers, handles.model);
else
    handles.auxData.can.X = handles.model.canard.Xca;
    set(handles.Xc_val,'String',num2str(handles.auxData.can.X,3));
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Xc_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Xh_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xh_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xh_val as text
%        str2double(get(hObject,'String')) returns contents of Xh_val as a double

handles.auxData.hor.X = get(handles.Xh_val,'Value'); 
if  (handles.auxData.hor.X <= handles.auxData.hor.Xmax) && ...
    (handles.auxData.hor.X >= handles.auxData.hor.Xmin)
    handles.model.horizontal.Xca =  handles.auxData.hor.X;
    set(handles.Xh_slider,'Value',handles.auxData.hor.X);
    handles = getResults(handles);
    refreshGraph(handles.surfObj, handles.centers, handles.model);
else
    handles.auxData.hor.X = handles.model.horizontal.Xca;
    set(handles.Xh_val,'String',num2str(handles.auxData.hor.X,3));
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Xh_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xh_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CLaw_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLaw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLaw_val as text
%        str2double(get(hObject,'String')) returns contents of CLaw_val as a double

handles.model.ala.CLa = str2double(get(hObject,'String'));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);
    
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CLaw_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLaw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CLah_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLah_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLah_val as text
%        str2double(get(hObject,'String')) returns contents of CLah_val as a double

handles.model.horizontal.CLa = str2double(get(hObject,'String'));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CLah_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLah_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CLac_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLac_val as text
%        str2double(get(hObject,'String')) returns contents of CLac_val as a double

handles.model.canard.CLa = str2double(get(hObject,'String'));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CLac_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLac_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Sw_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sw_val as text
%        str2double(get(hObject,'String')) returns contents of Sw_val as a double

handles.model.ala.S = str2double(get(hObject,'String'));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Sw_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Sh_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sh_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sh_val as text
%        str2double(get(hObject,'String')) returns contents of Sh_val as a double

handles.model.horizontal.S = str2double(get(hObject,'String'));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Sh_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sh_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Sc_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sc_val as text
%        str2double(get(hObject,'String')) returns contents of Sc_val as a double
handles.model.canard.S = str2double(get(hObject,'String'));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Sc_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function L_val_Callback(hObject, eventdata, handles)
% hObject    handle to L_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_val as text
%        str2double(get(hObject,'String')) returns contents of L_val as a double
if str2double(get(hObject,'String')) > 0
    handles.model.general.L = str2double(get(hObject,'String'));
    handles = getAuxData(handles);
    handles = refreshWindDat(handles);
    handles = getResults(handles);
    [handles.surfObj, handles.centers] = initGraph(handles.model,handles.dir);
    refreshGraph(handles.surfObj, handles.centers, handles.model);
else
    set(hObject,'String',num2Str(handles.model.general.L,3));
end

guidata(hObject,handles)

function Sref_val_Callback(hObject, eventdata, handles)
% hObject    handle to Sref_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sref_val as text
%        str2double(get(hObject,'String')) returns contents of Sref_val as a double
if str2double(get(hObject,'String')) > 0
    handles.model.general.Sref = str2double(get(hObject,'String'));
    handles = getResults(handles);
    refreshGraph(handles.surfObj, handles.centers, handles.model);
else
    set(hObject,'String',num2Str(handles.model.general.Sref,3));
end

guidata(hObject,handles)

function Xcg_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xcg_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xcg_val as text
%        str2double(get(hObject,'String')) returns contents of Xcg_val as a double

handles.auxData.Xcg.val = str2double(get(hObject,'String'));
if  (handles.auxData.Xcg.val > handles.auxData.Xcg.min) && ...
    (handles.auxData.Xcg.val < handles.auxData.Xcg.max) 
    handles.model.general.Xcg = handles.auxData.Xcg.val;
    set(handles.Xcg_slide, 'Value', handles.auxData.Xcg.val);
    handles = getResults(handles);
    refreshGraph(handles.surfObj, handles.centers, handles.model);
else
    handles.auxData.Xcg.val = handles.model.general.Xcg;
    set(hObject,'String',num2str(handles.auxData.Xcg.val));
end

guidata(hObject,handles);


% --- Executes on slider movement.
function Xcg_slide_Callback(hObject, eventdata, handles)
% hObject    handle to Xcg_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.auxData.Xcg.val = get(hObject,'Value');
handles.model.general.Xcg = handles.auxData.Xcg.val;
set(handles.Xcg_val, 'String', num2str(handles.auxData.Xcg.val));
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Xcg_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xcg_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in iw_button.
function iw_button_Callback(hObject, eventdata,handles)
% hObject    handle to iw_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iw_button
handles.auxData.iFix = 'w';
handles.model.ala.i = handles.auxData.iFix_val;
handles.model.canard.i = 0;
handles.model.horizontal.i = 0;
set(handles.iw_val,'String',num2str(handles.model.ala.i*180/pi));
set(handles.ih_val,'String','');
set(handles.ic_val,'String','');
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject,handles);


% --- Executes on button press in ih_button.
function ih_button_Callback(hObject, eventdata,handles)
% hObject    handle to ih_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ih_button
handles.auxData.iFix = 'h';
handles.model.horizontal.i = handles.auxData.iFix_val;
handles.model.canard.i = 0;
handles.model.ala.i = 0;
set(handles.iw_val,'String','');
set(handles.ih_val,'String',num2str(handles.model.horizontal.i*180/pi));
set(handles.ic_val,'String','');
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject,handles);

% --- Executes on button press in ic_button.
function ic_button_Callback(hObject,eventdata,handles)
% hObject    handle to ic_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ic_button
handles.auxData.iFix = 'c';
handles.model.canard.i = handles.auxData.iFix_val;
handles.model.horizontal.i = [];
handles.model.ala.i = [];
set(handles.iw_val,'String',num2str(handles.model.ala.i*180/pi));
set(handles.ih_val,'String',num2str(handles.model.horizontal.i*180/pi));
set(handles.ic_val,'String',num2str(handles.model.canard.i*180/pi));
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject,handles);

function fixInc_val_Callback(hObject, eventdata, handles)
% hObject    handle to fixInc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixInc_val as text
%        str2double(get(hObject,'String')) returns contents of fixInc_val as a double
handles.auxData.iFix_val = str2double(get(hObject,'String'))*pi/180;
set(hObject,'BackgroundColor',[255/255, 255/255, 255/255]);
switch handles.auxData.iFix
    case 'w'
        handles.model.ala.i = handles.auxData.iFix_val;
        handles.model.horizontal.i = [];
        handles.model.canard.i = [];
    case 'c'
        handles.model.canard.i = handles.auxData.iFix_val;
        handles.model.horizontal.i = [];
        handles.model.ala.i = [];
    case 'h'
        handles.model.horizontal.i = handles.auxData.iFix_val;
        handles.model.canard.i = [];
        handles.model.ala.i = [];
end

set(handles.iw_val,'String',num2str(handles.model.ala.i*180/pi));
set(handles.ih_val,'String',num2str(handles.model.horizontal.i*180/pi));
set(handles.ic_val,'String',num2str(handles.model.canard.i*180/pi));

refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function fixInc_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixInc_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcula_button.
function calcula_button_Callback(hObject, eventdata, handles)
% hObject    handle to calcula_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Resolucion del sistema de ecuaciones y obtencion de las incidencias
if strcmp(handles.model.conf,'convencional_canard')
    if (isfield(handles.auxData,'iFix_val') && ~isempty(handles.auxData.iFix_val))
        handles.model = checkModel(handles.model);
        cd(handles.dir.code);
        handles = pre_solveEquil(handles);
        cd(handles.dir.gui);

        handles = getAuxData(handles);
        handles = getResults(handles);
        handles = refreshWindDat(handles);

        refreshGraph(handles.surfObj, handles.centers, handles.model);
    else
        set(handles.fixInc_val,'BackgroundColor',[255/255, 153/255, 153/255]);
    end
else
    handles.model = checkModel(handles.model);
        cd(handles.dir.code);
        handles = pre_solveEquil(handles);
        cd(handles.dir.gui);

        handles = getAuxData(handles);
        handles = getResults(handles);
        handles = refreshWindDat(handles);

        refreshGraph(handles.surfObj, handles.centers, handles.model);
end

guidata(hObject,handles)


% --- Executes on button press in import_button.
function import_button_Callback(hObject, eventdata, handles)
% hObject    handle to import_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[handles.FileName, handles.PathName, ~] = uigetfile('*.mat','Seleccione modelo','../models');
[handles.FileName, handles.PathName, ~] = uigetfile('*.mat','Seleccione modelo',handles.dir.model);
if ischar(handles.FileName) && ischar(handles.PathName)
    handles.model2import = fullfile(handles.PathName, handles.FileName);
    loadedData = load(handles.model2import);
    modelName = char(fieldnames(loadedData));
    handles.model = loadedData.(modelName);
    [check,log] = checkImportModel(handles.model);
    if check
        handles.model = checkModel(handles.model);

        handles = getAuxData(handles);
        handles = getResults(handles);
        handles = refreshWindDat(handles);

        refreshGraph(handles.surfObj, handles.centers, handles.model);
    else
        mes = 'Next fields need to be specified to perform Preliminary Longitudinal Analysis:';
        insufficientData_dialog(handles.dir,mes,log);
    end
end
guidata(hObject,handles)


% --- Executes on button press in conv_button.
function conv_button_Callback(hObject, eventdata, handles)
% hObject    handle to conv_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conv_button

handles.model.conf = 'convencional';

handles = getAuxData(handles);
handles = refreshWindDat(handles);
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)


% --- Executes on button press in canard_button.
function canard_button_Callback(hObject, eventdata, handles)
% hObject    handle to canard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of canard_button
handles.model.conf = 'canard';

handles = getAuxData(handles);
handles = refreshWindDat(handles);
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)


% --- Executes on button press in convCanard_button.
function convCanard_button_Callback(hObject, eventdata, handles)
% hObject    handle to convCanard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convCanard_button
handles.model.conf = 'convencional_canard';

handles = getAuxData(handles);
handles = refreshWindDat(handles);
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes on button press in flyWing_button.
function flyWing_button_Callback(hObject, eventdata, handles)
% hObject    handle to flyWing_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flyWing_button
handles.model.conf = 'flyWing';

handles = getAuxData(handles);
handles = refreshWindDat(handles);
handles = getResults(handles);
refreshGraph(handles.surfObj, handles.centers, handles.model);

guidata(hObject,handles)

% --- Executes on button press in extraData_button.
function extraData_button_Callback(hObject, eventdata, handles)
% hObject    handle to extraData_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.model = longTool_ediMenut(handles.model,handles.dir);
% Se inicia auxData que contiene datos de los Sliders
handles.model = checkModel(handles.model);
handles = getAuxData(handles);
handles = getResults(handles);
handles = refreshWindDat(handles);

% Se inicia la interfaz una vez se tienen los datos del modelo
%handles = setWindowsConf(handles);

% Se inicia la grafica
[handles.surfObj, handles.centers] = initGraph(handles.model,handles.dir);
refreshGraph(handles.surfObj, handles.centers, handles.model);
guidata(hObject,handles)


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
guidata(hObject,handles)


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
modelo = handles.model;
[modelName,saveBol] = preLongAdj_selectName(modelo.name);
if saveBol
    fileName = [handles.dir.model,modelName,'.mat'];
    save(fileName, 'modelo');
end
guidata(hObject,handles)

function [check, log] = checkImportModel(model)
log = '';
check = 1;

if isempty(model.general.Sref)
    log = [log, 'Sref '];
    check = 0;
end

if isempty(model.general.mtow)
    log = [log,'MTOW '];
    check = 0;
end

if isempty(model.general.h)
    log = [log,'h '];
    check = 0;
end

if isempty(model.general.Minf)
    log = [log,'Minf '];
    check = 0;
end

if isempty(model.general.L)
    log = [log, 'L '];
    check = 0;
end

if isempty(model.ala.MAC)
    log = [log,'MAC_w '];
    check = 0;
end

if isempty(model.ala.eta)
    log = [log,'eta_w '];
    check = 0;
end

if isempty(model.ala.TR)
    log = [log,'TR_w '];
    check = 0;
end

if isempty(model.ala.AR)
    log = [log,'AR_w, '];
    check = 0;
end

if isempty(model.ala.b)
    log = [log,'b_w '];
    check = 0;
end

if isempty(model.ala.LAMc4)
    log = [log,'LAMc4_w '];
    check = 0;
end



function model = checkModel(model)

% GENERAL DATA
if isempty(model.general.rhoinf)
    model.general.Xcg = model.general.L/2;
end

% WING DATA
if isempty(model.ala.S)
    model.ala.S = 1;
end

if isempty(model.ala.CLa)
    model.ala.CLa = 5;
end

if isempty(model.ala.CL0)
    model.ala.CL0 = 0;
end

if isempty(model.ala.CM0)
    model.ala.CM0 = 0;
end

% CANARD DATA
if isempty(model.canard.S)
    model.canard.S = 1;
end

if isempty(model.canard.MAC)
    model.canard.MAC = model.ala.MAC*0.2;
end

if isempty(model.canard.CLa)
    model.canard.CLa = 5;
end

if isempty(model.canard.CL0)
    model.canard.CL0 = 0;
end

if isempty(model.canard.CM0)
    model.canard.CM0 = 0;
end

% HORIZONTAL DATA
if isempty(model.horizontal.S)
    model.horizontal.S = 1;
end

if isempty(model.horizontal.MAC)
    model.horizontal.MAC = model.ala.MAC*0.2;
end

if isempty(model.horizontal.CLa)
    model.horizontal.CLa = 5;
end

if isempty(model.horizontal.CL0)
    model.horizontal.CL0 = 0;
end

if isempty(model.horizontal.CM0)
    model.horizontal.CM0 = 0;
end
