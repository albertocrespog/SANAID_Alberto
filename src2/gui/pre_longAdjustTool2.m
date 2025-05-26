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

% Last Modified by GUIDE v2.5 04-Jan-2016 11:17:42

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
handles.codeDir = '../code';
handles.guiDir  = pwd;

% Se adquieren los datos del modelo y se reasignan en la estructura interna
% de la interfaz mdat

modelo          = getDefaultModel();

handles.mdat = getModelNfo(modelo);

% Se inicia la interfaz una vez se tienen los datos del modelo
handles = setWindowsConf(handles);

% Se inicia la grafica
[handles.surfObj, handles.centers] = initGraph(handles.mdat);
refreshGraph(handles.surfObj, handles.centers, handles.mdat);

guidata(hObject, handles);

uiwait(handles.figure1);



function modelData = getModelNfo(modelo)
% Se sacan los datos del modelo proporcionado y se almacenan en la
% estructura interna de la aplicacion
modelData.conf      = modelo.conf; 

% Wing data
modelData.w.CLa     = modelo.ala.CLa;
modelData.w.CL0     = modelo.ala.CL0;
modelData.w.CM      = modelo.ala.CM0;
modelData.w.eta     = modelo.ala.eta;
modelData.w.Zca     = modelo.ala.Zca;
modelData.w.Xca     = modelo.ala.Xca;
modelData.w.S       = modelo.ala.S;
modelData.w.b       = modelo.ala.b;
modelData.w.c       = modelo.ala.MAC;
modelData.w.AR      = modelo.ala.AR;
modelData.w.TR      = modelo.ala.TR;
modelData.w.LAM     = modelo.ala.LAM;
modelData.w.LAMc4   = modelo.ala.LAMc4;


switch modelData.conf
    case 'convencional'
    % Tail data
    modelData.t.CLa     = modelo.horizontal.CLa;
    modelData.t.CL0     = modelo.horizontal.CL0;
    modelData.t.CM      = modelo.horizontal.CM0;
    modelData.t.eta     = modelo.horizontal.eta;
    modelData.t.Zca     = modelo.horizontal.Zca;
    modelData.t.Xca     = modelo.horizontal.Xca;
    modelData.t.S       = modelo.horizontal.S;
    modelData.t.c       = modelo.horizontal.MAC;
    
    case 'ala_canard'
    % Canard data
    
    modelData.c.CLa     = modelo.canard.CLa;
    modelData.c.CL0     = modelo.canard.CL0;
    modelData.c.CM      = modelo.canard.CM0;
    modelData.c.eta     = modelo.canard.eta;
    modelData.c.Zca     = modelo.canard.Zca;
    modelData.c.Xca     = modelo.canard.Xca;
    modelData.c.S       = modelo.canard.S;
    modelData.c.c       = modelo.canard.MAC;
    
    case 'convencional_canard'
    
    modelData.t.CLa     = modelo.horizontal.CLa;
    modelData.t.CL0     = modelo.horizontal.CL0;
    modelData.t.CM      = modelo.horizontal.CM0;
    modelData.t.eta     = modelo.horizontal.eta;
    modelData.t.Zca     = modelo.horizontal.Zca;
    modelData.t.Xca     = modelo.horizontal.Xca;
    modelData.t.S       = modelo.horizontal.S;
    modelData.t.c       = modelo.horizontal.MAC;
    
    modelData.c.CLa     = modelo.canard.CLa;
    modelData.c.CL0     = modelo.canard.CL0;
    modelData.c.CM      = modelo.canard.CM0;
    modelData.c.eta     = modelo.canard.eta;
    modelData.c.Zca     = modelo.canard.Zca;
    modelData.c.Xca     = modelo.canard.Xca;
    modelData.c.S       = modelo.canard.S;
    modelData.c.c       = modelo.canard.MAC;
end

% Basic general information
modelData.nfo.mtow      = modelo.general.mtow;
modelData.nfo.w_w0      = modelo.general.w_w0;
modelData.nfo.W         = modelData.nfo.mtow*modelData.nfo.w_w0;
modelData.nfo.q         = modelo.general.qinf;
modelData.nfo.Sref      = modelo.general.Sref;
modelData.nfo.L_fus     = modelo.general.L;
modelData.nfo.Xcg       = modelo.general.Xcg;      

function handles = setWindowsConf(handles)
    % CG
    handles.mdat.nfo.Xcg_min = 0;
    handles.mdat.nfo.Xcg_max = handles.mdat.nfo.L_fus;
    set(handles.Xcg_min, 'String', num2str(handles.mdat.nfo.Xcg_min));
    set(handles.Xcg_max, 'String', num2str(handles.mdat.nfo.Xcg_max));
    set(handles.Xcg_val, 'String', num2str(handles.mdat.nfo.Xcg));
    set(handles.Xcg_slide, 'Min', handles.mdat.nfo.Xcg_min, 'Max', handles.mdat.nfo.Xcg_max, 'Value', handles.mdat.nfo.Xcg);
    
    set(handles.lfus_val, 'String', num2str(handles.mdat.nfo.L_fus));
    set(handles.c_val, 'String', num2str(handles.mdat.w.c));
    set(handles.Sref_val, 'String', num2str(handles.mdat.nfo.Sref));

    handles.mdat.w.Xmin = 0;
    handles.mdat.w.Xmax = handles.mdat.nfo.L_fus;
    if ~isfield(handles.mdat.w,'Xca') || isempty(handles.mdat.w.Xca)
        handles.mdat.w.Xca    = 0.5*(handles.mdat.w.Xmin + handles.mdat.w.Xmax);
    end
    set(handles.Xw_text, 'Visible', 'on', 'String', 'X_w');
    set(handles.Xw_slide,'Visible','on', 'Max', handles.mdat.w.Xmax,'Min',handles.mdat.w.Xmin,'Value',handles.mdat.w.Xca);
    set(handles.Xw_val, 'Visible', 'on', 'String', num2str(handles.mdat.w.Xca,3));
    set(handles.Xw_min_text, 'Visible', 'on', 'String', num2str(handles.mdat.w.Xmin));
    set(handles.Xw_max_text, 'Visible', 'on', 'String', num2str(handles.mdat.w.Xmax));
    
    set(handles.Sw_text, 'Visible', 'on', 'String', 'S_w');
    set(handles.Sw_val, 'Visible', 'on', 'String', num2str(handles.mdat.w.S));
    
    set(handles.CLaw_text, 'Visible', 'on', 'String', 'CLa_w');
    set(handles.CLaw_val, 'Visible', 'on', 'String', num2str(handles.mdat.w.CLa));
    
    set(handles.iw_text, 'Visible', 'on', 'String', 'i_w');
    set(handles.iw_val, 'Visible', 'on','Style','text');
    
switch handles.mdat.conf
    case 'convencional'
        % Posiciones
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        handles.mdat.t.Xmax = handles.mdat.w.Xmax;
        if ~isfield(handles.mdat.t,'Xca') || isempty(handles.mdat.t.Xca)
            handles.mdat.t.Xca = 0.5*(handles.mdat.t.Xmin + handles.mdat.t.Xmax);
        end
        set(handles.Xh_text, 'Visible', 'on','String', 'X_h');
        set(handles.Xh_slide,'Visible','on','Max',handles.mdat.t.Xmax,'Min',handles.mdat.t.Xmin,'Value',handles.mdat.t.Xca);
        set(handles.Xh_val, 'Visible', 'on','String',num2str(handles.mdat.t.Xca,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(handles.mdat.t.Xmax));
        
        set(handles.Xc_text, 'Visible', 'off');
        set(handles.Xc_slide,'Visible','off');
        set(handles.Xc_val, 'Visible', 'off');
        set(handles.Xc_min_text,'Visible','off');
        set(handles.Xc_max_text,'Visible','off');
        
        % Superficies
        set(handles.Sh_text, 'Visible', 'on', 'String', 'S_h');
        set(handles.Sh_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.S));
        set(handles.Sc_text,'Visible', 'off');
        set(handles.Sc_val,'Visible', 'off');
        
        % Pendientes de sustentanción
        set(handles.CLah_text, 'Visible', 'on', 'String', 'CLa_h');
        set(handles.CLah_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.CLa));
        set(handles.CLac_text,'Visible', 'off');
        set(handles.CLac_val,'Visible', 'off');
        
        
        % Incidencias
        set(handles.i_buttonGroup, 'Visible', 'off')
        set(handles.ih_text, 'Visible', 'on', 'String', 'i_t');
        set(handles.ih_val, 'Visible', 'on','Style','text');
        set(handles.ic_text, 'Visible', 'off');
        set(handles.ic_val, 'Visible', 'off');
        
        
    case 'ala_canard'
        
        % Posiciones
        handles.mdat.c.Xmin = handles.mdat.w.Xmin;
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if ~isfield(handles.mdat.c,'Xca') || isempty(handles.mdat.c.Xca)
            handles.mdat.c.Xca = 0.5*(handles.mdat.c.Xmin + handles.mdat.c.Xmax);
        end
        set(handles.Xh_text, 'Visible', 'on','String', 'X_c');
        set(handles.Xh_slide,'Visible','on','Max',handles.mdat.c.Xmax,'Min',handles.mdat.c.Xmin,'Value',handles.mdat.c.Xca);
        set(handles.Xh_val, 'Visible', 'on','String',num2str(handles.mdat.c.Xca,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(handles.mdat.c.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(handles.mdat.c.Xmax));
        
        
        set(handles.Xc_text, 'Visible', 'off');
        set(handles.Xc_slide,'Visible','off');
        set(handles.Xc_val, 'Visible', 'off');
        set(handles.Xc_min_text,'Visible','off');
        set(handles.Xc_max_text,'Visible','off');
        
        % Superficies
        set(handles.Sh_text, 'Visible', 'on', 'String', 'S_c');
        set(handles.Sh_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.S));
        set(handles.Sc_text,'Visible', 'off');
        set(handles.Sc_val,'Visible', 'off');
        
        % Pendientes de sustentanción
        set(handles.CLah_text, 'Visible', 'on', 'String', 'CLa_c');
        set(handles.CLah_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.CLa));
        set(handles.CLac_text,'Visible', 'off');
        set(handles.CLac_val,'Visible', 'off');
        
        
        % Incidencias
        set(handles.i_buttonGroup, 'Visible', 'off')
        set(handles.ih_text, 'Visible', 'on', 'String', 'i_c');
        set(handles.ih_val, 'Visible', 'on','Style','text');
        set(handles.ic_text, 'Visible', 'off');
        set(handles.ic_val, 'Visible', 'off');
        
    case 'convencional_canard'
        handles.mdat.conf
        % Posiciones
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        handles.mdat.t.Xmax = handles.mdat.w.Xmax;
        if ~isfield(handles.mdat.t,'Xca') || isempty(handles.mdat.t.Xca)
            handles.mdat.t.Xca = 0.5*(handles.mdat.t.Xmin + handles.mdat.t.Xmax);
        end
        
        handles.mdat.c.Xmin = handles.mdat.w.Xmin;
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if ~isfield(handles.mdat.c,'Xca') || isempty(handles.mdat.c.Xca)
            handles.mdat.c.Xca = 0.5*(handles.mdat.c.Xmin + handles.mdat.c.Xmax);
        end
        
        set(handles.Xh_text, 'Visible', 'on','String', 'X_h');
        set(handles.Xh_slide,'Visible','on','Max',handles.mdat.t.Xmax,'Min',handles.mdat.t.Xmin,'Value',handles.mdat.t.Xca);
        set(handles.Xh_val, 'Visible', 'on','String',num2str(handles.mdat.t.Xca,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(handles.mdat.t.Xmax));
        
        set(handles.Xc_text, 'Visible', 'on','String', 'X_c');
        set(handles.Xc_slide,'Visible','on','Max',handles.mdat.c.Xmax,'Min',handles.mdat.c.Xmin,'Value',handles.mdat.c.Xca);
        set(handles.Xc_val, 'Visible', 'on','String',num2str(handles.mdat.c.Xca,3));
        set(handles.Xc_min_text,'Visible','on','String', num2str(handles.mdat.c.Xmin));
        set(handles.Xc_max_text,'Visible','on','String', num2str(handles.mdat.c.Xmax));
        
        
        % Superficies
        set(handles.Sh_text, 'Visible', 'on', 'String', 'S_h');
        set(handles.Sh_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.S));
        
        set(handles.Sc_text, 'Visible', 'on', 'String', 'S_c');
        set(handles.Sc_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.S));
        
        
        % Pendientes de sustentanción
        set(handles.CLah_text, 'Visible', 'on', 'String', 'CLa_h');
        set(handles.CLah_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.CLa));
        
        set(handles.CLac_text, 'Visible', 'on', 'String', 'CLa_c');
        set(handles.CLac_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.CLa));
        
        
        % Incidencias
        set(handles.i_buttonGroup, 'Visible', 'on')
        set(handles.iw_button,'Value', 1);
        handles.mdat.nfo.i_fix = 'w';
        set(handles.iw_val,'Style', 'Edit');
        set(handles.ih_text, 'Visible', 'on', 'String', 'i_t');
        set(handles.ih_val, 'Visible', 'on','Style','text');
        set(handles.ic_text, 'Visible', 'on', 'String', 'i_c');
        set(handles.ic_val, 'Visible', 'on','Style','text');   
end


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

% Dependiendo de la configuracion este slider estara relacionado con el
% canard o con el horizontal de cola

% Actualizacion de los datos internos de la aplicacion
switch handles.mdat.conf % Switch entre las posibles configuraciones
    case {'convencional', 'convencional_canard'}
        handles.mdat.t.Xca = get(handles.Xh_slide,'Value');
        set(handles.Xh_val,'String', num2str(handles.mdat.t.Xca,3));
    case 'ala_canard'
        handles.mdat.c.Xca = get(handles.Xh_slide,'Value');
        set(handles.Xh_val,'String', num2str(handles.mdat.c.Xca,3));
end

% Ejecucion del calculo del punto neutro y almacenamiento de dicha variable
% en la estructura correspondiente. Actualizacion de los campos de texto de
% la aplicacion.
cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

% Actualizacion de la grafica
refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

handles.mdat.w.Xca = get(handles.Xw_slide,'Value');
set(handles.Xw_val,'String',num2str(handles.mdat.w.Xca,3));

switch handles.mdat.conf
    case 'convencional' % Tail
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        if handles.mdat.t.Xca < handles.mdat.w.Xca
            handles.mdat.t.Xca = handles.mdat.w.Xca;
        end
        set(handles.Xh_val, 'String', num2str(handles.mdat.t.Xca,3))
        set(handles.Xh_slide,'Value',handles.mdat.t.Xca,'Min',handles.mdat.t.Xmin);
        set(handles.Xh_min_text, 'String', num2str(handles.mdat.t.Xmin,2));
    
    case 'ala_canard' % Canard
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if handles.mdat.c.Xca > handles.mdat.w.Xca
            handles.mdat.c.Xca = handles.mdat.w.Xca;
        end
        
        set(handles.Xh_val, 'String', num2str(handles.mdat.c.Xca,3));
        set(handles.Xh_slide,'Value',handles.mdat.c.Xca,'Max',handles.mdat.c.Xmax);
        set(handles.Xh_max_text, 'String', num2str(handles.mdat.c.Xmax,2));
    
    case 'convencional_canard' % Canard + Tail
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        if handles.mdat.t.Xca < handles.mdat.w.Xca;
            handles.mdat.t.Xca = handles.mdat.w.Xca;
        end
        
        set(handles.Xh_val, 'String', num2str(handles.mdat.t.Xca,3));
        set(handles.Xh_slide,'Value',handles.mdat.t.Xca,'Min',handles.mdat.t.Xmin);
        set(handles.Xh_min_text, 'String', num2str(handles.mdat.t.Xmin,2));        

        
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if handles.mdat.c.Xca > handles.mdat.w.Xca
            handles.mdat.c.Xca = handles.mdat.w.Xca;
        end

        set(handles.Xc_val, 'String', num2str(handles.mdat.c.Xca,3));
        set(handles.Xc_slide,'Value',handles.mdat.c.Xca,'Max',handles.mdat.c.Xmax);
        set(handles.Xc_max_text, 'String', num2str(handles.mdat.c.Xmax,2));
        
end

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

% Este slider esta vinculado al canard en caso de que se tenga una
% configuracion convencional+canard
handles.mdat.c.Xca = get(handles.Xc_slide,'Value'); 
set(handles.Xc_val,'String', num2str(handles.mdat.c.Xca,3)); 

% Calculo del punto neutro y actualizacion de los campos de texto
% correspondientes
cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

% Refresco de la grafica
refreshGraph(handles.surfObj, handles.centers, handles.mdat);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Xc_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xc_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in calcula_button.
function calcula_button_Callback(hObject, eventdata, handles)
% hObject    handle to calcula_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Resolucion del sistema de ecuaciones y obtencion de las incidencias


switch handles.mdat.conf
    case 'convencional' % Convencional
        cd(handles.codeDir);
        [i_sol] = pre_equilibrio(handles.mdat);
        cd(handles.guiDir)
        
        handles.mdat.w.i = i_sol(1);
        handles.mdat.t.i = i_sol(2);
        set(handles.iw_val, 'String', num2str(handles.mdat.w.i*180/pi,4));
        set(handles.ih_val, 'String', num2str(handles.mdat.t.i*180/pi,4));
        
    case 'ala_canard' % Ala + canard
        cd(handles.codeDir);
        [i_sol] = pre_equilibrio(handles.mdat);
        cd(handles.guiDir)
        
        handles.mdat.w.i = i_sol(1);
        handles.mdat.c.i = i_sol(2);
        set(handles.iw_val, 'String', num2str(handles.mdat.w.i*180/pi,4));
        set(handles.ih_val, 'String', num2str(handles.mdat.c.i*180/pi,4));
               
    case 'convencional_canard' % Convencional + canard
        if isfield(handles.mdat.(handles.mdat.nfo.i_fix), 'i') && ~isempty(handles.mdat.(handles.mdat.nfo.i_fix).i)
            cd(handles.codeDir);
            [i_sol] = pre_equilibrio(handles.mdat);
            cd(handles.guiDir)
        switch handles.mdat.nfo.i_fix % Switch entre la incidencia fijada
            case 'w'
                handles.mdat.t.i = i_sol(1);
                handles.mdat.c.i = i_sol(2);   
            case 'c'
                handles.mdat.w.i = i_sol(1);
                handles.mdat.t.i = i_sol(2);
            case 't' 
                handles.mdat.w.i = i_sol(1);
                handles.mdat.c.i = i_sol(2);
        end
        set(handles.iw_val, 'String', num2str(handles.mdat.w.i*180/pi,4));
        set(handles.ih_val, 'String', num2str(handles.mdat.t.i*180/pi,4));
        set(handles.ic_val, 'String', num2str(handles.mdat.c.i*180/pi,4));
        end
        
end

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
guidata(hObject,handles)



function handles = Xw_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xw_val as text
%        str2double(get(hObject,'String')) returns contents of Xw_val as a double
%str2double(get(hObject,'String'))

% Actualizacion de la informacion interna de la aplicacion
handles.mdat.w.Xca = str2double(get(handles.Xw_val,'String'));

if handles.mdat.w.Xca < handles.mdat.w.Xmin || isnan(handles.mdat.w.Xca)
    handles.mdat.x.X = handles.mdat.w.Xmin; 
elseif handles.mdat.w.Xca > handles.mdat.w.Xmax
    handles.mdat.w.Xca = handles.mdat.w.Xmax;
end

set(handles.Xw_val,'String', num2str(handles.mdat.w.Xca));
set(handles.Xw_slide,'Value', handles.mdat.w.Xca);

switch handles.mdat.conf % Switch entre las configuraciones
    case 'convencional' % Tail
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        
        if handles.mdat.t.Xca <= handles.mdat.t.Xmin
            handles.mdat.t.Xca = handles.mdat.t.Xmin;
        end

        set(handles.Xh_min_text,'String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_slide, 'Value', handles.mdat.t.Xca,'Min',handles.mdat.t.Xmin);
        set(handles.Xh_val, 'String', num2str(handles.mdat.t.Xca));

    case 'ala_canard' % Canard
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        
        if handles.mdat.c.Xca > handles.mdat.c.Xmax
            handles.mdat.c.Xca = handles.mdat.c.Xmax;
        end

        set(handles.Xh_max_text,'String', num2str(handles.mdat.c.Xmax));
        set(handles.Xh_slide, 'Value', handles.mdat.c.Xca, 'Max', handles.mdat.c.Xmax);
        set(handles.Xh_val, 'String', num2str(handles.mdat.c.Xca));
        
    case 'convencional_canard' % Convencional+canard
        % Tail
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        
        if handles.mdat.t.Xca <= handles.mdat.t.Xmin
            handles.mdat.t.Xca = handles.mdat.t.Xmin;
        end

        set(handles.Xh_min_text,'String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_slide, 'Value', handles.mdat.t.Xca, 'Min', handles.mdat.t.Xmin);
        set(handles.Xh_val, 'String', num2str(handles.mdat.t.Xca));
        
        % Canard
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if handles.mdat.c.Xca > handles.mdat.c.Xmax
            handles.mdat.c.Xca = handles.mdat.c.Xmax;
        end

        set(handles.Xc_max_text,'String', num2str(handles.mdat.c.Xmax));
        set(handles.Xc_slide, 'Value', handles.mdat.c.Xca, 'Max', handles.mdat.c.Xmax);
        set(handles.Xc_val, 'String', num2str(handles.mdat.c.Xca));
end

% Calculo del centro aerodinamico de la aeronave y actualizacion de los
% campos
cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

% Actualizacion de las graficas
refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

handles.mdat.c.Xca = str2double(get(handles.Xc_val,'String'));

if handles.mdat.c.Xca < handles.mdat.c.Xmin || isnan(handles.mdat.c.Xca)
    handles.mdat.c.Xca = handles.mdat.c.Xmin;
elseif handles.mdat.c.Xca > handles.mdat.c.Xmax
    handles.mdat.c.Xca = handles.mdat.c.Xmax;
end

set(handles.Xc_val,'String', num2str(handles.mdat.c.Xca));
set(handles.Xc_slide,'Value', handles.mdat.c.Xca);

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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


function handles = Xh_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xh_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xh_val as text
%        str2double(get(hObject,'String')) returns contents of Xh_val as a double

switch handles.mdat.conf 
    case {'convencional', 'convencional_canard'} % Tail || Canard+tail
        handles.mdat.t.Xca = str2double(get(handles.Xh_val,'String'));
        
        if handles.mdat.t.Xca < handles.mdat.t.Xmin || isnan(handles.mdat.t.Xca)
            handles.mdat.t.Xca = handles.mdat.t.Xmin;
        elseif handles.mdat.t.Xca > handles.mdat.t.Xmax
            handles.mdat.t.Xca = handles.mdat.t.Xmax;
        end
        
        set(handles.Xh_val,'String', num2str(handles.mdat.t.Xca));
        set(handles.Xh_slide,'Value', handles.mdat.t.Xca);
        
    case 'ala_canard' % Canard
        handles.mdat.c.Xca = str2double(get(handles.Xh_val,'String'));
        
        if handles.mdat.c.Xca < handles.mdat.c.Xmin || isnan(handles.mdat.c.Xca)
            handles.mdat.c.Xca = handles.mdat.c.Xmin;
        elseif handles.mdat.c.Xca > handles.mdat.c.Xmax
            handles.mdat.c.Xca = handles.mdat.c.Xmax;
        end

        set(handles.Xh_val,'String', num2str(handles.mdat.c.Xca));
        set(handles.Xh_slide,'Value', handles.mdat.c.Xca);

end

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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


% --- Executes on button press in iw_button.
function iw_button_Callback(hObject, eventdata,handles)
% hObject    handle to iw_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iw_button
set(handles.iw_val,'Style','edit','String','');
set(handles.ih_val,'Style','text','String','');
set(handles.ic_val,'Style','text','String','');
handles.mdat.nfo.i_fix = 'w';

guidata(hObject,handles);


% --- Executes on button press in ih_button.
function ih_button_Callback(hObject, eventdata,handles)
% hObject    handle to ih_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ih_button
set(handles.iw_val,'Style','text','String','');
set(handles.ih_val,'Style','edit','String','');
set(handles.ic_val,'Style','text','String','');
handles.mdat.nfo.i_fix = 't';

guidata(hObject,handles)


% --- Executes on button press in ic_button.
function ic_button_Callback(hObject,eventdata,handles)
% hObject    handle to ic_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ic_button
set(handles.iw_val,'Style','text','String','');
set(handles.ih_val,'Style','text','String','');
set(handles.ic_val,'Style','edit','String','');
handles.mdat.nfo.i_fix = 'c';

guidata(hObject,handles)



function iw_val_Callback(hObject, eventdata, handles)
% hObject    handle to iw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iw_val as text
%        str2double(get(hObject,'String')) returns contents of iw_val as a double
handles.mdat.w.i = str2double(get(handles.iw_val,'String'))*pi/180;

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function iw_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ih_val_Callback(hObject, eventdata, handles)
% hObject    handle to ih_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ih_val as text
%        str2double(get(hObject,'String')) returns contents of ih_val as a double
switch handles.mdat.conf
    case {'convencional', 'convencional_canard'}
        handles.mdat.t.i = str2double(get(handles.ih_val,'String'))*pi/180;
    case 'ala_canard'
        handles.mdat.c.i = str2double(get(handles.ih_val,'String'))*pi/180;        
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ih_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ih_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ic_val_Callback(hObject, eventdata, handles)
% hObject    handle to ic_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ic_val as text
%        str2double(get(hObject,'String')) returns contents of ic_val as a double

handles.mdat.c.i = str2double(get(handles.ic_val,'String'))*pi/180;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ic_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ic_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lfus_val_Callback(hObject, eventdata, handles)
% hObject    handle to lfus_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lfus_val as text
%        str2double(get(hObject,'String')) returns contents of lfus_val as a double

handles.mdat.nfo.L_fus = str2double(get(hObject,'String'));
handles.mdat.nfo.Xcg_max = handles.mdat.nfo.L_fus;
if handles.mdat.nfo.Xcg > handles.mdat.nfo.Xcg_max
    handles.mdat.nfo.Xcg = handles.mdat.nfo.Xcg_max;
end
set(handles.Xcg_slide, 'Max', handles.mdat.nfo.Xcg_max, 'Value', handles.mdat.nfo.Xcg);
set(handles.Xcg_val,'String', num2str(handles.mdat.nfo.Xcg));
set(handles.Xcg_max, 'String', num2str(handles.mdat.nfo.Xcg_max));

handles.mdat.w.Xmax = handles.mdat.nfo.L_fus;
if handles.mdat.w.Xca > handles.mdat.w.Xmax
    handles.mdat.w.Xca = handles.mdat.w.Xmax;
end
set(handles.Xw_slide, 'Max', handles.mdat.w.Xmax, 'Value', handles.mdat.w.Xca);
set(handles.Xw_val,'String', num2str(handles.mdat.w.Xca));
set(handles.Xw_max_text, 'String', num2str(handles.mdat.w.Xmax));

switch handles.mdat.conf
    case 'convencional'
        if handles.mdat.t.Xca > handles.mdat.w.Xmax
            handles.mdat.t.Xca = handles.mdat.w.Xmax;
        end
        handles.mdat.t.Xmax = handles.mdat.w.Xmax;
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        set(handles.Xh_min_text, 'String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_max_text, 'String', num2str(handles.mdat.t.Xmax));
        
    case 'ala_canard'
        if handles.mdat.c.Xca > handles.mdat.w.Xca
            handles.mdat.c.Xca = handles.mdat.w.Xca;
        end
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        set(handles.Xh_max_text, 'String', num2str(handles.mdat.c.Xmax));
    case 'convencional_canard'
        if handles.mdat.t.Xca > handles.mdat.w.Xmax
            handles.mdat.t.Xca = handles.mdat.w.Xmax;
        end
        handles.mdat.t.Xmax = handles.mdat.w.Xmax;
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        set(handles.Xh_min_text, 'String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_slide, 'Min', handles.mdat.t.Xmin);
        set(handles.Xh_max_text, 'String', num2str(handles.mdat.t.Xmax));
        set(handles.Xh_slide, 'Max', handles.mdat.t.Xmax, 'Value', handles.mdat.t.Xca);
        set(handles.Xh_val, 'String', num2str(handles.mdat.t.Xca));
        
        if handles.mdat.c.Xca > handles.mdat.w.Xca
            handles.mdat.c.Xca = handles.mdat.w.Xca;
        end
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        set(handles.Xc_max_text, 'String', num2str(handles.mdat.c.Xmax));
        set(handles.Xc_slide, 'Value', handles.mdat.c.Xca, 'Max', handles.mdat.c.Xmax);
        set(handles.Xc_val, 'String', num2str(handles.mdat.c.Xca));
        
end

if ~isfield(handles.mdat.nfo, 'Xna')
    cd(handles.codeDir);
    handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
    cd(handles.guiDir);
end
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
guidata(hObject,handles)


function Xcg_val_Callback(hObject, eventdata, handles)
% hObject    handle to Xcg_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xcg_val as text
%        str2double(get(hObject,'String')) returns contents of Xcg_val as a double

handles.mdat.nfo.Xcg = str2double(get(hObject,'String'));
if handles.mdat.nfo.Xcg > handles.mdat.nfo.Xcg_max
    handles.mdat.nfo.Xcg = handles.mdat.nfo.Xcg_max;
    set(hObject, 'String', num2str(handles.mdat.nfo.Xcg));
elseif handles.mdat.nfo.Xcg < handles.mdat.nfo.Xcg_min
    handles.mdat.nfo.Xcg = handles.mdat.nfo.Xcg_min;
    set(hObject, 'String', num2str(handles.mdat.nfo.Xcg));
end
set(handles.Xcg_slide, 'Value', handles.mdat.nfo.Xcg);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);

guidata(hObject,handles);


% --- Executes on slider movement.
function Xcg_slide_Callback(hObject, eventdata, handles)
% hObject    handle to Xcg_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.mdat.nfo.Xcg = get(hObject,'Value');
set(handles.Xcg_val,'String', num2str(handles.mdat.nfo.Xcg));


handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);

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


% --- Executes on button press in import_button.
function import_button_Callback(hObject, eventdata, handles)
% hObject    handle to import_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.FileName, handles.PathName, ~] = uigetfile('*.mat','Seleccione modelo','../models');
handles.model2import = fullfile(handles.PathName, handles.FileName);
loadedData = load(handles.model2import);
modelName = char(fieldnames(loadedData));
handles.model = loadedData.(modelName);

handles.mdat = getModelNfo(handles.model);
handles = setWindowsConf(handles);

guidata(hObject,handles)

function CLaw_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLaw_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLaw_val as text
%        str2double(get(hObject,'String')) returns contents of CLaw_val as a double

handles.mdat.w.CLa = str2double(get(hObject,'String'));
cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

switch handles.mdat.conf
    case {'convencional', 'convencional_canard'}
        handles.mdat.t.CLa = str2double(get(hObject,'String'));
    case 'ala_canard'
        handles.mdat.c.CLa = str2double(get(hObject,'String'));
end

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

handles.mdat.c.CLa = str2double(get(hObject,'String'));

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

handles.mdat.w.S = str2double(get(hObject,'String'));
cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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

switch handles.mdat.conf
    case {'convencional', 'convencional_canard'}
        handles.mdat.t.S = str2double(get(hObject,'String'));
    case 'ala_canard'
        handles.mdat.c.S = str2double(get(hObject,'String'));
end

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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
handles.mdat.c.S = str2double(get(hObject,'String'));

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
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


% --- Executes on button press in conv_button.
function conv_button_Callback(hObject, eventdata, handles)
% hObject    handle to conv_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conv_button

handles.mdat.conf = 'convencional';

handles = refreshWindowsConf(handles);

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);
% 
handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
guidata(hObject,handles)


% --- Executes on button press in canard_button.
function canard_button_Callback(hObject, eventdata, handles)
% hObject    handle to canard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of canard_button
handles.mdat.conf = 'ala_canard';

handles = refreshWindowsConf(handles);

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
guidata(hObject,handles)


% --- Executes on button press in convCanard_button.
function convCanard_button_Callback(hObject, eventdata, handles)
% hObject    handle to convCanard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convCanard_button
handles.mdat.conf = 'convencional_canard';

handles = refreshWindowsConf(handles);

cd(handles.codeDir);
handles.mdat.nfo.Xna = pre_aeroCenter(handles.mdat);
cd(handles.guiDir);

handles.mdat.nfo.ME = 1/handles.mdat.w.c*(handles.mdat.nfo.Xna - handles.mdat.nfo.Xcg);
set(handles.Xna_val,'String', num2str(handles.mdat.nfo.Xna,3));
set(handles.ME_val,'String', num2str(handles.mdat.nfo.ME,3));

refreshGraph(handles.surfObj, handles.centers, handles.mdat);
guidata(hObject,handles)


% --- Executes on button press in extraData_button.
function extraData_button_Callback(hObject, eventdata, handles)
% hObject    handle to extraData_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
guidata(hObject,handles)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
guidata(hObject,handles)



function handles = refreshWindowsConf(handles)
    
switch handles.mdat.conf
    case 'convencional'
        % Posiciones
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        handles.mdat.t.Xmax = handles.mdat.w.Xmax;
        
        handles.mdat.t.Xca = 0.5*(handles.mdat.t.Xmin + handles.mdat.t.Xmax);
        
        set(handles.Xh_text, 'Visible', 'on','String', 'X_h');
        set(handles.Xh_slide,'Visible','on','Max',handles.mdat.t.Xmax,'Min',handles.mdat.t.Xmin,'Value',handles.mdat.t.Xca);
        set(handles.Xh_val, 'Visible', 'on','String',num2str(handles.mdat.t.Xca,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(handles.mdat.t.Xmax));
        
        set(handles.Xc_text, 'Visible', 'off');
        set(handles.Xc_slide,'Visible','off');
        set(handles.Xc_val, 'Visible', 'off');
        set(handles.Xc_min_text,'Visible','off');
        set(handles.Xc_max_text,'Visible','off');
        
        % Superficies
        set(handles.Sh_text, 'Visible', 'on', 'String', 'S_h');
        if isfield(handles.mdat,'t') && isfield(handles.mdat.t,'S')  && ~isempty(handles.mdat.t.S)
            set(handles.Sh_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.S));
        else
            set(handles.Sh_val, 'Visible', 'on', 'String', '')
        end
        set(handles.Sc_text,'Visible', 'off');
        set(handles.Sc_val,'Visible', 'off');
        
        % Pendientes de sustentanción
        set(handles.CLah_text, 'Visible', 'on', 'String', 'CLa_h');
        if isfield(handles.mdat,'t') && isfield(handles.mdat.t,'CLa')  && ~isempty(handles.mdat.t.CLa)
            set(handles.CLah_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.CLa));
        else
            set(handles.CLah_val, 'Visible', 'on', 'String', '')
        end
        set(handles.CLac_text,'Visible', 'off');
        set(handles.CLac_val,'Visible', 'off');
        
        
        % Incidencias
        set(handles.i_buttonGroup, 'Visible', 'off')
        set(handles.ih_text, 'Visible', 'on', 'String', 'i_t');
        set(handles.ih_val, 'Visible', 'on','Style','text');
        set(handles.ic_text, 'Visible', 'off');
        set(handles.ic_val, 'Visible', 'off');
        
        
    case 'ala_canard'
        
        % Posiciones
        handles.mdat.c.Xmin = handles.mdat.w.Xmin;
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if ~isfield(handles.mdat.c,'Xca') || isempty(handles.mdat.c.Xca)
            handles.mdat.c.Xca = 0.5*(handles.mdat.c.Xmin + handles.mdat.c.Xmax);
        end
        set(handles.Xc_text, 'Visible', 'on','String', 'X_c');
        set(handles.Xc_slide,'Visible','on','Max',handles.mdat.c.Xmax,'Min',handles.mdat.c.Xmin,'Value',handles.mdat.c.Xca);
        set(handles.Xc_val, 'Visible', 'on','String',num2str(handles.mdat.c.Xca,3));
        set(handles.Xc_min_text,'Visible','on','String', num2str(handles.mdat.c.Xmin));
        set(handles.Xc_max_text,'Visible','on','String', num2str(handles.mdat.c.Xmax));
        
        
        set(handles.Xh_text, 'Visible', 'off');
        set(handles.Xh_slide,'Visible','off');
        set(handles.Xh_val, 'Visible', 'off');
        set(handles.Xh_min_text,'Visible','off');
        set(handles.Xh_max_text,'Visible','off');
        
        % Superficies
        set(handles.Sh_text, 'Visible', 'on', 'String', 'S_c');
        if isfield(handles.mdat,'c') && isfield(handles.mdat.c,'CLa')  && ~isempty(handles.mdat.c.CLa)
            set(handles.Sh_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.S));
        else
            set(handles.Sh_val, 'Visible', 'on', 'String', '')
        end
        set(handles.Sc_text,'Visible', 'off');
        set(handles.Sc_val,'Visible', 'off');
        
        % Pendientes de sustentanción
        set(handles.CLah_text, 'Visible', 'on', 'String', 'CLa_c');
        if isfield(handles.mdat,'c') && isfield(handles.mdat.c,'CLa')  && ~isempty(handles.mdat.c.CLa)
            set(handles.CLah_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.CLa));
        else
            set(handles.CLah_val, 'Visible', 'on', 'String', '')
        end
        
        set(handles.CLac_text,'Visible', 'off');
        set(handles.CLac_val,'Visible', 'off');
        
        
        % Incidencias
        set(handles.i_buttonGroup, 'Visible', 'off')
        set(handles.ih_text, 'Visible', 'on', 'String', 'i_c');
        set(handles.ih_val, 'Visible', 'on','Style','text');
        set(handles.ic_text, 'Visible', 'off');
        set(handles.ic_val, 'Visible', 'off');
        
    case 'convencional_canard'
        
        % Posiciones
        handles.mdat.t.Xmin = handles.mdat.w.Xca;
        handles.mdat.t.Xmax = handles.mdat.w.Xmax;
        if ~isfield(handles.mdat.t,'Xca') || isempty(handles.mdat.t.Xca)
            handles.mdat.t.Xca = 0.5*(handles.mdat.t.Xmin + handles.mdat.t.Xmax);
        end
        
        handles.mdat.c.Xmin = handles.mdat.w.Xmin;
        handles.mdat.c.Xmax = handles.mdat.w.Xca;
        if ~isfield(handles.mdat.c,'Xca') || isempty(handles.mdat.c.Xca)
            handles.mdat.c.Xca = 0.5*(handles.mdat.c.Xmin + handles.mdat.c.Xmax);
        end
        
        
        
        set(handles.Xh_text, 'Visible', 'on','String', 'X_h');
        set(handles.Xh_slide,'Visible','on','Max',handles.mdat.t.Xmax,'Min',handles.mdat.t.Xmin,'Value',handles.mdat.t.Xca);
        set(handles.Xh_val, 'Visible', 'on','String',num2str(handles.mdat.t.Xca,3));
        set(handles.Xh_min_text,'Visible','on','String', num2str(handles.mdat.t.Xmin));
        set(handles.Xh_max_text,'Visible','on','String', num2str(handles.mdat.t.Xmax));
        
        set(handles.Xc_text, 'Visible', 'on','String', 'X_c');
        set(handles.Xc_slide,'Visible','on','Max',handles.mdat.c.Xmax,'Min',handles.mdat.c.Xmin,'Value',handles.mdat.c.Xca);
        set(handles.Xc_val, 'Visible', 'on','String',num2str(handles.mdat.c.Xca,3));
        set(handles.Xc_min_text,'Visible','on','String', num2str(handles.mdat.c.Xmin));
        set(handles.Xc_max_text,'Visible','on','String', num2str(handles.mdat.c.Xmax));
        
        
        % Superficies
        set(handles.Sh_text, 'Visible', 'on', 'String', 'S_h');
        if isfield(handles.mdat,'t') && isfield(handles.mdat.t,'S')  && ~isempty(handles.mdat.t.S)
            set(handles.Sh_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.S));
        else
            set(handles.Sh_val, 'Visible', 'on', 'String', '')
        end
        
        set(handles.Sc_text, 'Visible', 'on', 'String', 'S_c');
        if isfield(handles.mdat,'c') && isfield(handles.mdat.c,'CLa')  && ~isempty(handles.mdat.c.CLa)
            set(handles.Sc_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.S));
        else
            set(handles.Sc_val, 'Visible', 'on', 'String', '')
        end
        
        % Pendientes de sustentanción
        set(handles.CLah_text, 'Visible', 'on', 'String', 'CLa_h');
        if isfield(handles.mdat,'t') && isfield(handles.mdat.t,'CLa')  && ~isempty(handles.mdat.t.CLa)
            set(handles.CLah_val, 'Visible', 'on', 'String', num2str(handles.mdat.t.CLa));
        else
            set(handles.CLah_val, 'Visible', 'on', 'String', '')
        end
        
        set(handles.CLac_text, 'Visible', 'on', 'String', 'CLa_c');
        if isfield(handles.mdat,'c') && isfield(handles.mdat.c,'CLa')  && ~isempty(handles.mdat.c.CLa)
            set(handles.CLac_val, 'Visible', 'on', 'String', num2str(handles.mdat.c.CLa));
        else
            set(handles.CLac_val, 'Visible', 'on', 'String', '')
        end
        
        % Incidencias
        set(handles.i_buttonGroup, 'Visible', 'on')
        set(handles.iw_button,'Value', 1);
        handles.mdat.nfo.i_fix = 'w';
        set(handles.iw_val,'Style', 'Edit');
        set(handles.ih_text, 'Visible', 'on', 'String', 'i_t');
        set(handles.ih_val, 'Visible', 'on','Style','text');
        set(handles.ic_text, 'Visible', 'on', 'String', 'i_c');
        set(handles.ic_val, 'Visible', 'on','Style','text');   
end


% --- Executes on button press in flyWing_button.
function flyWing_button_Callback(hObject, eventdata, handles)
% hObject    handle to flyWing_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flyWing_button



function fixInc_text_Callback(hObject, eventdata, handles)
% hObject    handle to fixInc_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixInc_text as text
%        str2double(get(hObject,'String')) returns contents of fixInc_text as a double


% --- Executes during object creation, after setting all properties.
function fixInc_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixInc_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
