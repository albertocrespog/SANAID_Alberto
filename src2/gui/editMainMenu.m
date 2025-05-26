function varargout = editMainMenu(varargin)
% EDITMAINMENU MATLAB code for editMainMenu.fig
%      EDITMAINMENU, by itself, creates a new EDITMAINMENU or raises the existing
%      singleton*.
%
%      H = EDITMAINMENU returns the handle to a new EDITMAINMENU or the handle to
%      the existing singleton*.
%
%      EDITMAINMENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDITMAINMENU.M with the given input arguments.
%
%      EDITMAINMENU('Property','Value',...) creates a new EDITMAINMENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before editMainMenu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to editMainMenu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help editMainMenu

% Last Modified by GUIDE v2.5 07-Jan-2016 00:52:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editMainMenu_OpeningFcn, ...
                   'gui_OutputFcn',  @editMainMenu_OutputFcn, ...
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


% --- Executes just before editMainMenu is made visible.
function editMainMenu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to editMainMenu (see VARARGIN)

% Choose default command line output for editMainMenu


handles.model       = varargin{1};
handles.dir         = varargin{2};

set(hObject, 'Name', 'Edition Menu');
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

handles.changes     = 0; % No changes in the model when opening the menu
handles = button_layout(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes editMainMenu wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = editMainMenu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.model;
varargout{2} = handles.changes;
delete(handles.figure1);


function handles = button_layout(handles)
switch handles.model.conf
    case 'convencional'
        set(handles.canardEdit,'Visible','off');
        set(handles.horizEdit,'Visible','on','String', 'HORIZONTAL STABILIZER');
        %set(handles.horizEdit,'Position', [0.767 6.006 5.847 1.085]);
        set(handles.horizEdit,'Position', [0.767 3.678 5.847 1.085]);
        
        
    case 'canard'
        set(handles.horizEdit,'Visible','off');
        set(handles.canardEdit,'Visible','on');
        %set(handles.canardEdit,'Position', [0.767 6.006 5.847 1.085]);
        set(handles.canardEdit,'Position', [0.767 3.678 5.847 1.085]);
        

    case 'convencional_canard'
        set(handles.horizEdit,'Visible','on', 'String', 'HORIZONTAL STABILIZER');
        set(handles.canardEdit,'Visible','on');
%         set(handles.horizEdit,'Position', [0.767 4.419 5.847 1.058]);
%         set(handles.canardEdit,'Position', [0.767 6.006 5.847 1.085]);
        set(handles.horizEdit,'Position', [0.767 2.09 5.847 1.085]);
        set(handles.canardEdit,'Position', [0.767 3.678 5.847 1.085]);
    
    case 'flyWing'
        set(handles.horizEdit,'Visible','on', 'String', 'VIRTUAL HORIZONTAL','Position', [0.767 3.678 5.847 1.085]);
        set(handles.canardEdit,'Visible','off');
        if strcmp(handles.model.confVert,'no_vert')
            set(handles.vertEdit,'Visible','off');
        end
end



% --- Executes on button press in genEdit.
function genEdit_Callback(hObject, eventdata, handles)
% hObject    handle to genEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.model,handles.changes]   = generalData_edition(handles.model, handles.dir, handles.changes);

handles = button_layout(handles);
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in propulsionEdit.
function propulsionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to propulsionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
[handles.model.propulsion, handles.changes] = propulsiveProperties(handles.model, handles.dir, handles.changes);

guidata(hObject, handles);

% --- Executes on button press in weightEdit.
function weightEdit_Callback(hObject, eventdata, handles)
% hObject    handle to weightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.model.general,handles.changes] = weightProperties(handles.model.general,handles.dir, handles.changes);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in fusEdit.
function fusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.model.fuselaje,handles.changes] = fuselage_AutoDef(handles.model, handles.dir, handles.changes);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in wingEdit.
function wingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to wingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.model.ala = aeroSurf_prop(handles.model.name, 'ALA', handles.model.ala);
[handles.model.ala,handles.changes] = wingProperties(handles.model,handles.dir,handles.changes);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in vertEdit.
function vertEdit_Callback(hObject, eventdata, handles)
% hObject    handle to vertEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.model.vertical,handles.changes] = verticalProperties(handles.model, handles.dir, handles.changes);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in canardEdit.
function canardEdit_Callback(hObject, eventdata, handles)
% hObject    handle to canardEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.model.canard = aeroSurf_prop(handles.model.name, 'CANARD', handles.model.canard);
[handles.model.canard, handles.changes] = hor_canProperties(handles.model, 'canard', handles.dir, handles.changes);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in horizEdit.
function horizEdit_Callback(hObject, eventdata, handles)
% hObject    handle to horizEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.model.horizontal = aeroSurf_prop(handles.model.name, 'ESTABILIZADOR HORIZONTAL', handles.model.horizontal);
[handles.model.horizontal,handles.changes] = hor_canProperties(handles.model, 'horizontal', handles.dir, handles.changes);

% Update handles structure
guidata(hObject, handles);


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
