function varargout = noName_dialog(varargin)
% NONAME_DIALOG MATLAB code for noName_dialog.fig
%      NONAME_DIALOG by itself, creates a new NONAME_DIALOG or raises the
%      existing singleton*.
%
%      H = NONAME_DIALOG returns the handle to a new NONAME_DIALOG or the handle to
%      the existing singleton*.
%
%      NONAME_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NONAME_DIALOG.M with the given input arguments.
%
%      NONAME_DIALOG('Property','Value',...) creates a new NONAME_DIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before noName_dialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to noName_dialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help noName_dialog

% Last Modified by GUIDE v2.5 06-Jan-2016 19:21:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @noName_dialog_OpeningFcn, ...
                   'gui_OutputFcn',  @noName_dialog_OutputFcn, ...
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

% --- Executes just before noName_dialog is made visible.
function noName_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to noName_dialog (see VARARGIN)

% Update handles structure
guidata(hObject, handles);

set(hObject, 'Name', 'The model needs a name');
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

% % Determine the position of the dialog - centered on the callback figure
% % if available, else, centered on the screen
% FigPos=get(0,'DefaultFigurePosition');
% OldUnits = get(hObject, 'Units');
% set(hObject, 'Units', 'pixels');
% OldPos = get(hObject,'Position');
% FigWidth = OldPos(3);
% FigHeight = OldPos(4);
% if isempty(gcbf)
%     ScreenUnits=get(0,'Units');
%     set(0,'Units','pixels');
%     ScreenSize=get(0,'ScreenSize');
%     set(0,'Units',ScreenUnits);
% 
%     FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
%     FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
% else
%     GCBFOldUnits = get(gcbf,'Units');
%     set(gcbf,'Units','pixels');
%     GCBFPos = get(gcbf,'Position');
%     set(gcbf,'Units',GCBFOldUnits);
%     FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
%                    (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
% end
% FigPos(3:4)=[FigWidth FigHeight];
% set(hObject, 'Position', FigPos);
% set(hObject, 'Units', OldUnits);
% 
% % Show a question icon from dialogicons.mat - variables questIconData
% % and questIconMap
% load dialogicons.mat
% 
% IconData=questIconData;
% questIconMap(256,:) = get(handles.figure1, 'Color');
% IconCMap=questIconMap;
% 
% Img=image(IconData, 'Parent', handles.axes1);
% set(handles.figure1, 'Colormap', IconCMap);
% 
% set(handles.axes1, ...
%     'Visible', 'off', ...
%     'YDir'   , 'reverse'       , ...
%     'XLim'   , get(Img,'XData'), ...
%     'YLim'   , get(Img,'YData')  ...
%     );

handles.dir = varargin{1};
handles.output = 1;

handles.warning_img=imread([handles.dir.images,'warningImage.png']);
axes(handles.axes1);
axis off
imshow(handles.warning_img);
pos = get(0,'PointerLocation');
set(handles.figure1,'Position',[pos(1)-380,pos(2)-25,450,110]);
% Update handles structure
guidata(hObject, handles);
    

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes noName_dialog wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = noName_dialog_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in yesButton.
function yesButton_Callback(hObject, eventdata, handles)
% hObject    handle to yesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = 1;

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in noButton.
function noButton_Callback(hObject, eventdata, handles)
% hObject    handle to noButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = 0;

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'enter')
    % User said no by hitting enter
    handles.output = 1;
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return') || isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 0;
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
