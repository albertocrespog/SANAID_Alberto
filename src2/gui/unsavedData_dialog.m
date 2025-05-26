function varargout = unsavedData_dialog(varargin)
% UNSAVEDDATA_DIALOG Application M-file for untitled.fig
%   UNSAVEDDATA_DIALOG, by itself, creates a new UNSAVEDDATA_DIALOG or raises the existing
%   singleton*.
%
%   H = UNSAVEDDATA_DIALOG returns the handle to a new UNSAVEDDATA_DIALOG or the handle to
%   the existing singleton*.
%
%   UNSAVEDDATA_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in UNSAVEDDATA_DIALOG.M with the given input arguments.
%
%   UNSAVEDDATA_DIALOG('Property','Value',...) creates a new UNSAVEDDATA_DIALOG or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before modaldlg_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to unsavedData_dialog_OpeningFcn via varargin.
%
%   *See GUI Options - GUI allows only one instance to run (singleton).
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Copyright 2000-2006 The MathWorks, Inc.

% Last Modified by GUIDE v2.5 04-Aug-2015 00:10:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',          mfilename, ...
                   'gui_Singleton',     gui_Singleton, ...
                   'gui_OpeningFcn',    @unsavedData_dialog_OpeningFcn, ...
                   'gui_OutputFcn',     @unsavedData_dialog_OutputFcn, ...
                   'gui_LayoutFcn',     [], ...
                   'gui_Callback',      []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before unsavedData_dialog is made visible.
function unsavedData_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to unsavedData_dialog (see VARARGIN)

% Choose default command line output for unsavedData_dialog
handles.output = 'CANCEL';
handles.dir = varargin{1};
handles.warning_img=imread([handles.dir.images,'warningImage.png']);
axes(handles.axes1);
axis off
imshow(handles.warning_img);
pos = get(0,'PointerLocation');
set(handles.figure1,'Position',[pos(1)-380,pos(2)-25,450,110]);
% Update handles structure
guidata(hObject, handles);
    
% UIWAIT makes unsavedData_dialog wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = unsavedData_dialog_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in discard_button.
function discard_button_Callback(hObject, eventdata, handles)
% hObject    handle to discard_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = get(hObject,'String');

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

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" - do uiresume if we get it
get(hObject,'CurrentKey')
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    


% --- Executes during object creation, after setting all properties.
function question_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to question_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



