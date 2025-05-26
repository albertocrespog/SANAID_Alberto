function varargout = stabilityDerivatives_calc(varargin)
% STABILITYDERIVATIVES_CALC MATLAB code for stabilityDerivatives_calc.fig
%      STABILITYDERIVATIVES_CALC, by itself, creates a new STABILITYDERIVATIVES_CALC or raises the existing
%      singleton*.
%
%      H = STABILITYDERIVATIVES_CALC returns the handle to a new STABILITYDERIVATIVES_CALC or the handle to
%      the existing singleton*.
%
%      STABILITYDERIVATIVES_CALC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STABILITYDERIVATIVES_CALC.M with the given input arguments.
%
%      STABILITYDERIVATIVES_CALC('Property','Value',...) creates a new STABILITYDERIVATIVES_CALC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stabilityDerivatives_calc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stabilityDerivatives_calc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stabilityDerivatives_calc

% Last Modified by GUIDE v2.5 05-Oct-2015 16:14:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stabilityDerivatives_calc_OpeningFcn, ...
                   'gui_OutputFcn',  @stabilityDerivatives_calc_OutputFcn, ...
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


% --- Executes just before stabilityDerivatives_calc is made visible.
function stabilityDerivatives_calc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stabilityDerivatives_calc (see VARARGIN)

% Choose default command line output for stabilityDerivatives_calc
% handles.modelo = varargin{1};
% handles.output = modelo;
handles.output = 0;

set(handles.longStabDer, 'data', zeros(4,3));
set(handles.longControlDer, 'data', zeros(2,3));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stabilityDerivatives_calc wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stabilityDerivatives_calc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.alpha = get(hObject,'String');
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.W_W0 = get(hObject,'String');
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_button.
function calc_button_Callback(hObject, eventdata, handles)
% hObject    handle to calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
