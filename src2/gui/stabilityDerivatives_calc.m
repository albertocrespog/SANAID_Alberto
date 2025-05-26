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

% Last Modified by GUIDE v2.5 06-Oct-2015 18:35:51

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
handles.model = varargin{1};
handles.output = handles.model;

switch(handles.model.confVert)
    case 'no_vert'
        vertField = fieldnames(handles.model.vertical);
        nField = length(vertField);
        for k = 1:nField
            handles.model.vertical.(vertField{k}) = 0;
        end
    case 'convencional' 
        handles.vertFact = 1;
        
    case 'twin_vertical'
        handles.vertFact = 2;
        handles.model.vertical.S = handles.model.vertical.S*handles.vertFact;
end

handles.alpha = 0;

set(handles.alpha_edit, 'String', num2str(handles.alpha));
set(handles.wW0_edit, 'String', num2str(handles.model.general.w_w0));
set(handles.h_edit, 'String', num2str(handles.model.general.h));
set(handles.M_edit, 'String', num2str(handles.model.general.Minf));

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
delete(handles.figure1);



function alpha_edit_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_edit as text
%        str2double(get(hObject,'String')) returns contents of alpha_edit as a double
handles.alpha = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function alpha_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wW0_edit_Callback(hObject, eventdata, handles)
% hObject    handle to wW0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wW0_edit as text
%        str2double(get(hObject,'String')) returns contents of wW0_edit as a double
handles.model.general.w_w0 = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function wW0_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wW0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in exit_button.
function calc_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.model = calc_derivadasEstabilidad(handles.model, handles.alpha);
deriv = handles.model.derivadas;

long_mat        =   [deriv.CL_a,        deriv.CD_a,         deriv.CM_a; 
                    deriv.CL_u,         deriv.CD_u,         deriv.CM_u;
                    deriv.CL_q,         deriv.CD_q,         deriv.CM_q; 
                    deriv.CL_alphaDot,  deriv.CD_alphaDot,  deriv.CM_alphaDot];


longCont_mat    = [deriv.CL_de,deriv.CD_de,deriv.CM_de; deriv.CL_dc,deriv.CD_dc,deriv.CM_dc];


lat_mat         = [deriv.Cy_beta,deriv.Cl_beta,deriv.Cn_beta; deriv.Cy_p,deriv.Cl_p,deriv.Cn_p;deriv.Cy_r,deriv.Cl_r,deriv.Cn_r; deriv.Cy_betaDot,deriv.Cl_betaDot,deriv.Cn_betaDot];


latCont_mat     = [deriv.Cy_dr,deriv.Cl_dr,deriv.Cn_dr;deriv.Cy_da,deriv.Cl_da,deriv.Cn_da];


prop_mat        = [deriv.Cy_Tbeta,deriv.CM_Ta,deriv.CM_Tu,deriv.CM_T1,deriv.CT_x1,deriv.CT_xu,deriv.CT_xa,deriv.Cn_Tbeta];
% long_mat = rand(4,3);
% longCont_mat = rand(2,3);
% lat_mat = rand(4,3);
% latCont_mat = rand(2,3);
% prop_mat    = rand(1,8);

set(handles.longStabDer, 'Data', long_mat);
set(handles.longControlDer, 'Data', longCont_mat);
set(handles.latStabDer, 'Data', lat_mat);
set(handles.latControlDer, 'Data', latCont_mat);
set(handles.propStabDer, 'Data', prop_mat);

guidata(hObject, handles);



function h_edit_Callback(hObject, eventdata, handles)
% hObject    handle to h_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h_edit as text
%        str2double(get(hObject,'String')) returns contents of h_edit as a double
handles.model.general.h = str2double(get(hObject, 'String'));

[handles.model.general.Tinf, handles.model.general.rhoinf, handles.model.general.pinf, handles.model.general.ainf] = atmos_inter(handles.model.general.h);


handles.model.general.Vinf = handles.model.general.Minf*handles.model.general.ainf;
handles.model.general.qinf = handles.model.general.rhoinf*(handles.model.general.Vinf)^2;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function h_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function M_edit_Callback(hObject, eventdata, handles)
% hObject    handle to M_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_edit as text
%        str2double(get(hObject,'String')) returns contents of M_edit as a double
handles.model.general.Minf = str2double(get(hObject, 'String'));
handles.model.general.Vinf = handles.model.general.Minf*handles.model.general.ainf;
handles.model.general.qinf = handles.model.general.rhoinf*(handles.model.general.Vinf)^2;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function M_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = handles.model;
guidata(hObject, handles);
