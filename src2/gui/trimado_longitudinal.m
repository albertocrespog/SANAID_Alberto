function varargout = trimado_longitudinal(varargin)
% TRIMADO_LONGITUDINAL MATLAB code for trimado_longitudinal.fig
%      TRIMADO_LONGITUDINAL, by itself, creates a new TRIMADO_LONGITUDINAL or raises the existing
%      singleton*.
%
%      H = TRIMADO_LONGITUDINAL returns the handle to a new TRIMADO_LONGITUDINAL or the handle to
%      the existing singleton*.
%
%      TRIMADO_LONGITUDINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIMADO_LONGITUDINAL.M with the given input arguments.
%
%      TRIMADO_LONGITUDINAL('Property','Value',...) creates a new TRIMADO_LONGITUDINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trimado_longitudinal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trimado_longitudinal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trimado_longitudinal

% Last Modified by GUIDE v2.5 27-Mar-2016 17:12:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trimado_longitudinal_OpeningFcn, ...
                   'gui_OutputFcn',  @trimado_longitudinal_OutputFcn, ...
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


% --- Executes just before trimado_longitudinal is made visible.
function trimado_longitudinal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trimado_longitudinal (see VARARGIN)
handles.modelo  = varargin{1};
handles.dir     = varargin{2};

[check,log] = checkModel(handles.modelo);
if check
    if strcmp(handles.modelo.conf,'convencional_canard') 
        set(handles.movilSurf_txt,'Visible','on');
        set(handles.movilSurf_val,'Visible','on');
        handles.fixSurf = 'canard';
    else
        set(handles.movilSurf_txt,'Visible','off');
        set(handles.movilSurf_val,'Visible','off');
        handles.fixSurf = [];
        if strcmp(handles.modelo.conf,'convencional') || strcmp(handles.modelo.conf,'flyWing')
            set(handles.CLdef_txt,'String','CLde (1/rad)');
            set(handles.def_txt,'String','de (º)');
            set(handles.CMdef_txt,'String','CMde (1/rad)');
            set(handles.deflec_radioButton,'String','de (º)');
            set(handles.CMdef_radioButton,'String','CMdc (1/rad)');
        else
            set(handles.def_txt,'String','dc (º)');
            set(handles.CLdef_txt,'String','CLdc (1/rad)');
            set(handles.CMdef_txt,'String','CMde (1/rad)');
            set(handles.deflec_radioButton,'String','dc (º)');
            set(handles.CMdef_radioButton,'String','CMdc (1/rad)');
        end
    end
    handles.resolMode   = 'basic';
    handles.XCG_evol    = [];
    handles.W_W0_evol   = [];

    set(handles.surfDef_val,'Enable','inactive');
    set(handles.alpha_val,'Enable','inactive');
    set(handles.graphSave_button,'Enable','inactive');

    cla(handles.axes1);
    grid on;
    xlabel('W/W_0');
    ylabel('X_{CG}/l_{fus}');

    cla(handles.axes2);
    xlabel('W/W_0');
    grid on;

% Choose default command line output for trimado_longitudinal
%handles.output = hObject;

% Update handles structure
    guidata(hObject, handles);

% UIWAIT makes trimado_longitudinal wait for user response (see UIRESUME)
    uiwait(handles.figure1);
else
    mes = 'Next fields need to be specified to perform Longitudinal Analysis:';
    insufficientData_dialog(handles.dir,mes,log);
    delete(hObject);
end


% --- Outputs from this function are returned to the command line.
function varargout = trimado_longitudinal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

delete(hObject)



function WW0_val_Callback(hObject, eventdata, handles)
% hObject    handle to WW0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WW0_val as text
%        str2double(get(hObject,'String')) returns contents of WW0_val as a double
handles.W_W0_evol = str2num(get(hObject,'String'));
if ~isempty(handles.XCG_evol) && length(handles.XCG_evol) == length(handles.W_W0_evol) && length(handles.W_W0_evol) > 1
    axes(handles.axes1);
    cla(handles.axes1);
    hold on; grid on;
    plot(handles.W_W0_evol,handles.XCG_evol,'b','linewidth',1.5);
    xlabel('W/W_0');
    ylabel('X_{CG}/l_{fus}');
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function WW0_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WW0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XCG_val_Callback(hObject, eventdata, handles)
% hObject    handle to XCG_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XCG_val as text
%        str2double(get(hObject,'String')) returns contents of XCG_val as a double
handles.XCG_evol = str2num(get(hObject,'String'));
if ~isempty(handles.W_W0_evol) && length(handles.XCG_evol) == length(handles.W_W0_evol) && length(handles.W_W0_evol) > 1
    axes(handles.axes1);
    cla(handles.axes1);
    hold on; grid on;
    plot(handles.W_W0_evol,handles.XCG_evol,'b','linewidth',1.5);
    xlabel('W/W_0');
    ylabel('X_{CG}/l_{fus}');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function XCG_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XCG_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[matFileName,matPathName,~] = uigetfile('*.mat','Select the .mat file',handles.dir.main);
if ~isequal(matPathName, 0) && ~isequal(matFileName, 0)
    handles.datos = load([matPathName, matFileName]);
    obj_name = fieldnames(handles.datos);
    handles.datos = handles.datos.(obj_name{1});
    handles.datosSize = size(handles.datos);
    
    if handles.datosSize(1) == 2
        handles.W_W0_evol   = handles.datos(1,:);
        set(handles.WW0_val,'String',num2str(handles.W_W0_evol,2));
        handles.XCG_evol    = handles.datos(2,:);
        set(handles.XCG_val,'String',num2str(handles.XCG_evol,2));
        
        cla(handles.axes1);
        hold on; grid on;
        plot(handles.W_W0_evol,handles.XCG_evol,'linewidth',1.5);
        xlabel('W/W_0');
        ylabel('X_{CG}/l_{fus}');
    end
end
guidata(hObject, handles);


% --- Executes on selection change in movilSurf_val.
function movilSurf_val_Callback(hObject, eventdata, handles)
% hObject    handle to movilSurf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns movilSurf_val contents as cell array
%        contents{get(hObject,'Value')} returns selected item from movilSurf_val
if get(hObject,'Value') == 1
    handles.fixSurf = 'canard';
    set(handles.def_txt,'String','de (º)');
    set(handles.deflec_radioButton,'String','de (º)');
elseif get(hObject,'Value') == 2
    handles.fixSurf = 'horizontal';
    set(handles.def_txt,'String','dc (º)');
    set(handles.deflec_radioButton,'String','dc (º)');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function movilSurf_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movilSurf_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trim_button.
function trim_button_Callback(hObject, eventdata, handles)
% hObject    handle to trim_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
cd(handles.dir.code);
handles.results = trim_longitudinal(handles.modelo, handles.XCG_evol, handles.W_W0_evol, handles.resolMode, handles.fixSurf); 
cd(handles.dir.gui);
handles.alphaTrim   = handles.results(1,:);
handles.deflecTrim  = handles.results(2,:);
handles.CDiTrim     = handles.results(3,:);
handles.CL0Trim     = handles.results(4,:);
handles.CM0Trim     = handles.results(5,:);
handles.CLaTrim     = handles.results(6,:);
handles.CLdefTrim   = handles.results(7,:);
handles.CMaTrim     = handles.results(8,:);
handles.CMdefTrim   = handles.results(9,:);

if length(handles.XCG_evol) == 1
    set(handles.alpha_val,'String',num2str(handles.alphaTrim*180/pi,3));
    set(handles.surfDef_val,'String',num2str(handles.deflecTrim*180/pi,3));
    set(handles.CLa_val,'String', num2str(handles.CLaTrim,3));
    set(handles.CLdef_val,'String', num2str(handles.CLdefTrim,3));
    set(handles.CLa_val,'String', num2str(handles.CLaTrim,3));
    set(handles.CL0_val,'String', num2str(handles.CL0Trim,3));
    set(handles.CM0_val,'String', num2str(handles.CM0Trim,3));
    set(handles.CMa_val,'String', num2str(handles.CMaTrim,3));
    set(handles.CMdef_val,'String', num2str(handles.CMdefTrim,3));
else
    set(handles.alpha_val,'String',num2str(handles.alphaTrim(1)*180/pi,3));
    set(handles.surfDef_val,'String',num2str(handles.deflecTrim(1)*180/pi,3));
    set(handles.CLa_val,'String', num2str(handles.CLaTrim(1),3));
    set(handles.CLdef_val,'String', num2str(handles.CLdefTrim(1),3));
    set(handles.CLa_val,'String', num2str(handles.CLaTrim(1),3));
    set(handles.CL0_val,'String', num2str(handles.CL0Trim(1),3));
    set(handles.CM0_val,'String', num2str(handles.CM0Trim(1),3));
    set(handles.CMa_val,'String', num2str(handles.CMaTrim(1),3));
    set(handles.CMdef_val,'String', num2str(handles.CMdefTrim(1),3));
    
    axes(handles.axes2);
    cla;
    hold on;
    grid on;
    alpha   = plot(handles.W_W0_evol,handles.alphaTrim*180/pi,'r','linewidth',2);
    handles.alphaGraph = hgtransform('Parent',handles.axes2);
    set(alpha,'Parent',handles.alphaGraph);
    
    def     = plot(handles.W_W0_evol,handles.deflecTrim*180/pi,'b','linewidth',2);
    handles.deflecGraph = hgtransform('Parent',handles.axes2);
    set(def,'Parent',handles.deflecGraph);
    
    CDi     = plot(handles.W_W0_evol,handles.CDiTrim,'g','linewidth',2);
    handles.CDiGraph = hgtransform('Parent',handles.axes2);
    set(CDi,'Parent',handles.CDiGraph);
    
    CMa     = plot(handles.W_W0_evol,handles.CMaTrim,'m','linewidth',2);
    handles.CMaGraph = hgtransform('Parent',handles.axes2);
    set(CMa,'Parent',handles.CMaGraph);
    
    CMdef     = plot(handles.W_W0_evol,handles.CMdefTrim,'c','linewidth',2);
    handles.CMdefGraph = hgtransform('Parent',handles.axes2);
    set(CMdef,'Parent',handles.CMdefGraph);
    
    CM0     = plot(handles.W_W0_evol,handles.CM0Trim,'k','linewidth',2);
    handles.CM0Graph = hgtransform('Parent',handles.axes2);
    set(CM0,'Parent',handles.CM0Graph);
    
    
    if get(handles.deflec_radioButton,'Value')
        set(handles.deflecGraph,'Visible','on');
    else
        set(handles.deflecGraph,'Visible','off');
    end
    
    if get(handles.alpha_radioButton,'Value')
        set(handles.alphaGraph,'Visible','on');
    else
        set(handles.alphaGraph,'Visible','off');
    end
    
    if get(handles.CDi_radioButton,'Value')
        set(handles.CDiGraph,'Visible','on');
    else
        set(handles.CDiGraph,'Visible','off');
    end
    
    if get(handles.CMa_radioButton,'Value')
        set(handles.CMaGraph,'Visible','on');
    else
        set(handles.CMaGraph,'Visible','off');
    end
    
    if get(handles.CMdef_radioButton,'Value')
        set(handles.CMdefGraph,'Visible','on');
    else
        set(handles.CMdefGraph,'Visible','off');
    end
    
    if get(handles.CM0_radioButton,'Value')
        set(handles.CM0Graph,'Visible','on');
    else
        set(handles.CM0Graph,'Visible','off');
    end
end
    
guidata(hObject, handles);


function alpha_val_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_val as text
%        str2double(get(hObject,'String')) returns contents of alpha_val as a double


% --- Executes during object creation, after setting all properties.
function alpha_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function surfDef_val_Callback(hObject, eventdata, handles)
% hObject    handle to surfDef_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of surfDef_val as text
%        str2double(get(hObject,'String')) returns contents of surfDef_val as a double


% --- Executes during object creation, after setting all properties.
function surfDef_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to surfDef_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveResults_button.
function saveResults_button_Callback(hObject, eventdata, handles)
% hObject    handle to saveResults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'results') && ~isempty(handles.results)
    [writefname, writepname] = uiputfile('*.mat','Save results as');
    writepfname = fullfile(writepname, writefname);
    results = handles.results;
    save(writepfname,'results');
    %uisave('handles.results','Trim_Results');
end


% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);
guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
guidata(hObject, handles);


% --- Executes on button press in graphSave_button.
function graphSave_button_Callback(hObject, eventdata, handles)
% hObject    handle to graphSave_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in alpha_radioButton.
function alpha_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha_radioButton
if ~isempty(handles.results)
    if get(hObject,'Value')
        set(handles.alphaGraph,'Visible','on');
    else
        set(handles.alphaGraph,'Visible','off');
    end
end
guidata(hObject, handles);


% --- Executes on button press in deflec_radioButton.
function deflec_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to deflec_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deflec_radioButton
if ~isempty(handles.results)
    if get(hObject,'Value')
        set(handles.deflecGraph,'Visible','on');
    else
        set(handles.deflecGraph,'Visible','off');
    end
end
guidata(hObject, handles);


% --- Executes on button press in CDi_radioButton.
function CDi_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to CDi_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CDi_radioButton

if ~isempty(handles.results)
    if get(hObject,'Value')
        set(handles.CDiGraph,'Visible','on');
    else
        set(handles.CDiGraph,'Visible','off');
    end
end
guidata(hObject, handles)


% --- Executes on button press in basic_radioButton.
function basic_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to basic_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of basic_radioButton
handles.resolMode = 'basic';
guidata(hObject, handles);


% --- Executes on button press in adv_radioButton.
function adv_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to adv_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adv_radioButton
handles.resolMode = 'advan';
% Not available yet
set(handles.basic_radioButton,'Value',1);
handles.resolMode = 'basic';
guidata(hObject, handles);


% --- Executes on button press in CMa_radioButton.
function CMa_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to CMa_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CMa_radioButton
if ~isempty(handles.results)
    if get(hObject,'Value')
        set(handles.CMaGraph,'Visible','on');
    else
        set(handles.CMaGraph,'Visible','off');
    end
end

guidata(hObject, handles)


% --- Executes on button press in CMdef_radioButton.
function CMdef_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to CMdef_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CMdef_radioButton
if ~isempty(handles.results)
    if get(hObject,'Value')
        set(handles.CMdefGraph,'Visible','on');
    else
        set(handles.CMdefGraph,'Visible','off');
    end
end

guidata(hObject, handles)



function CLa_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLa_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLa_val as text
%        str2double(get(hObject,'String')) returns contents of CLa_val as a double


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



function CLdef_val_Callback(hObject, eventdata, handles)
% hObject    handle to CLdef_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLdef_val as text
%        str2double(get(hObject,'String')) returns contents of CLdef_val as a double


% --- Executes during object creation, after setting all properties.
function CLdef_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLdef_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CMa_val_Callback(hObject, eventdata, handles)
% hObject    handle to CMa_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CMa_val as text
%        str2double(get(hObject,'String')) returns contents of CMa_val as a double


% --- Executes during object creation, after setting all properties.
function CMa_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CMa_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CMdef_val_Callback(hObject, eventdata, handles)
% hObject    handle to CMdef_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CMdef_val as text
%        str2double(get(hObject,'String')) returns contents of CMdef_val as a double


% --- Executes during object creation, after setting all properties.
function CMdef_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CMdef_val (see GCBO)
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


function CM0_val_Callback(hObject, eventdata, handles)
% hObject    handle to CM0_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CM0_val as text
%        str2double(get(hObject,'String')) returns contents of CM0_val as a double


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


% --- Executes on button press in CM0_radioButton.
function CM0_radioButton_Callback(hObject, eventdata, handles)
% hObject    handle to CM0_radioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CM0_radioButton

if ~isempty(handles.results)
    if get(hObject,'Value')
        set(handles.CM0Graph,'Visible','on');
    else
        set(handles.CM0Graph,'Visible','off');
    end
end

guidata(hObject, handles)


function [check,log] = checkModel(model)
log = '';
check = 1;

% GENERAL DATA
if isempty(model.general.Sref) || isnan(model.general.Sref)
    log = [log, 'Sref '];
    check = 0;
end

if isempty(model.general.mtow) || isnan(model.general.mtow)
    log = [log,'MTOW '];
    check = 0;
end

if isempty(model.general.h) || isnan(model.general.h)
    log = [log,'h '];
    check = 0;
end

if isempty(model.general.Minf) || isnan(model.general.Minf)
    log = [log,'Minf '];
    check = 0;
end

if isempty(model.general.L) || isnan(model.general.L)
    log = [log,'L '];
    check = 0;
end

% WING DATA
if isempty(model.ala.S) || isnan(model.ala.S)
    log = [log,'S_w '];
    check = 0;
end

if isempty(model.ala.CLa) || isnan(model.ala.CLa)
    log = [log,'CLa_w '];
    check = 0;
end

if isempty(model.ala.eta) || isnan(model.ala.eta)
    log = [log,'eta_w '];
    check = 0;
end

if isempty(model.ala.CL0) || isnan(model.ala.CL0)
    log = [log,'CL0_w '];
    check = 0;
end

if isempty(model.ala.CM0) || isnan(model.ala.CM0)
    log = [log,'CM0_w '];
    check = 0;
end


if isempty(model.ala.AR) || isnan(model.ala.AR)
    log = [log,'AR_w '];
    check = 0;
end

if isempty(model.ala.b) || isnan(model.ala.b)
    log = [log,'b_w '];
    check = 0;
end

if isempty(model.ala.TR) || isnan(model.ala.TR)
    log = [log,'TR_w '];
    check = 0;
end

if isempty(model.ala.Xca) || isnan(model.ala.Xca)
    log = [log,'Xca_w '];
    check = 0;
end

if isempty(model.ala.i) || isnan(model.ala.i)
    log = [log,'i_w '];
    check = 0;
end

if isempty(model.ala.cr) || isnan(model.ala.cr)
    log = [log,'c_root_w '];
    check = 0;
end

if isempty(model.ala.LAMc4) || isnan(model.ala.LAMc4)
    log = [log,'LAMc4_w '];
    check = 0;
end

if isempty(model.ala.Zca) || isnan(model.ala.Zca)
    log = [log,'Zca_w '];
    check = 0;
end

switch model.conf
    case 'canard'
        % CANARD DATA
        if isempty(model.canard.S) || isnan(model.canard.S)
            log = [log,'S_c '];
            check = 0;
        end

        if isempty(model.canard.CLa) || isnan(model.canard.CLa)
            log = [log,'CLa_c '];
            check = 0;
        end

        if isempty(model.canard.eta) || isnan(model.canard.eta)
            log = [log,'eta_c '];
            check = 0;
        end

        if isempty(model.canard.CL0) || isnan(model.canard.CL0)
            log = [log,'CL0_c '];
            check = 0;
        end

        if isempty(model.canard.CM0) || isnan(model.canard.CM0)
            log = [log,'CM0_c '];
            check = 0;
        end

        if isempty(model.canard.AR) || isnan(model.canard.AR)
            log = [log,'AR_c '];
            check = 0;
        end

        if isempty(model.canard.b) || isnan(model.canard.b)
            log = [log,'b_c '];
            check = 0;
        end

        if isempty(model.canard.TR) || isnan(model.canard.TR)
            log = [log,'TR_c '];
            check = 0;
        end

        if isempty(model.canard.Xca) || isnan(model.canard.Xca)
            log = [log,'Xca_c '];
            check = 0;
        end

        if isempty(model.canard.i) || isnan(model.canard.i) 
            log = [log,'i_c '];
            check = 0;
        end

        if isempty(model.canard.y1_b2) || isnan(model.canard.y1_b2)
            log = [log,'y1_b2_c '];
            check = 0;
        end

        if isempty(model.canard.y0_b2) || isnan(model.canard.y0_b2)
            log = [log,'y0_b2_c '];
            check = 0;
        end

        if isempty(model.canard.cm_c) || isnan(model.canard.cm_c)
            log = [log,'cm_c_c '];
            check = 0;
        end
        
        if isempty(model.canard.t_c) || isnan(model.canard.t_c)
            log = [log,'t_c_c '];
            check = 0;
        end

    case 'convencional'
        % HORIZONTAL DATA
        if isempty(model.horizontal.S) || isnan(model.horizontal.S)
            log = [log,'S_h '];
            check = 0;
        end

        if isempty(model.horizontal.CLa) || isnan(model.horizontal.CLa)
            log = [log,'CLa_h '];
            check = 0;
        end

        if isempty(model.horizontal.eta) || isnan(model.horizontal.eta)
            log = [log,'eta_h '];
            check = 0;
        end

        if isempty(model.horizontal.CL0) || isnan(model.horizontal.CL0)
            log = [log,'CL0_h '];
            check = 0;
        end

        if isempty(model.horizontal.CM0) || isnan(model.horizontal.CM0)
            log = [log,'CM0_h '];
            check = 0;
        end

        if isempty(model.horizontal.AR) || isnan(model.horizontal.AR)
            log = [log,'AR_h '];
            check = 0;
        end

        if isempty(model.horizontal.b) || isnan(model.horizontal.b)
            log = [log,'b_h '];
            check = 0;
        end

        if isempty(model.horizontal.TR) || isnan(model.horizontal.TR)
            log = [log,'TR_h '];
            check = 0;
        end

        if isempty(model.horizontal.Xca) || isnan(model.horizontal.Xca)
            log = [log,'Xca_h '];
            check = 0;
        end

        if isempty(model.horizontal.i) || isnan(model.horizontal.i)
            log = [log,'i_h '];
            check = 0;
        end

        if isempty(model.horizontal.y1_b2) || isnan(model.horizontal.y1_b2)
            log = [log,'y1_b2_h '];
            check = 0;
        end

        if isempty(model.horizontal.y0_b2) || isnan(model.horizontal.y0_b2)
            log = [log,'y0_b2_h '];
            check = 0;
        end

        if isempty(model.horizontal.cm_c) || isnan(model.horizontal.cm_c)
            log = [log,'cm_c_h '];
            check = 0;
        end
        
        if isempty(model.horizontal.Zca) || isnan(model.horizontal.Zca)
            log = [log,'Zca_h '];
            check = 0;
        end
        
        if isempty(model.horizontal.t_c) || isnan(model.horizontal.t_c)
            log = [log,'t_c_h '];
            check = 0;
        end
        
    case 'convencional_canard'
         % CANARD DATA
        if isempty(model.canard.S) || isnan(model.canard.S)
            log = [log,'S_c '];
            check = 0;
        end

        if isempty(model.canard.CLa) || isnan(model.canard.CLa)
            log = [log,'CLa_c '];
            check = 0;
        end

        if isempty(model.canard.eta) || isnan(model.canard.eta)
            log = [log,'eta_c '];
            check = 0;
        end

        if isempty(model.canard.CL0) || isnan(model.canard.CL0)
            log = [log,'CL0_c '];
            check = 0;
        end

        if isempty(model.canard.CM0) || isnan(model.canard.CM0)
            log = [log,'CM0_c '];
            check = 0;
        end

        if isempty(model.canard.AR) || isnan(model.canard.AR)
            log = [log,'AR_c '];
            check = 0;
        end

        if isempty(model.canard.b) || isnan(model.canard.b)
            log = [log,'b_c '];
            check = 0;
        end

        if isempty(model.canard.TR) || isnan(model.canard.TR)
            log = [log,'TR_c '];
            check = 0;
        end

        if isempty(model.canard.Xca) || isnan(model.canard.Xca)
            log = [log,'Xca_c '];
            check = 0;
        end

        if isempty(model.canard.i) || isnan(model.canard.i) 
            log = [log,'i_c '];
            check = 0;
        end

        if isempty(model.canard.y1_b2) || isnan(model.canard.y1_b2)
            log = [log,'y1_b2_c '];
            check = 0;
        end

        if isempty(model.canard.y0_b2) || isnan(model.canard.y0_b2)
            log = [log,'y0_b2_c '];
            check = 0;
        end

        if isempty(model.canard.cm_c) || isnan(model.canard.cm_c)
            log = [log,'cm_c_c '];
            check = 0;
        end
        
        if isempty(model.canard.t_c) || isnan(model.canard.t_c)
            log = [log,'t_c_c '];
            check = 0;
        end

        % HORIZONTAL DATA
        if isempty(model.horizontal.S) || isnan(model.horizontal.S)
            log = [log,'S_h '];
            check = 0;
        end

        if isempty(model.horizontal.CLa) || isnan(model.horizontal.CLa)
            log = [log,'CLa_h '];
            check = 0;
        end

        if isempty(model.horizontal.eta) || isnan(model.horizontal.eta)
            log = [log,'eta_h '];
            check = 0;
        end

        if isempty(model.horizontal.CL0) || isnan(model.horizontal.CL0)
            log = [log,'CL0_h '];
            check = 0;
        end

        if isempty(model.horizontal.CM0) || isnan(model.horizontal.CM0)
            log = [log,'CM0_h '];
            check = 0;
        end

        if isempty(model.horizontal.AR) || isnan(model.horizontal.AR)
            log = [log,'AR_h '];
            check = 0;
        end

        if isempty(model.horizontal.b) || isnan(model.horizontal.b)
            log = [log,'b_h '];
            check = 0;
        end

        if isempty(model.horizontal.TR) || isnan(model.horizontal.TR)
            log = [log,'TR_h '];
            check = 0;
        end

        if isempty(model.horizontal.Xca) || isnan(model.horizontal.Xca)
            log = [log,'Xca_h '];
            check = 0;
        end

        if isempty(model.horizontal.i) || isnan(model.horizontal.i)
            log = [log,'i_h '];
            check = 0;
        end

        if isempty(model.horizontal.y1_b2) || isnan(model.horizontal.y1_b2)
            log = [log,'y1_b2_h '];
            check = 0;
        end

        if isempty(model.horizontal.y0_b2) || isnan(model.horizontal.y0_b2)
            log = [log,'y0_b2_h '];
            check = 0;
        end

        if isempty(model.horizontal.cm_c) || isnan(model.horizontal.cm_c)
            log = [log,'cm_c_h '];
            check = 0;
        end
        
        if isempty(model.horizontal.Zca) || isnan(model.horizontal.Zca)
            log = [log,'Zca_h '];
            check = 0;
        end
        
        if isempty(model.horizontal.t_c) || isnan(model.horizontal.t_c)
            log = [log,'t_c_h '];
            check = 0;
        end
       
end
