function varargout = caracteristicas_avanz(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @caracteristicas_avanz_OpeningFcn, ...
                   'gui_OutputFcn',  @caracteristicas_avanz_OutputFcn, ...
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

function caracteristicas_avanz_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);


scrsz = get(0, 'ScreenSize'); 
pos_act=get(gcf,'Position');   
xr=scrsz(3) - pos_act(3); 
xp=round(xr/2);   
yr=scrsz(4) - pos_act(4); 
yp=round(yr/2);   
set(handles.figure1,'Position',[xp yp pos_act(3) pos_act(4)]);

global modificar_modelo archivo_modelo nombre_modelo boton_previo

if modificar_modelo == 1,
    set(handles.edit31,'String',nombre_modelo,'Enable','off');
    load(archivo_modelo);
    if propul(1) == 1,
    propul(3) = propul(3) * 1/4.448221615255; %[N]
    propul(4) = propul(4) * 1/(2.832546065 * 10^-5); %[kg/(N*s)]
    else
    propul(3) = propul(3) * 1/745.699872; %[W]
    propul(4) = propul(4) * 1/(1.68965941 * 10^-7); %[kg/(W*s)]
    end
    
    handles.propul(1) = propul(1); % Type of engine
    handles.propul(2) = propul(2); % Number of engines
    handles.propul(3) = propul(3); % Thrust (lbf) or Power (shp) per engine
    handles.propul(4) = propul(4); % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
    handles.propul(5) = propul(5); % Normativa
    handles.propul(6) = propul(6); % Prop efficiency
    handles.propul(7) = propul(7); % By-pass
    
    handles.aerodinamica(1) = aerodinamica(1); % Superficie alar
    handles.aerodinamica(2) = aerodinamica(2); % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
    handles.aerodinamica(3) = aerodinamica(3); % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
    handles.aerodinamica(4) = aerodinamica(4); % Coeficiente de resistencia inducida K2: CD = CD0 + K1*CL^2 - K2*CL
    handles.aerodinamica(5) = aerodinamica(5); % Coeficiente de sustentación máxima en limpio
    
    handles.aero_despegue(1) = aero_despegue(1); % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
    handles.aero_despegue(2) = aero_despegue(2); % Coeficiente de resistencia inducida K1: CD = CD0 + K1*CL^2 - K2*CL
    handles.aero_despegue(3) = aero_despegue(3); % Coeficiente de sustentación para ángulo de ataque nulo
    handles.aero_despegue(4) = aero_despegue(4); % Coeficiente de sustentación máxima en sucio Take Off
    
    handles.aero_aterrizaje(1) = aero_aterrizaje(1); % Coeficiente de resistencia parasitaria CD0: CD = CD0 + K1*CL^2 - K2*CL
    handles.aero_aterrizaje(2) = aero_aterrizaje(2); % Coeficiente de resistencia inducida K2: CD = CD0 + K1*CL^2 - K2*CL
    handles.aero_aterrizaje(3) = aero_aterrizaje(3); % Coeficiente de sustentación para ángulo de ataque nulo
    handles.aero_aterrizaje(4) = aero_aterrizaje(4); % Coeficiente de sustentación máxima en sucio Landing
    
    handles.pesos(1) = pesos(1); % Peso en vacío
    handles.pesos(2) = pesos(2); % Carga de pago
    handles.pesos(3) = pesos(3); % Peso de tripulación
    handles.pesos(4) = pesos(4); % % Fuel restante al final
  

else
    set(handles.edit31,'String','','Enable','on');
    handles.aerodinamica = [0 0 0 0 0];
    handles.aero_despegue = [0 0 0 0];
    handles.aero_aterrizaje = [0 0 0 0];
    handles.propul = [1 0 0 0 1 0 1];
    handles.pesos = [0 0 0 0];
end
set(handles.popupmenu1,'Value',handles.propul(1));
set(handles.popupmenu2,'Value',handles.propul(5));
set(handles.popupmenu3,'Value',handles.propul(7));


posicion = get(handles.pushbutton1,'Position');
if abs(posicion(4) - 0.8)<10^-1,
posicion(4) = posicion(4)*1.2;
set(handles.pushbutton1,'BackgroundColor',[1 1 1]);
set(handles.pushbutton1,'Position',posicion);
set(handles.pushbutton2,'Position',[3.6 13.933 2.1 0.8]);
set(handles.pushbutton2,'BackgroundColor',[0.941 0.941 0.941]);
set(handles.pushbutton3,'Position',[5.65 13.933 2.1 0.8]);
set(handles.pushbutton3,'BackgroundColor',[0.941 0.941 0.941]);
set(handles.pushbutton1,'Visible','off');
set(handles.pushbutton2,'Visible','off');
set(handles.pushbutton3,'Visible','off');
drawnow;
set(handles.pushbutton1,'Visible','on');
set(handles.pushbutton2,'Visible','on');
set(handles.pushbutton3,'Visible','on');
end

set(handles.edit1,'String',num2str(handles.aero_despegue(1)));    %cd0
set(handles.edit2,'String',num2str(handles.aero_despegue(2)));    %k1
set(handles.edit3,'String',num2str(handles.aero_despegue(3)));    %cl
set(handles.edit4,'String',num2str(handles.aero_despegue(4)));    %clmax

set(handles.edit5,'String',num2str(handles.aero_aterrizaje(1)));    %cd0
set(handles.edit6,'String',num2str(handles.aero_aterrizaje(2)));    %k1
set(handles.edit7,'String',num2str(handles.aero_aterrizaje(3)));    %cl
set(handles.edit8,'String',num2str(handles.aero_aterrizaje(4)));    %clmax

set(handles.edit9,'String',num2str(handles.aerodinamica(1)));    %s
set(handles.edit10,'String',num2str(handles.aerodinamica(2)));    %cd0
set(handles.edit11,'String',num2str(handles.aerodinamica(3)));    %k1
set(handles.edit12,'String',num2str(handles.aerodinamica(4)));    %k2
set(handles.edit13,'String',num2str(handles.aerodinamica(5)));    %clmax

boton_previo = 1;
guidata(hObject,handles);

function varargout = caracteristicas_avanz_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function pushbutton1_Callback(hObject, eventdata, handles)
global archivo_modelo modificar_modelo boton_previo 

color = get(hObject,'BackgroundColor');
if color(1) == 1,return; end;

set(handles.uipanel2,'Visible','Off');
set(handles.uipanel3,'Visible','Off');
set(handles.uipanel4,'Visible','Off');
try 
    drawnow; 
catch
end;


set(handles.uipanel2,'Title','Aerodinamica en el despegue');
set(handles.uipanel4,'Title','Aerodinamica general');

set(handles.text1,'String','Coeficiente de resistencia parasitaria');
set(handles.text2,'String','Coeficiente de resistencia inducida k1');
set(handles.text3,'String','Coeficiente de sustentacion para AoA=0');
set(handles.text4,'String','Coeficiente de sustentacion max');

set(handles.text9,'String','Superficie alar');
set(handles.text10,'String','Coeficiente de resistencia parasitaria');
set(handles.text11,'String','Coeficiente de resistencia inducida k1');
set(handles.text12,'String','Coeficiente de resistencia inducida k2');

set(handles.edit1,'String',num2str(handles.aero_despegue(1)));    %cd0
set(handles.edit2,'String',num2str(handles.aero_despegue(2)));    %k1
set(handles.edit3,'String',num2str(handles.aero_despegue(3)));    %cl
set(handles.edit4,'String',num2str(handles.aero_despegue(4)));    %clmax

set(handles.edit5,'String',num2str(handles.aero_aterrizaje(1)));    %cd0
set(handles.edit6,'String',num2str(handles.aero_aterrizaje(2)));    %k1
set(handles.edit7,'String',num2str(handles.aero_aterrizaje(3)));    %cl
set(handles.edit8,'String',num2str(handles.aero_aterrizaje(4)));    %clmax

set(handles.edit9,'String',num2str(handles.aerodinamica(1)));    %s
set(handles.edit10,'String',num2str(handles.aerodinamica(2)));    %cd0
set(handles.edit11,'String',num2str(handles.aerodinamica(3)));    %k1
set(handles.edit12,'String',num2str(handles.aerodinamica(4)));    %k2
set(handles.edit13,'String',num2str(handles.aerodinamica(5)));    %clmax

set(handles.text13,'Visible','on');
set(handles.edit13,'Visible','on');
set(handles.text12,'Visible','on');
set(handles.edit12,'Visible','on');
set(handles.edit11,'Visible','on');
set(handles.text31,'Visible','off');
set(handles.text32,'Visible','off');
set(handles.text33,'Visible','off');
set(handles.text34,'Visible','off');
set(handles.text35,'Visible','on');

set(handles.edit10,'Enable','on');
set(handles.text10,'Enable','on');
set(handles.text11,'Enable','on');

set(handles.popupmenu1,'Visible','off');
set(handles.popupmenu2,'Visible','off');
set(handles.popupmenu3,'Visible','off');
set(handles.edit1,'Visible','on');
set(handles.edit9,'Visible','on');

set(handles.uipanel2,'Visible','On');
set(handles.uipanel3,'Visible','On');
set(handles.uipanel4,'Visible','On');

posicion = get(hObject,'Position');
if abs(posicion(4) - 0.8)<10^-2,
posicion(4) = posicion(4)*1.2;
set(hObject,'BackgroundColor',[1 1 1]);
set(hObject,'Position',posicion);
set(handles.pushbutton2,'Position',[3.6 13.933 2.1 0.8]);
set(handles.pushbutton2,'BackgroundColor',[0.941 0.941 0.941]);
set(handles.pushbutton3,'Position',[5.65 13.933 2.1 0.8]);
set(handles.pushbutton3,'BackgroundColor',[0.941 0.941 0.941]);
set(hObject,'Visible','off');
set(handles.pushbutton2,'Visible','off');
set(handles.pushbutton3,'Visible','off');
drawnow;
set(hObject,'Visible','on');
set(handles.pushbutton2,'Visible','on');
set(handles.pushbutton3,'Visible','on');
end

guidata(hObject,handles);


function pushbutton2_Callback(hObject, eventdata, handles)

color = get(hObject,'BackgroundColor');
if color(1) == 1,return; end;

if handles.propul(1) == 1,
set(handles.text3,'String','Empuje a nivel del mar (por motor)');
set(handles.text33,'String','lbf');
set(handles.text34,'String','lb/(lbf*h)');
set(handles.popupmenu3,'Enable','on');
set(handles.edit10,'Enable','off');
set(handles.text10,'Enable','off');
set(handles.text11,'Enable','on');
else
set(handles.text3,'String','Potencia a nivel del mar (por motor)');
set(handles.text33,'String','shp');
set(handles.text34,'String','lb/(shp*h)');
set(handles.popupmenu3,'Enable','off');
set(handles.edit10,'Enable','on');
set(handles.text10,'Enable','on');
set(handles.text11,'Enable','off');
end

set(handles.uipanel2,'Visible','Off');
set(handles.uipanel3,'Visible','Off');
set(handles.uipanel4,'Visible','Off');
try 
    drawnow;
catch
end

set(handles.uipanel2,'Title','Propulsion');
set(handles.uipanel4,'Title','Propulsion');

set(handles.text1,'String','Tipo de motor');
set(handles.text2,'String','Numero de motores');
motor = get(handles.popupmenu1,'Value');
if motor == 1,
set(handles.text3,'String','Empuje a nivel del mar (por motor)');
set(handles.text33,'String','lbf');
set(handles.text34,'String','lb/(lbf*h)');
else
set(handles.text3,'String','Potencia a nivel del mar (por motor)');
set(handles.text33,'String','shp');
set(handles.text34,'String','lb/(shp*h)');
end
set(handles.text4,'String','Consumo esp. a nivel del mar');
set(handles.text9,'String','Normativa');
set(handles.text10,'String','Rendimiento de la helice');
set(handles.text11,'String','Relacion de derivacion');

set(handles.popupmenu1,'Visible','on');
set(handles.popupmenu2,'Visible','on');
set(handles.popupmenu3,'Visible','on');
set(handles.popupmenu1,'Value',handles.propul(1));
set(handles.popupmenu2,'Value',handles.propul(5));
set(handles.popupmenu3,'Value',handles.propul(7));
set(handles.edit1,'Visible','off');
set(handles.edit9,'Visible','off');

   
set(handles.edit2,'String',num2str(handles.propul(2)));    
set(handles.edit3,'String',num2str(handles.propul(3)));    
set(handles.edit4,'String',num2str(handles.propul(4)));    

set(handles.edit10,'String',num2str(handles.propul(6)));      

set(handles.text13,'Visible','off');
set(handles.edit13,'Visible','off');
set(handles.text12,'Visible','off');
set(handles.edit12,'Visible','off');
set(handles.edit11,'Visible','off');
set(handles.text31,'Visible','off');
set(handles.text32,'Visible','off');
set(handles.text33,'Visible','on');
set(handles.text34,'Visible','on');
set(handles.text35,'Visible','off');

set(handles.uipanel2,'Visible','On');
set(handles.uipanel4,'Visible','On');

posicion = get(hObject,'Position');
if abs(posicion(4) - 0.8)<10^-2,
posicion(4) = posicion(4)*1.2;
set(hObject,'BackgroundColor',[1 1 1]);
set(hObject,'Position',posicion);
set(handles.pushbutton1,'Position',[1.56 13.933 2.1 0.8]);
set(handles.pushbutton1,'BackgroundColor',[0.941 0.941 0.941]);
set(handles.pushbutton3,'Position',[5.65 13.933 2.1 0.8]);
set(handles.pushbutton3,'BackgroundColor',[0.941 0.941 0.941]);
set(hObject,'Visible','off');
set(handles.pushbutton1,'Visible','off');
set(handles.pushbutton3,'Visible','off');
drawnow;
set(hObject,'Visible','on');
set(handles.pushbutton1,'Visible','on');
set(handles.pushbutton3,'Visible','on');
end

guidata(hObject,handles);

function pushbutton3_Callback(hObject, eventdata, handles)

color = get(hObject,'BackgroundColor');
if color(1) == 1,return; end;


set(handles.uipanel2,'Visible','Off');
set(handles.uipanel3,'Visible','Off');
set(handles.uipanel4,'Visible','Off');
try 
drawnow;
catch
end;

set(handles.uipanel2,'Title','Pesos');

set(handles.text1,'String','Peso en vacio');
set(handles.text2,'String','Carga de pago al inicio');
set(handles.text3,'String','Peso de la tripulacion');
set(handles.text4,'String','% Fuel restante al final');

set(handles.edit1,'String',num2str(handles.pesos(1)));    
set(handles.edit2,'String',num2str(handles.pesos(2)));    
set(handles.edit3,'String',num2str(handles.pesos(3)));    
set(handles.edit4,'String',num2str(handles.pesos(4)));    

set(handles.text31,'Visible','on');
set(handles.text31,'String','kg');
set(handles.text32,'Visible','on');
set(handles.text32,'String','kg');
set(handles.text33,'Visible','on');
set(handles.text33,'String','kg');
set(handles.text34,'Visible','off');

set(handles.popupmenu1,'Visible','off');
set(handles.popupmenu2,'Visible','off');
set(handles.edit1,'Visible','on');
set(handles.edit9,'Visible','on');

set(handles.uipanel2,'Visible','On');

posicion = get(hObject,'Position');
if abs(posicion(4) - 0.8)<10^-2,
posicion(4) = posicion(4)*1.2;
set(hObject,'BackgroundColor',[1 1 1]);
set(hObject,'Position',posicion);
set(handles.pushbutton1,'Position',[1.56 13.933 2.1 0.8]);
set(handles.pushbutton1,'BackgroundColor',[0.941 0.941 0.941]);
set(handles.pushbutton2,'Position',[3.6 13.933 2.1 0.8]);
set(handles.pushbutton2,'BackgroundColor',[0.941 0.941 0.941]);
set(hObject,'Visible','off');
set(handles.pushbutton1,'Visible','off');
set(handles.pushbutton2,'Visible','off');
drawnow;
set(hObject,'Visible','on');
set(handles.pushbutton1,'Visible','on');
set(handles.pushbutton2,'Visible','on');
end

set(handles.uipanel2,'Visible','On');



function edit5_Callback(hObject, eventdata, handles)
valor = str2double(get(hObject,'String'));
if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    handles.aero_aterrizaje(1) = str2double(valor);
    return
end
handles.aero_aterrizaje(1) = valor;
guidata(hObject,handles);


function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
valor = str2double(get(hObject,'String'));
if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    handles.aero_aterrizaje(2) = str2double(valor);
    return
end
handles.aero_aterrizaje(2) = valor;
guidata(hObject,handles);


function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
valor = str2double(get(hObject,'String'));
if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    handles.aero_aterrizaje(3) = str2double(valor);
    return
end
handles.aero_aterrizaje(3) = valor;
guidata(hObject,handles);


function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
valor = str2double(get(hObject,'String'));
if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    handles.aero_aterrizaje(4) = str2double(valor);
    return
end
handles.aero_aterrizaje(4) = valor;
guidata(hObject,handles);


function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aero_despegue(1) = valor;
end
if boton3(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.pesos(1) = valor;
end
guidata(hObject,handles);


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
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aero_despegue(2) = valor;
end
if boton2(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.propul(2) = valor;
end
if boton3(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.pesos(2) = valor;
end
guidata(hObject,handles);


function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aero_despegue(3) = valor;
end
if boton2(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.propul(3) = valor;
end
if boton3(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.pesos(3) = valor;
end
guidata(hObject,handles);
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aero_despegue(4) = valor;
end
if boton2(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.propul(4) = valor;
end
if boton3(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.pesos(4) = valor;
end
guidata(hObject,handles);
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aerodinamica(1) = valor;

guidata(hObject,handles);

function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aerodinamica(2) = valor;
end
if boton2(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.propul(6) = valor;
end
guidata(hObject,handles);

function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aerodinamica(3) = valor;
end
if boton2(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.propul(7) = valor;
end
guidata(hObject,handles);
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
boton1 = get(handles.pushbutton1,'BackgroundColor');
boton2 = get(handles.pushbutton2,'BackgroundColor');
boton3 = get(handles.pushbutton3,'BackgroundColor');
if boton1(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor),
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.aerodinamica(4) = valor;
end
if boton2(1) == 1,
    valor = str2double(get(hObject,'String'));
    if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
    end
handles.propul(8) = valor;
end
guidata(hObject,handles);
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
valor = str2double(get(hObject,'String'));
if isnan(valor) || valor<0,
    msgbox('El dato introducido no es valido o esta fuera de rango','Perfil de modelo','Warn');
    valor = '0';
    set(hObject,'String',valor);
    return
end
handles.aerodinamica(5) = valor;
guidata(hObject,handles);

function pushbutton4_Callback(hObject, eventdata, handles)
global modificar_modelo archivo_modelo nombre_modelo 

handles.dir.main            = pwd;
addpath(handles.dir.main);

handles.dir.gui             = [handles.dir.main,'\gui\'];
addpath(genpath(handles.dir.gui));

handles.dir.data            = [handles.dir.gui,'\data\'];

nombre_modelo = get(handles.edit31,'String');

if isempty(nombre_modelo),
    msgbox('El nombre introducido no es válido','Perfil de modelo','Warn'); 
    return
end

modif = get(handles.edit31,'Enable');
if strcmp(modif,'on'),
load('lista_modelos_avanz.mat');
post = strcmp(nombre_modelo,lista);
if sum(post) > 0
    msgbox('El nombre introducido ya está creado','Perfil de modelo','Warn'); 
    return
end
end

propul(1) = handles.propul(1);
propul(2) = handles.propul(2);
empuje = handles.propul(3);
consumo = handles.propul(4);

if handles.propul(1) == 1,
    propul(3) = empuje * 4.448221615255; %[N]
    propul(4) = consumo * 2.832546065 * 10^-5; %[kg/(N*s)]
else
    propul(3) = empuje * 745.699872; %[W]
    propul(4) = consumo * 1.68965941 * 10^-7; %[kg/(W*s)]
end

propul(5) = handles.propul(5);
propul(6) = handles.propul(6);
propul(7) = handles.propul(7);


aerodinamica(1) = handles.aerodinamica(1);
aerodinamica(2) = handles.aerodinamica(2);
aerodinamica(3) = handles.aerodinamica(3);
aerodinamica(4) = handles.aerodinamica(4);
aerodinamica(5) = handles.aerodinamica(5);

aero_despegue(1) = handles.aero_despegue(1);
aero_despegue(2) = handles.aero_despegue(2);
aero_despegue(3) = handles.aero_despegue(3);
aero_despegue(4) = handles.aero_despegue(4);

aero_aterrizaje(1) = handles.aero_aterrizaje(1);
aero_aterrizaje(2) = handles.aero_aterrizaje(2);
aero_aterrizaje(3) = handles.aero_aterrizaje(3);
aero_aterrizaje(4) = handles.aero_aterrizaje(4);

pesos(1) = handles.pesos(1);
pesos(2) = handles.pesos(2);
pesos(3) = handles.pesos(3);
pesos(4) = handles.pesos(4);

archivo_modelo = strcat(nombre_modelo,'.mat');
fileName = fullfile(handles.dir.data,archivo_modelo);
save(fileName,'propul','aerodinamica','pesos','aero_despegue','aero_aterrizaje');

if modificar_modelo == 0,
load('lista_modelos_avanz.mat');
initial_name=cellstr(lista);
lista = [initial_name;{nombre_modelo}];
fileName2 = fullfile(handles.dir.data,'lista_modelos_avanz.mat');
save(fileName2,'lista');
end

guidata(hObject,handles);

close gcf
seleccion_avanzado;


function edit31_Callback(hObject, eventdata, handles)

function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
handles.propul(5) = get(hObject,'Value');
guidata(hObject,handles);


function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
handles.propul(1) = get(hObject,'Value');
if handles.propul(1) == 1,
set(handles.text3,'String','Empuje a nivel del mar');
set(handles.text33,'String','lbf');
set(handles.text34,'String','lb/(lbf*h)');
set(handles.popupmenu3,'Enable','on');
set(handles.edit10,'Enable','off');
set(handles.text10,'Enable','off');
set(handles.text11,'Enable','on');
else
set(handles.text3,'String','Potencia a nivel del mar');
set(handles.text33,'String','shp');
set(handles.text34,'String','lb/(shp*h)');
set(handles.popupmenu3,'Enable','off');
set(handles.edit10,'Enable','on');
set(handles.text10,'Enable','on');
set(handles.text11,'Enable','off');
end

guidata(hObject,handles);



function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
handles.propul(7) = get(hObject,'Value');
guidata(hObject,handles);

function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dir.main            = pwd;
addpath(handles.dir.main);

%Ruta para la carpeta de las help%
handles.dir.help          = [handles.dir.main,'\help\'];  %new%
addpath(handles.dir.help);  %new%

winopen([handles.dir.help,'help_caracteristicas_avanz.pdf']);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close gcf
seleccion_avanzado;
