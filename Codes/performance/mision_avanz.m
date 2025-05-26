function varargout = mision_avanz(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mision_avanz_OpeningFcn, ...
                   'gui_OutputFcn',  @mision_avanz_OutputFcn, ...
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

function mision_avanz_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
import uiextras.jTree.*
%----------------- CENTRADO EN LA PANTALLA ------------------------------%
scrsz = get(0, 'ScreenSize'); 
pos_act=get(gcf,'Position');   
xr=scrsz(3) - pos_act(3); 
xp=round(xr/2);   
yr=scrsz(4) - pos_act(4); 
yp=round(yr/2);   
set(handles.figure1,'Position',[xp yp pos_act(3) pos_act(4)]);

handles.dir.main            = pwd;
addpath(handles.dir.main);

handles.dir.images          = [handles.dir.main,'\img\'];  
addpath(handles.dir.images);  

handles.dir.gui            = [handles.dir.main,'\gui\'];
addpath(genpath(handles.dir.gui));

handles.dir.data             = [handles.dir.gui,'\data\'];
addpath(handles.dir.data);

set(handles.pushbutton3, 'CData', double(imread('suma.jpg'))/255);
set(handles.pushbutton4, 'CData', double(imread('resta.jpg'))/255);
set(handles.pushbutton5, 'CData', double(imread('llave.jpg'))/255);
lista = get(handles.popupmenu2,'String');
if isempty(lista{1})
else
lista = [{''};lista];
set(handles.popupmenu2,'String',lista);
end

global modificar_mision nombre_mision archivo_mision configurar
configurar.valor = 0;

if modificar_mision == 1
     set(handles.edit1,'String',nombre_mision,'Enable','off');
     set(handles.popupmenu2,'Enable','off');
     set(handles.pushbutton1,'Enable','off');
     set(handles.pushbutton3,'Enable','off');
     set(handles.pushbutton4,'Enable','off');
     set(handles.pushbutton5,'Enable','off');
     set(handles.pushbutton2,'Enable','on');
     modifik = 1;
     archivik = archivo_mision;
 else
     set(handles.edit1,'String','','Enable','on');
     set(handles.popupmenu2,'Enable','off');
     set(handles.pushbutton1,'Enable','off');
     set(handles.pushbutton2,'Enable','off');
     modifik = 0;
end

set(handles.checkbox1,'Visible','off');
set(handles.checkbox2,'Visible','off');

f = gcf;
handles.t = Tree('Parent',f,'FontName','Courier New');
set(handles.t,'Units','pixels','Position',[30 99 330 421]);
handles.patriarca = TreeNode('Name','MISION','Parent',handles.t.Root);
Icon1 = 'puzzle2.png';
setIcon(handles.patriarca,Icon1);
handles.t.Enable = true;
handles.t.DndEnabled = false;
handles.t.RootVisible = false; 

  if modifik == 1, %SI VENIMOS DESDE EL MODIFICADOR DE MISIONES
      load(archivik);
      for i = 1:tramos
      nodo(i) = TreeNode('Name',seg(i).nombre,'Parent',handles.patriarca,'UserData',seg(i).datos);
      setIcon(nodo(i),seg(i).Icon);
        for j = 1:(length(seg(i).nietos)),
            if ~isempty(seg(i).nietos{j}),
            nodos_nietos(j) = TreeNode('Name',seg(i).nietos{j},'Parent',nodo(i));
            setIcon(nodos_nietos(j),'documento.png');
            end
        end
      end
      handles.patriarca.expand();
  end

guidata(hObject,handles);

function varargout = mision_avanz_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global configurar

handles.dir.main            = pwd;
addpath(handles.dir.main);

handles.dir.images          = [handles.dir.main,'\img\'];  
addpath(handles.dir.images);  

handles.dir.gui            = [handles.dir.main,'\gui\'];
addpath(genpath(handles.dir.gui));

handles.dir.data             = [handles.dir.gui,'\data\'];
addpath(handles.dir.data);

set(hObject,'Enable','off');
set(handles.pushbutton3,'Enable','on');
set(handles.pushbutton4,'Enable','on');
set(handles.pushbutton5,'Enable','on');

set(handles.checkbox1,'Visible','off');

valor = get(handles.popupmenu2,'Value');
x = get(handles.popupmenu4,'Value');
y = get(handles.popupmenu4,'Visible');
if x == 1 && strcmp(y,'on')== 1, 
    msgbox('Seleccione el subtipo de segmento','Perfil de la mision','Warn');
    return
end

if valor == 4,
    if x == 9 || x == 10,
        msgbox('Acceso restringido a este subtipo de segmento','Perfil de la mision','Warn');
    return 
    end
end

if valor == 8,
    if x == 9 || x == 10,
        msgbox('Acceso restringido a este subtipo de segmento','Perfil de la mision','Warn');
    return 
    end
end

if configurar.valor == 1,
    delete(configurar.nodos_junior);   
end

import uiextras.jTree.*
valor = get(handles.popupmenu2,'Value');
handles.valor = valor;
lista = get(handles.popupmenu2,'String');
if valor == 4 || valor == 5 || valor == 7 || valor == 8,
    valorz = get(handles.popupmenu4,'Value');
    handles.valorz = valorz;
end

switch valor
    case 2
        handles.taxi(1) = str2double(get(handles.edit10,'String')); % 1: TEMPERATURA LOCAL (K)
        handles.taxi(2) = str2double(get(handles.edit11,'String')); % 2: ALTURA LOCAL (m)
        handles.taxi(3) = str2double(get(handles.edit12,'String')); % 3: PRESION LOCAL (Pa)
        handles.taxi(4) = 0.05;                                     % 4: PALANCA DE RALENTI EN TAXI = 0.05
        handles.taxi(5) = str2double(get(handles.edit13,'String')); % 5: VELOCIDAD A LA QUE HACE EL TAXI (m/s)
        handles.taxi(6) = str2double(get(handles.edit14,'String')); % 6: TIEMPO DE ESPERA EN TAXI (s)
    case 3
        handles.despegue(1) = str2double(get(handles.edit10,'String')); % 1: TEMPERATURA LOCAL (K)
        handles.despegue(2) = str2double(get(handles.edit11,'String')); % 2: ALTURA LOCAL (m)
        handles.despegue(3) = str2double(get(handles.edit12,'String')); % 3: PRESION LOCAL (Pa)
        handles.despegue(4) = str2double(get(handles.edit13,'String')); % 4: COEFICIENTE DE FRICCION CON LA PISTA (MU)
        handles.despegue(5) = 10;                                       % 5: ALTURA DE OBSTACULO (35 FT CIVIL, 50 FT MILITAR)
        handles.despegue(6) = 0.02;                                     % 6: GAMMA DE SUBIDA MINIMO
        handles.despegue(7) = str2double(get(handles.edit14,'String')); % 7: PALANCA DE GASES PARA DESPEGUE
    case 4
        handles.subida.h_inicial = str2double(get(handles.edit10,'String')); % 1: ALTURA INICIAL
        handles.subida.h_final = str2double(get(handles.edit11,'String')); % 1: ALTURA FINAL
        switch valorz
            case 2 % 'Subida dados M y gamma';
                handles.subida.Mach = str2double(get(handles.edit12,'String')); % 3: MACH DE VUELO
                handles.subida.gamma = str2double(get(handles.edit13,'String')); % 2: GAMMA DE SUBIDA
            case 3 % 'Subida dados EAS y gamma';
                handles.subida.EAS = str2double(get(handles.edit12,'String')); % 5: VELOCIDAD EAS
                handles.subida.gamma = str2double(get(handles.edit13,'String')); % 2: GAMMA DE SUBIDA
            case 4 % 'Subida dados TAS y gamma';
                handles.subida.TAS = str2double(get(handles.edit12,'String')); % 4: VELOCIDAD TAS
                handles.subida.gamma = str2double(get(handles.edit13,'String')); % 2: GAMMA DE SUBIDA
            case 5 % 'Subida dados M y palanca';
                handles.subida.Mach = str2double(get(handles.edit12,'String')); % 3: MACH DE VUELO
                handles.subida.palanca = str2double(get(handles.edit13,'String')); % 6: PALANCA DE GASES
            case 6 % 'Subida dados EAS y palanca';
                handles.subida.EAS = str2double(get(handles.edit12,'String')); % 5: VELOCIDAD EAS
                handles.subida.palanca = str2double(get(handles.edit13,'String')); % 6: PALANCA DE GASES
            case 7 % 'Subida dados TAS y palanca';
                handles.subida.TAS = str2double(get(handles.edit12,'String')); % 4: VELOCIDAD TAS
                handles.subida.palanca = str2double(get(handles.edit13,'String')); % 6: PALANCA DE GASES
            case 8 % 'Subida dados V inicial,final y gamma';
                handles.subida.V_ini = str2double(get(handles.edit12,'String')); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA)
                handles.subida.V_fin = str2double(get(handles.edit13,'String')); % 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)
                handles.subida.gamma = str2double(get(handles.edit14,'String')); % 2: GAMMA DE SUBIDA
            case 9 % 'Subida steppest climb';
                handles.subida.palanca = str2double(get(handles.edit12,'String')); % 6: PALANCA DE GASES
            case 10 % 'Subida fastest climb';
                handles.subida.palanca = str2double(get(handles.edit12,'String')); % 6: PALANCA DE GASES
            case 11 % 'Subida dados V inicial,final y palanca'
                handles.subida.V_ini = str2double(get(handles.edit12,'String')); % 7: VELOCIDAD INICIAL (PARA SUBIDA ACELERADA) 
                handles.subida.V_fin = str2double(get(handles.edit13,'String'));% 8: VELOCIDAD FINAL (PARA SUBIDA ACELERADA)
                handles.subida.palanca = str2double(get(handles.edit14,'String')); % 6: PALANCA DE GASES
        end
    case 5
        switch valorz
            case 2
                handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));
                handles.crucero.dist_final = str2double(get(handles.edit11,'String'));
                handles.crucero.Mach = str2double(get(handles.edit12,'String'));
            case 3
                handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));
                handles.crucero.dist_final = str2double(get(handles.edit11,'String'));
                handles.crucero.CL = str2double(get(handles.edit12,'String'));
            case 4
                handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));
                handles.crucero.dist_final = str2double(get(handles.edit11,'String'));
                handles.crucero.V_ini = str2double(get(handles.edit12,'String'));
                handles.crucero.V_fin = str2double(get(handles.edit13,'String'));
                handles.crucero.palanca = str2double(get(handles.edit14,'String'));
            case 5
                handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));
                handles.crucero.dist_final = str2double(get(handles.edit11,'String'));
                handles.crucero.Mach = str2double(get(handles.edit12,'String'));
                handles.crucero.Cd0 = str2double(get(handles.edit13,'String'));
                handles.crucero.k1 = str2double(get(handles.edit14,'String'));
                handles.crucero.k2 = str2double(get(handles.edit15,'String'));
            case 6
                handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));
                handles.crucero.fuel = str2double(get(handles.edit11,'String'));
            case 7
                handles.crucero.h_inicial = str2double(get(handles.edit10,'String'));
                handles.crucero.fuel = str2double(get(handles.edit11,'String'));
        end
    case 6
        handles.soltar_carga.carga = str2double(get(handles.edit10,'String'));
    case 7
        handles.viraje.h_inicial = str2double(get(handles.edit10,'String'));
        handles.viraje.tiempo_final = str2double(get(handles.edit11,'String'));
        switch valorz
            case 2
                handles.viraje.velocidad = str2double(get(handles.edit12,'String'));
                handles.viraje.palanca = str2double(get(handles.edit13,'String'));
            case 3
                handles.viraje.velocidad = str2double(get(handles.edit12,'String'));
                handles.viraje.CL = str2double(get(handles.edit13,'String'));
            case 4
                handles.viraje.velocidad = str2double(get(handles.edit12,'String'));
                handles.viraje.balance = str2double(get(handles.edit13,'String'));
            case 5
                handles.viraje.velocidad = str2double(get(handles.edit12,'String'));
                handles.viraje.n = str2double(get(handles.edit13,'String'));
            case 6
                handles.viraje.velocidad = str2double(get(handles.edit12,'String'));
                handles.viraje.radio = str2double(get(handles.edit13,'String'));
            case 7
                handles.viraje.velocidad = str2double(get(handles.edit12,'String'));
                handles.viraje.vel_guiniada = str2double(get(handles.edit13,'String'));
            case 8
                handles.viraje.palanca = str2double(get(handles.edit12,'String'));
            case 9
                handles.viraje.palanca = str2double(get(handles.edit12,'String'));
            case 10
                handles.viraje.palanca = str2double(get(handles.edit12,'String'));
        end
    case 8
        handles.descenso.h_inicial = str2double(get(handles.edit10,'String'));
        handles.descenso.h_final = str2double(get(handles.edit11,'String'));
        switch valorz
            case 2
                handles.descenso.Mach = str2double(get(handles.edit12,'String')); 
                handles.descenso.gamma = str2double(get(handles.edit13,'String'));
            case 3
                handles.descenso.EAS = str2double(get(handles.edit12,'String')); 
                handles.descenso.gamma = str2double(get(handles.edit13,'String'));
            case 4
                handles.descenso.TAS = str2double(get(handles.edit12,'String')); 
                handles.descenso.gamma = str2double(get(handles.edit13,'String'));
            case 5
                handles.descenso.Mach = str2double(get(handles.edit12,'String')); 
                handles.descenso.palanca = str2double(get(handles.edit13,'String'));
            case 6
                handles.descenso.EAS = str2double(get(handles.edit12,'String')); 
                handles.descenso.palanca = str2double(get(handles.edit13,'String'));
            case 7
                handles.descenso.TAS = str2double(get(handles.edit12,'String')); 
                handles.descenso.palanca = str2double(get(handles.edit13,'String'));
            case 8
                handles.descenso.V_ini = str2double(get(handles.edit12,'String')); 
                handles.descenso.V_fin = str2double(get(handles.edit13,'String'));
                handles.descenso.gamma = str2double(get(handles.edit14,'String'));
            case 9
                handles.descenso.palanca = str2double(get(handles.edit12,'String'));
            case 10
                handles.descenso.palanca = str2double(get(handles.edit12,'String'));
            case 11
                handles.descenso.V_ini = str2double(get(handles.edit12,'String')); 
                handles.descenso.V_fin = str2double(get(handles.edit13,'String'));
                handles.descenso.palanca = str2double(get(handles.edit14,'String'));
        end
    case 9
        handles.aterrizaje(1) = str2double(get(handles.edit10,'String'));
        handles.aterrizaje(2) = str2double(get(handles.edit11,'String'));
        handles.aterrizaje(3) = str2double(get(handles.edit12,'String'));
        handles.aterrizaje(4) = str2double(get(handles.edit13,'String'));
        handles.aterrizaje(5) = str2double(get(handles.edit14,'String'));
        handles.aterrizaje(6) = str2double(get(handles.edit15,'String'));
end

visiblez = get(handles.popupmenu4,'Visible');
visible0 = get(handles.edit10,'Visible');
visible1 = get(handles.edit11,'Visible');
visible2 = get(handles.edit12,'Visible');
visible3 = get(handles.edit13,'Visible');
visible4 = get(handles.edit14,'Visible');
visible5 = get(handles.edit15,'Visible');

if strcmp(visible0,'on') == 1,
    valor0 = get(handles.edit10,'String');
    handles.valor0 = valor0;
else
    handles.valor0 = -100;
end
if strcmp(visible1,'on') == 1,
    valor1 = get(handles.edit11,'String');
    handles.valor1 = valor1;
else
    handles.valor1 = -100;
end
if strcmp(visible2,'on') == 1,
    valor2 = get(handles.edit12,'String');
    handles.valor2 = valor2;
else
    handles.valor2 = -100;
end
if strcmp(visible3,'on') == 1,
    valor3 = get(handles.edit13,'String');
    handles.valor3 = valor3;
else
    handles.valor3 = -100;
end
if strcmp(visible4,'on') == 1,
    valor4 = get(handles.edit14,'String');
    handles.valor4 = valor4;
else
    handles.valor4 = -100;
end
if strcmp(visible5,'on') == 1,
    valor5 = get(handles.edit15,'String');
    handles.valor5 = valor5;
else
    handles.valor5 = -100;
end

guidata(hObject,handles);

Icon = 'documento.png';
nombre = lista{valor};
if configurar.valor == 0
nodo = TreeNode('Name',nombre,'Parent',handles.patriarca,'UserData',handles);
else
nodo = configurar.nodo;
set(nodo,'UserData',handles);
end

if strcmp(visiblez,'on') == 1,
    listaz = get(handles.popupmenu4,'String');
    nodoz = TreeNode('Name',listaz{valorz},'Parent',nodo,'Value',valorz);
    setIcon(nodoz,Icon);
end
if strcmp(visible0,'on') == 1,
    nombre0 = get(handles.text10,'String');
    valor0 = get(handles.edit10,'String');
    handles.valor0 = valor0;
    unidad0 = get(handles.text20,'String');
    nombre0 = [nombre0,':',blanks(20 - length(nombre0)),valor0,blanks(1),unidad0];
    nodo0 = TreeNode('Name',nombre0,'Parent',nodo,'Value',valor0);
    setIcon(nodo0,Icon);
end
if strcmp(visible1,'on') == 1,
    nombre1 = get(handles.text11,'String');
    valor1 = get(handles.edit11,'String');
    handles.valor1 = valor1;
    unidad1 = get(handles.text21,'String');
    nombre1 = [nombre1,':',blanks(20 - length(nombre1)),valor1,blanks(1),unidad1];
    nodo1 = TreeNode('Name',nombre1,'Parent',nodo,'Value',valor1);
    setIcon(nodo1,Icon);
end
if strcmp(visible2,'on') == 1,
    nombre2 = get(handles.text12,'String');
    valor2 = get(handles.edit12,'String');
    handles.valor2 = valor2;
    unidad2 = get(handles.text22,'String');
    nombre2 = [nombre2,':',blanks(20 - length(nombre2)),valor2,blanks(1),unidad2];
    nodo2 = TreeNode('Name',nombre2,'Parent',nodo,'Value',valor2);
    setIcon(nodo2,Icon);
end
if strcmp(visible3,'on') == 1,
    nombre3 = get(handles.text13,'String');
    valor3 = get(handles.edit13,'String');
    handles.valor3 = valor3;
    unidad3 = get(handles.text23,'String');
    nombre3 = [nombre3,':',blanks(20 - length(nombre3)),valor3,blanks(1),unidad3];
    nodo3 = TreeNode('Name',nombre3,'Parent',nodo,'Value',valor3);
    setIcon(nodo3,Icon);
end
if strcmp(visible4,'on') == 1,
    nombre4 = get(handles.text14,'String');
    valor4 = get(handles.edit14,'String');
    handles.valor4 = valor4;
    unidad4 = get(handles.text24,'String');
    nombre4 = [nombre4,':',blanks(20 - length(nombre4)),valor4,blanks(1),unidad4];
    nodo4 = TreeNode('Name',nombre4,'Parent',nodo,'Value',valor4);
    setIcon(nodo4,Icon);
end
if strcmp(visible5,'on') == 1,
    nombre5 = get(handles.text15,'String');
    valor5 = get(handles.edit15,'String');
    handles.valor5 = valor5;
    unidad5 = get(handles.text25,'String');
    nombre5 = [nombre5,':',blanks(20 - length(nombre5)),valor5,blanks(1),unidad5];
    nodo5 = TreeNode('Name',nombre5,'Parent',nodo,'Value',valor5);
    setIcon(nodo5,Icon);
end
    

switch valor
    case 2
        Icon = 'taxi.png';
    case 3
        Icon = 'takeoff.png';
    case 4
        Icon = 'climb.png';
    case 5
        Icon = 'crucero.png';
    case 6
        Icon = 'paquete.png';
    case 7
        Icon = 'viraje.png';
    case 8
        Icon = 'descend.png';
    case 9
        Icon = 'landing.png';
end
setIcon(nodo,Icon);
handles.patriarca.expand();
handles.t.Enable = true;
handles.t.DndEnabled = true;
guidata(hObject,handles);


set(handles.text10,'Visible','off'); set(handles.edit10,'Visible','off');  set(handles.text20,'Visible','off');
set(handles.text11,'Visible','off'); set(handles.edit11,'Visible','off');  set(handles.text21,'Visible','off');
set(handles.text12,'Visible','off'); set(handles.edit12,'Visible','off');  set(handles.text22,'Visible','off');
set(handles.text13,'Visible','off'); set(handles.edit13,'Visible','off');  set(handles.text23,'Visible','off');
set(handles.text14,'Visible','off'); set(handles.edit14,'Visible','off');  set(handles.text24,'Visible','off');
set(handles.text15,'Visible','off'); set(handles.edit15,'Visible','off');  set(handles.text25,'Visible','off');
set(handles.text98,'Visible','off');
set(handles.popupmenu4,'Visible','off');
set(handles.popupmenu4,'Value',1);
set(handles.popupmenu2,'Value',1);
set(handles.popupmenu2,'Enable','off');
set(handles.pushbutton2,'Enable','on');
set(handles.checkbox1,'Visible','off');



set(handles.pushbutton5,'Enable','on');
 configurar.valor = 0;
 
%------------------------------------------------------------------------% 

function pushbutton2_Callback(hObject, eventdata, handles)
global modificar_mision nombre_mision archivo_mision
if modificar_mision == 0,    
nodos_segmentos = get(handles.patriarca,'Children');
if numel(nodos_segmentos) <= 1,
     msgbox('El numero de segmentos introducido es insuficiente','Perfil de mision','Warn');
     return;
end
tramos = length(nodos_segmentos);
nodos_nietos = get(nodos_segmentos,'Children');
for kk = 1:length(nodos_nietos),
    n(kk) = length(nodos_nietos{kk,1});
end
nombre_nietos = cell(tramos,max(n));
for i = 1:tramos,
    for j = 1:n(i)
        nombre_nietos{i,j} = nodos_nietos{i,1}(1,j).Name;
    end
end
datos_segmentos = get(nodos_segmentos,'UserData');
nombre_segmentos = get(nodos_segmentos,'Name');

for i = 1:tramos,
    if strcmp(nombre_segmentos(i),'Taxi'),
        seg(i).nombre = 'Taxi';
        seg(i).datos = datos_segmentos{i}.taxi;
        seg(i).Icon = 'taxi.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strcmp(nombre_segmentos(i),'Despegue'),
        seg(i).nombre = 'Despegue';
        seg(i).datos = datos_segmentos{i}.despegue;
        seg(i).Icon = 'takeoff.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strcmp(nombre_segmentos(i),'Subida'),
        seg(i).nombre = 'Subida';
        seg(i).opcion = (datos_segmentos{i}.valorz)-1;
        seg(i).datos = datos_segmentos{i}.subida;
        seg(i).Icon = 'climb.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strcmp(nombre_segmentos(i),'Crucero'),
        seg(i).nombre = 'Crucero';
        seg(i).opcion = (datos_segmentos{i}.valorz)-1;
        seg(i).datos = datos_segmentos{i}.crucero;
        seg(i).Icon = 'crucero.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strcmp(nombre_segmentos(i),'Soltar carga'),
        seg(i).nombre = 'Soltar carga';
        seg(i).datos = datos_segmentos{i}.soltar_carga;
        seg(i).Icon = 'paquete.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strcmp(nombre_segmentos(i),'Viraje'),
        seg(i).nombre = 'Viraje';
        seg(i).opcion = (datos_segmentos{i}.valorz)-1;
        seg(i).datos = datos_segmentos{i}.viraje;
        seg(i).Icon = 'viraje.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strcmp(nombre_segmentos(i),'Descenso'),
        seg(i).nombre = 'Descenso';
        seg(i).opcion = (datos_segmentos{i}.valorz)-1;
        seg(i).datos = datos_segmentos{i}.descenso;
        seg(i).Icon = 'descend.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
    if strncmp(nombre_segmentos{i},'Aterrizaje',3),
        seg(i).nombre = 'Aterrizaje';
        seg(i).datos = datos_segmentos{i}.aterrizaje;
        seg(i).Icon = 'landing.png';
        seg(i).nietos = nombre_nietos(i,:);
    end
end

nombre_mision = get(handles.edit1,'String');
if isempty(nombre_mision),
    msgbox('El nombre introducido no es válido','Perfil de la mision','Warn'); 
    return
end

handles.dir.main            = pwd;
addpath(handles.dir.main);

handles.dir.gui             = [handles.dir.main,'\gui\'];
addpath(genpath(handles.dir.gui));

handles.dir.data            = [handles.dir.gui,'\data\'];

modif = get(handles.edit1,'Enable');
if strcmp(modif,'on'),
mision = get(handles.edit1,'String');
load('lista_misiones_avanz.mat');
post = strcmp(mision,lista);
if sum(post) > 0
    msgbox('El nombre introducido ya está creado','Perfil de mision','Warn'); 
    return
end
end

archivo_mision = strcat(nombre_mision,'.mat');
arbol = handles.t;
fileName = fullfile(handles.dir.data,archivo_mision);
save(fileName,'seg','tramos');


load('lista_misiones_avanz.mat');
initial_name=cellstr(lista);
lista = [initial_name;{nombre_mision}];
fileName2 = fullfile(handles.dir.data,'lista_misiones_avanz.mat');
save(fileName2,'lista');


guidata(hObject,handles);

end

close gcf
seleccion_avanzado;


function popupmenu2_Callback(hObject, eventdata, handles)
set(handles.text10,'Visible','off'); set(handles.edit10,'Visible','off');  set(handles.text20,'Visible','off');
set(handles.text11,'Visible','off'); set(handles.edit11,'Visible','off');  set(handles.text21,'Visible','off');
set(handles.text12,'Visible','off'); set(handles.edit12,'Visible','off');  set(handles.text22,'Visible','off');
set(handles.text13,'Visible','off'); set(handles.edit13,'Visible','off');  set(handles.text23,'Visible','off');
set(handles.text14,'Visible','off'); set(handles.edit14,'Visible','off');  set(handles.text24,'Visible','off');
set(handles.text15,'Visible','off'); set(handles.edit15,'Visible','off');  set(handles.text25,'Visible','off');
set(handles.text98,'Visible','off');
set(handles.popupmenu4,'Visible','off');
set(handles.popupmenu4,'Value',1);
set(handles.pushbutton1,'Enable','on');
set(handles.checkbox1,'Visible','off');
set(handles.checkbox1,'Value',0);
set(handles.checkbox2,'Visible','off');
set(handles.checkbox2,'Value',0);

        set(handles.uipanel2,'Visible','off');
        set(handles.text40,'Visible','off');
        set(handles.text41,'Visible','off');
        set(handles.text42,'Visible','off');
        set(handles.edit8,'Visible','off');
        set(handles.edit9,'Visible','off');
        set(handles.edit40,'Visible','off');
        

valor = get(hObject,'Value');
switch valor
    case 2
        set(handles.text10,'String','Temperatura local'); set(handles.edit10,'String','0');  set(handles.text20,'String','[K]');
        set(handles.text11,'String','Altura local');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Presion local');     set(handles.edit12,'String','0');  set(handles.text22,'String','[Pa]');
        set(handles.text13,'String','Velocidad de taxi'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Tiempo de espera');  set(handles.edit14,'String','0');  set(handles.text24,'String','[s]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
        
        set(handles.checkbox1,'Visible','off');
        set(handles.checkbox2,'Visible','off');
        set(handles.text40,'Visible','off');
        set(handles.text41,'Visible','off');
        set(handles.text42,'Visible','off');
        set(handles.edit8,'Visible','off');
        set(handles.edit9,'Visible','off');
        set(handles.edit40,'Visible','off');
        set(handles.uipanel2,'Visible','off');
       
    case 3   
        set(handles.text10,'String','Temperatura local'); set(handles.edit10,'String','0');  set(handles.text20,'String','[K]');
        set(handles.text11,'String','Altura local');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Presion local');     set(handles.edit12,'String','0');  set(handles.text22,'String','[Pa]');
        set(handles.text13,'String','Coef de friccion'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        set(handles.text14,'String','Palanca de gases');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
        
        set(handles.checkbox1,'Visible','on');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Visible','off');
        set(handles.checkbox2,'Value',0);
        set(handles.text40,'Visible','off');
        set(handles.text41,'Visible','off');
        set(handles.text42,'Visible','off');
        set(handles.edit8,'Visible','off');
        set(handles.edit9,'Visible','off');
        set(handles.edit40,'Visible','off');
        set(handles.uipanel2,'Visible','off');
        
    case 4
               
        set(handles.text98,'Visible','on');
        set(handles.popupmenu4,'Visible','on');
        lista = {'';'Subida dados M y gamma';'Subida dados EAS y gamma';'Subida dados TAS y gamma';...
                 'Subida dados M y palanca';'Subida dados EAS y palanca';'Subida dados TAS y palanca';...
                 'Subida dados V inicial,final y gamma';'Subida steppest climb';'Subida fastest climb';...
                 'Subida dados V inicial,final y palanca'};
        set(handles.popupmenu4,'String',lista);
        set(handles.checkbox1,'Visible','on');        
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Visible','on');
        set(handles.checkbox2,'Value',0);

        
        
    case 5        
        set(handles.text98,'Visible','on');
        set(handles.popupmenu4,'Visible','on');
        lista = {'';'Crucero dado M y distancia';'Crucero dado CL y distancia';'Crucero dados V inicial,final y palanca';...
                 'Crucero con polar en funcion de M';'Crucero de max alcance dado peso final';'Crucero de max autonomia dado peso final'};       
        set(handles.popupmenu4,'String',lista);
        
        set(handles.checkbox1,'Visible','on');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Visible','off');
        set(handles.checkbox2,'Value',0);
        set(handles.uipanel2,'Visible','off');
        
    case 6
        
        set(handles.text10,'String','Carga lanzada'); set(handles.edit10,'String','0');  set(handles.text20,'String','kg');
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
    
        set(handles.checkbox1,'Visible','off');
        set(handles.checkbox2,'Visible','off');
        set(handles.checkbox2,'Value',0);
        set(handles.uipanel2,'Visible','off');
    case 7
        set(handles.text98,'Visible','on');
        set(handles.popupmenu4,'Visible','on');
        lista = {'';'Viraje horizontal dado V y palanca';'Viraje horizontal dado V y CL';'Viraje horizontal dado V y balance';...
                 'Viraje horizontal dado V y n';'V.H dado V y radio de giro ';'V.H dado V y velocidad de guiñada';...
                 'V.H dado palanca y a factor de carga max';'V.H dado palanca y a v de guiñada max';'V.H dado palanca y a radio min'};
        set(handles.popupmenu4,'String',lista);
        
        set(handles.checkbox1,'Visible','on');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Visible','on');
        set(handles.checkbox2,'Value',0);
       
    case 8
        set(handles.text98,'Visible','on');
        set(handles.popupmenu4,'Visible','on');
        lista = {'';'Descenso dados M y gamma';'Descenso dados EAS y gamma';'Descenso dados TAS y gamma';...
                 'Descenso dados M y palanca';'Descenso dados EAS y palanca';'Descenso dados TAS y palanca';...
                 'Descenso dados V inicial,final y gamma';'Descenso a minimo gamma';'Descenso "slowest sink"';...
                 'Descenso dados V inicial,final y palanca'};
        set(handles.popupmenu4,'String',lista);
        
        set(handles.checkbox1,'Visible','on');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Visible','on');
        set(handles.checkbox2,'Value',0);
        
    case 9
        set(handles.text10,'String','Temperatura local'); set(handles.edit10,'String','0');  set(handles.text20,'String','[K]');
        set(handles.text11,'String','Altura local');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Presion local');     set(handles.edit12,'String','0');  set(handles.text22,'String','[Pa]');
        set(handles.text13,'String','Coef de friccion'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        set(handles.text14,'String','Palanca de reversa');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');
        set(handles.text15,'String','Tiempo en activar frenos');  set(handles.edit15,'String','0');  set(handles.text25,'String','[s]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
        set(handles.text15,'Visible','on'); set(handles.edit15,'Visible','on');  set(handles.text25,'Visible','on');
        
        set(handles.checkbox1,'Visible','on');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Visible','off');
        set(handles.checkbox2,'Value',0);
        set(handles.uipanel2,'Visible','off');
end



function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.

function pushbutton3_Callback(hObject, eventdata, handles)
set(handles.popupmenu2,'Enable','on');

set(handles.pushbutton4,'Enable','off');
set(handles.pushbutton5,'Enable','off');
global configurar
configurar.valor = 0;




function pushbutton4_Callback(hObject, eventdata, handles)
nodo = get(handles.t,'SelectedNodes');
try
logic = nodo.Parent;
if strcmp(logic.Name,'Root') == 1,
    return;
else
delete(nodo);
end
catch
    msgbox('No se ha seleccionado ningun segmento valido','Seleccion de mision','Warn');
    return
end



function popupmenu4_Callback(hObject, eventdata, handles)

set(handles.pushbutton1,'Enable','on');
set(handles.text10,'Visible','off'); set(handles.edit10,'Visible','off');  set(handles.text20,'Visible','off');
set(handles.text11,'Visible','off'); set(handles.edit11,'Visible','off');  set(handles.text21,'Visible','off');
set(handles.text12,'Visible','off'); set(handles.edit12,'Visible','off');  set(handles.text22,'Visible','off');
set(handles.text13,'Visible','off'); set(handles.edit13,'Visible','off');  set(handles.text23,'Visible','off');
set(handles.text14,'Visible','off'); set(handles.edit14,'Visible','off');  set(handles.text24,'Visible','off');
set(handles.text15,'Visible','off'); set(handles.edit15,'Visible','off');  set(handles.text25,'Visible','off');

valor = get(handles.popupmenu2,'Value');
valor2 = get(hObject,'Value');
switch valor
    case 4 %Subida
        switch valor2
            case 2 %M & gamma
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Mach de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        set(handles.text13,'String','Gamma de subida'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 3 %EAS y gamma
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','EAS de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Gamma de subida'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');       
            case 4 %TAS Y GAMMA
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Gamma de subida'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
            case 5 %M y palanca
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Mach de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
            case 6 %EAS y palanca
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','EAS de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 7 %V Y PALANCA
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 8 %V0 vf Y gamma
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad inicial');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Velocidad final'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Gamma de subida');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');

        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
            case 9 %Steppest climb
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Palanca de gases'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 10 %Fastest climb
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Palanca de gases'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 11 %VO VF PALANCA
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad inicial');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Velocidad final'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Palanca de gases');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');

        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
        end
        
    case 5 %Crucero
        switch valor2
            case 2 %M & distancia
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Distancia final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Mach de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 3 %Cl y distancia
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Distancia final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','CL de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 4 %V0 vf Y palanca
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Distancia final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad inicial');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Velocidad final'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Palanca de gases');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');

        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
            case 5 %Polar = f(M)
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Distancia final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Mach de vuelo'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        set(handles.text13,'String','Cd0 a Mach de vuelo'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        set(handles.text14,'String','k1 a Mach de vuelo');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');
        set(handles.text15,'String','k2 a Mach de vuelo');  set(handles.edit15,'String','0');  set(handles.text25,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
        set(handles.text15,'Visible','on'); set(handles.edit15,'Visible','on');  set(handles.text25,'Visible','on');
        
            case 6 %Max alcance dado fuel
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Fuel a consumir');      set(handles.edit11,'String','0');  set(handles.text21,'String','[kg]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
            case 7 %Max autonomia dado fuel
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Fuel a consumir');      set(handles.edit11,'String','0');  set(handles.text21,'String','[kg]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        end
    case 7 %Viraje
        switch valor2
            case 2 %V y palanca
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Velocidad de viraje');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 3 %V y CL
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Velocidad de viraje');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','CL en el viraje'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');       
            case 4 %V y mu
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Velocidad de viraje');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Angulo de balance'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
            case 5 %V y n
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Velocidad de viraje');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Factor de carga'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
            case 6 %V y radio
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Velocidad de viraje');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Radio de giro'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 7 %V Y velocidad de guiniada
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Velocidad de viraje');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Velocidad de guiñada'); set(handles.edit13,'String','0');  set(handles.text23,'String','[rad/s]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 8 %Palanca y nmax
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Palanca de gases');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 9 %Palanca y velocidad de guiñada max
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Palanca de gases'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 10 %Palanca y radio minimo
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Tiempo final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[s]');
        set(handles.text12,'String','Palanca de gases'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        end
    case 8  %Descenso
        switch valor2
            case 2 %M & gamma
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Mach de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        set(handles.text13,'String','Gamma de descenso'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 3 %EAS y gamma
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','EAS de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Gamma de descenso'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');       
            case 4 %TAS Y GAMMA
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Gamma de descenso'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
            case 5 %M y palanca
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Mach de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
            case 6 %EAS y palanca
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','EAS de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 7 %V Y PALANCA
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad de vuelo');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Palanca de gases'); set(handles.edit13,'String','0');  set(handles.text23,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on');
            case 8 %V0 vf Y gamma
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad inicial');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Velocidad final'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Gamma de descenso');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');

        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
            case 9 %gamma min
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Palanca de gases'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 10 %slowest sink
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Palanca de gases'); set(handles.edit12,'String','0');  set(handles.text22,'String','[-]');
        
        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
            case 11 %VO VF PALANCA
        set(handles.text10,'String','Altura inicial'); set(handles.edit10,'String','0');  set(handles.text20,'String','[m]');
        set(handles.text11,'String','Altura final');      set(handles.edit11,'String','0');  set(handles.text21,'String','[m]');
        set(handles.text12,'String','Velocidad inicial');     set(handles.edit12,'String','0');  set(handles.text22,'String','[m/s]');
        set(handles.text13,'String','Velocidad final'); set(handles.edit13,'String','0');  set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Palanca de gases');  set(handles.edit14,'String','0');  set(handles.text24,'String','[-]');

        set(handles.text10,'Visible','on'); set(handles.edit10,'Visible','on');  set(handles.text20,'Visible','on');
        set(handles.text11,'Visible','on'); set(handles.edit11,'Visible','on');  set(handles.text21,'Visible','on');
        set(handles.text12,'Visible','on'); set(handles.edit12,'Visible','on');  set(handles.text22,'Visible','on');
        set(handles.text13,'Visible','on'); set(handles.edit13,'Visible','on');  set(handles.text23,'Visible','on'); 
        set(handles.text14,'Visible','on'); set(handles.edit14,'Visible','on');  set(handles.text24,'Visible','on');
        end
end
        
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)

try    
nodo = get(handles.t,'SelectedNodes');
nodos_junior = get(nodo,'Children');
set(handles.popupmenu2,'Value',nodo.UserData.valor,'Visible','on');
catch
    msgbox('No se ha seleccionado ningun segmento valido','Seleccion de mision','Warn');
    return
end

set(handles.pushbutton3,'Enable','off');
set(handles.pushbutton4,'Enable','off');
set(handles.pushbutton1,'Enable','on');

set(handles.text10,'Visible','off'); set(handles.edit10,'Visible','off');  set(handles.text20,'Visible','off');
set(handles.text11,'Visible','off'); set(handles.edit11,'Visible','off');  set(handles.text21,'Visible','off');
set(handles.text12,'Visible','off'); set(handles.edit12,'Visible','off');  set(handles.text22,'Visible','off');
set(handles.text13,'Visible','off'); set(handles.edit13,'Visible','off');  set(handles.text23,'Visible','off');
set(handles.text14,'Visible','off'); set(handles.edit14,'Visible','off');  set(handles.text24,'Visible','off');
set(handles.text15,'Visible','off'); set(handles.edit15,'Visible','off');  set(handles.text25,'Visible','off');
set(handles.text98,'Visible','off');
set(handles.popupmenu4,'Visible','off');

if nodo.UserData.valor == 2,
        set(handles.text10,'String','Temperatura local');   set(handles.text20,'String','[K]');
        set(handles.text11,'String','Altura local');        set(handles.text21,'String','[m]');
        set(handles.text12,'String','Presion local');       set(handles.text22,'String','[Pa]');
        set(handles.text13,'String','Velocidad de taxi');   set(handles.text23,'String','[m/s]');
        set(handles.text14,'String','Tiempo de espera');    set(handles.text24,'String','[s]');
end
if nodo.UserData.valor == 3,
        set(handles.text10,'String','Temperatura local');   set(handles.text20,'String','[K]');
        set(handles.text11,'String','Altura local');        set(handles.text21,'String','[m]');
        set(handles.text12,'String','Presion local');      set(handles.text22,'String','[Pa]');
        set(handles.text13,'String','Coef de friccion');   set(handles.text23,'String','[-]');
        set(handles.text14,'String','Palanca de gases');    set(handles.text24,'String','[-]');
end        
if nodo.UserData.valor == 4,
    lista = {'';'Subida dados M y gamma';'Subida dados EAS y gamma';'Subida dados TAS y gamma';...
                 'Subida dados M y palanca';'Subida dados EAS y palanca';'Subida dados TAS y palanca';...
                 'Subida dados V inicial,final y gamma';'Subida steppest climb';'Subida fastest climb';...
                 'Subida dados V inicial,final y palanca'};
set(handles.popupmenu4,'String',lista);
set(handles.popupmenu4,'Value',nodo.UserData.valorz,'Visible','on');
mision_avanz('popupmenu4_Callback',handles.popupmenu4,'',handles)
set(handles.text98,'Visible','on');             
end
if nodo.UserData.valor == 5,
    lista = {'';'Crucero dado M y distancia';'Crucero dado CL y distancia';'Crucero dados V inicial,final y palanca';...
                 'Crucero con polar en funcion de M';'Crucero de max alcance dado peso final';'Crucero de max autonomia dado peso final'};
set(handles.popupmenu4,'String',lista);
set(handles.popupmenu4,'Value',nodo.UserData.valorz,'Visible','on');
mision_avanz('popupmenu4_Callback',handles.popupmenu4,'',handles)
set(handles.text98,'Visible','on');             
end
if nodo.UserData.valor == 6,
set(handles.text10,'String','Carga lanzada');   set(handles.text20,'String','kg');
end
if nodo.UserData.valor == 7,
lista = {'';'Viraje horizontal dado V y palanca';'Viraje horizontal dado V y CL';'Viraje horizontal dado V y balance';...
                 'Viraje horizontal dado V y n';'V.H dado V y radio de giro ';'V.H dado V y velocidad de guiñada';...
                 'V.H dado palanca y a factor de carga max';'V.H dado palanca y a v de guiñada max';'V.H dado palanca y a radio min'};
set(handles.popupmenu4,'String',lista);
set(handles.popupmenu4,'Value',nodo.UserData.valorz,'Visible','on');
mision_avanz('popupmenu4_Callback',handles.popupmenu4,'',handles)
set(handles.text98,'Visible','on');             
end
if nodo.UserData.valor == 8,
lista = {'';'Descenso dados M y gamma';'Descenso dados EAS y gamma';'Descenso dados TAS y gamma';...
                 'Descenso dados M y palanca';'Descenso dados EAS y palanca';'Descenso dados TAS y palanca';...
                 'Descenso dados V inicial,final y gamma';'Descenso a minimo gamma';'Descenso "slowest sink"';...
                 'Descenso dados V inicial,final y palanca'};
set(handles.popupmenu4,'String',lista);
set(handles.popupmenu4,'Value',nodo.UserData.valorz,'Visible','on');
mision_avanz('popupmenu4_Callback',handles.popupmenu4,'',handles)
set(handles.text98,'Visible','on');
end
if nodo.UserData.valor == 9,
        set(handles.text10,'String','Temperatura local');   set(handles.text20,'String','[K]');
        set(handles.text11,'String','Altura local');      set(handles.text21,'String','[m]');
        set(handles.text12,'String','Presion local');       set(handles.text22,'String','[Pa]');
        set(handles.text13,'String','Coef de friccion');   set(handles.text23,'String','[-]');
        set(handles.text14,'String','Palanca de reversa');    set(handles.text24,'String','[-]');
        set(handles.text15,'String','Tiempo en activar frenos');   set(handles.text25,'String','[s]');        
end


try
if nodo.UserData.valor0 > -10,
set(handles.edit10,'String',nodo.UserData.valor0,'Visible','on');
set(handles.text10,'Visible','on'); set(handles.text20,'Visible','on');
end
if nodo.UserData.valor1 > -10,
set(handles.edit11,'String',nodo.UserData.valor1,'Visible','on');
set(handles.text11,'Visible','on'); set(handles.text21,'Visible','on');
end
if nodo.UserData.valor2 > -10,
set(handles.edit12,'String',nodo.UserData.valor2,'Visible','on');
set(handles.text12,'Visible','on'); set(handles.text22,'Visible','on');
end
if nodo.UserData.valor3 > -10,
set(handles.edit13,'String',nodo.UserData.valor3,'Visible','on');
set(handles.text13,'Visible','on'); set(handles.text23,'Visible','on');
end
if nodo.UserData.valor4 > -10,
set(handles.edit14,'String',nodo.UserData.valor4,'Visible','on');
set(handles.text14,'Visible','on'); set(handles.text24,'Visible','on');
end
if nodo.UserData.valor5 > -10,
set(handles.edit15,'String',nodo.UserData.valor5,'Visible','on');
set(handles.text15,'Visible','on'); set(handles.text25,'Visible','on');
end
catch 
end

global configurar
configurar.valor = 1;
configurar.nodo = nodo;
configurar.nodos_junior = nodos_junior;
set(hObject,'Enable','inactive');


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

winopen([handles.dir.help,'help_mision_avanz.pdf']);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close gcf
seleccion_avanzado;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


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



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value=get(handles.checkbox2,'Value');
if value==1,
        set(handles.uipanel2,'Visible','on');
        set(handles.text40,'Visible','on');
        set(handles.text41,'Visible','on');
        set(handles.text42,'Visible','on');
        set(handles.edit8,'Visible','on');
        set(handles.edit9,'Visible','on');
        set(handles.edit40,'Visible','on');
else
        set(handles.uipanel2,'Visible','off');
        set(handles.text40,'Visible','off');
        set(handles.text41,'Visible','off');
        set(handles.text42,'Visible','off');
        set(handles.edit8,'Visible','off');
        set(handles.edit9,'Visible','off');
        set(handles.edit40,'Visible','off');
end
    
% Hint: get(hObject,'Value') returns toggle state of checkbox2