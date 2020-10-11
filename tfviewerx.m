function varargout = tfviewerx(varargin)
% TFVIEWERX MATLAB code for tfviewerx.fig
%    USAGE:
% 	tfviewerx(time,freq,DATA,chanlocs[,title]);
%
%	   time     : vector of time points
%	   freq     : vector of frequencies (e.g., of wavelet convolution)
%	   DATA     : channels X frequency X time matrix
%	   chanlocs : eeglab chanlocs structure
%	   title    : (optional) string title for figure
%
% written by mikexcohen@gmail.com

%      TFVIEWERX, by itself, creates a new TFVIEWERX or raises the existing
%      singleton*.
%
%      H = TFVIEWERX returns the handle to a new TFVIEWERX or the handle to
%      the existing singleton*.
%
%      TFVIEWERX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TFVIEWERX.M with the given input arguments.
%
%      TFVIEWERX('Property','Value',...) creates a new TFVIEWERX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tfviewerx_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tfviewerx_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 05-Jun-2012 15:37:29

%% 
if strcmp(get(0,'DefaultFigureWindowStyle'),'docked'), set(0,'DefaultFigureWindowStyle','normal'); end


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tfviewerx_OpeningFcn, ...
                   'gui_OutputFcn',  @tfviewerx_OutputFcn, ...
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


% --- Executes just before tfviewerx is made visible.
function tfviewerx_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tfviewerx (see VARARGIN)

% Choose default command line output for tfviewerx
handles.output = hObject;

% check inputs
if nargin<4
    help tfviewerx
    warning('No data provided; using pre-loaded data')
    
    tempdata      = get(handles.figure1,'UserData');
    data.t        = tempdata.times;
    data.f        = tempdata.frex;
    data.ctf      = tempdata.tf;
    data.chanlocs = tempdata.chanlocs;
else
    data.t        = varargin{1};
    data.f        = varargin{2};
    data.ctf      = varargin{3};
    data.chanlocs = varargin{4};
end;


if length(data.chanlocs)>100
    set(handles.topo_plotrad,'string','.75');
else
    set(handles.topo_plotrad,'string','.55');
end

set(handles.save_filepath,'string',[ pwd '/' ]);

if length(data.t)~=size(data.ctf,3)
    error('Something wrong with time...')
end
if length(data.f)~=size(data.ctf,2)                  
    error('Something wrong with frequency...')
end
if ~isfield(data.chanlocs,'X')
    error('Something wrong with chanlocs...')
end



set(handles.figure1,'UserData',data);
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;

%% plot data

set(handles.tf_pickelectrode,'String',{data.chanlocs.labels});

% plot TF
yskip=5;
[junk,handles.tfobject]=contourf(handles.tf_axes,t,f,squeeze(ctf(get(handles.tf_pickelectrode,'Value'),:,:)),40,'linecolor','none');
set(handles.tf_axes,'yscale','log','ytick',round(logspace(log10(f(1)),log10(f(end)),length(f(1:yskip:end)))*10)/10);
hold(handles.tf_axes,'on')

handles.tflineV=plot(handles.tf_axes,[str2double(get(handles.tf_time2plot,'String')) str2double(get(handles.tf_time2plot,'String'))],get(handles.tf_axes,'ylim'),'k');
handles.tflineH=plot(handles.tf_axes,get(handles.tf_axes,'xlim'),[str2double(get(handles.tf_freq2plot,'String')) str2double(get(handles.tf_freq2plot,'String'))],'k');

% update xlim
set(handles.tf_xlim,'string',[ num2str(round(t(1))) ' ' num2str(round(t(end))) ]);

% plot topo
handles=update_topomap(handles); guidata(handles.figure1,handles); 
update_tfmap(handles)

% Update handles structure
guidata(handles.figure1,handles);
set(handles.tfobject,'buttondownfcn',{@tf_axes_ButtonDownFcn,handles});

% title, maybe
if length(varargin)==5
    set(handles.figure1,'name',[ 'tfviewerx:   ' varargin{5} ]);
end

%%

% UIWAIT makes tfviewerx wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tfviewerx_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function tf_clim_Callback(hObject, eventdata, handles)
update_tfmap(handles);


% --- Executes during object creation, after setting all properties.
function tf_clim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_clim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savefigure.
function savefigure_Callback(hObject, eventdata, handles)

fpath = get(handles.save_filepath,'String');
if ~any(fpath(end)=='\/'), fpath=[fpath '/']; end
fname = get(handles.save_basefilename,'String');
fname(strfind(fname,' '))='_'; % no spaces!

button   = get(handles.save_fileformat,'Value');
withcbar = get(handles.save_colorbar,'Value');

temph=figure; copyobj(handles.tf_axes,temph); set(gca,'units','normalized','position',[0.1 0.1 .8 .8])
if withcbar, colorbar; end
if button==1
    print(temph,[fpath fname '_TFmap.png' ],'-dpng');
    ext='png';
else print(temph,[fpath fname '_TFmap.eps' ],'-depsc'); ext='eps';
end; close(temph);

temph=figure; copyobj(handles.topo_axes,temph); set(gca,'units','normalized','position',[0.1 0.1 .8 .8])
if withcbar, colorbar; end
if button==1
    print(temph,[fpath fname '_topomap.png' ],'-dpng');
else print(temph,[fpath fname '_topomap.eps' ],'-depsc');
end; close(temph);

disp([ 'Saved TF plot and topomap to ' fpath fname '_TF*.' ext ]);

% --- Executes on selection change in topo_plottype.
function topo_plottype_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 

% --- Executes during object creation, after setting all properties.
function topo_plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topo_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function topo_numcontour_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 

% --- Executes during object creation, after setting all properties.
function topo_numcontour_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on button press in tf_showcrosshairs.
function tf_showcrosshairs_Callback(hObject, eventdata, handles)
if get(handles.tf_showcrosshairs,'Value')==1
    set(handles.tflineV,'visible','on');
    set(handles.tflineH,'visible','on');
else
    set(handles.tflineV,'visible','off');
    set(handles.tflineH,'visible','off');
end

%% --- Executes on selection change in tf_plottype.
function tf_plottype_Callback(hObject, eventdata, handles)

data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
cla(handles.tf_axes,'reset');

if get(handles.tf_plottype,'Value')==2
    handles.tfobject=imagesc(t,[],squeeze(ctf(get(handles.tf_pickelectrode,'Value'),:,:)),'Parent',handles.tf_axes); axis xy
    yskip=ceil(length(f)/5);
    set(handles.tf_axes,'ydir','normal','ytick',1:yskip:length(f),'yticklabel',round(f(1:yskip:end)*10)/10);
elseif any(get(handles.tf_plottype,'Value')==[1 3])
    [junk,handles.tfobject]=contourf(handles.tf_axes,t,f,squeeze(ctf(get(handles.tf_pickelectrode,'Value'),:,:)),40,'linecolor','none');
end

% update lines
hold(handles.tf_axes,'on')
handles.tflineV=plot(handles.tf_axes,[str2double(get(handles.tf_time2plot,'String')) str2double(get(handles.tf_time2plot,'String'))],get(handles.tf_axes,'ylim'),'k');
handles.tflineH=plot(handles.tf_axes,get(handles.tf_axes,'xlim'),[str2double(get(handles.tf_freq2plot,'String')) str2double(get(handles.tf_freq2plot,'String'))],'k');
if get(handles.tf_showcrosshairs,'Value')==1, linevisibility='on'; else linevisibility='off'; end
set(handles.tflineV,'visible',linevisibility);
set(handles.tflineH,'visible',linevisibility);

% Update handles structure
guidata(handles.figure1,handles);
set(handles.tfobject,'buttondownfcn',{@tf_axes_ButtonDownFcn,handles});

% update axis
update_tfmap(handles);

%% --- Executes during object creation, after setting all properties.
function tf_plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tf_freqscaling.
function tf_freqscaling_Callback(hObject, eventdata, handles)
update_tfmap(handles);
handles=update_topomap(handles); guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function tf_freqscaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_freqscaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on mouse press over axes background.
function tf_axes_ButtonDownFcn(hObject, eventdata, handles)
xy = get(handles.tf_axes,'Currentpoint');
x=xy(1,1); y=xy(1,2);

set(handles.tflineV,'xdata',[x x],'ydata',get(handles.tf_axes,'ylim'),'hittest','off');
set(handles.tflineH,'xdata',get(handles.tf_axes,'xlim'),'ydata',[y y],'hittest','off');

% update TF displays
set(handles.tf_time2plot,'String',num2str(round(x)));
set(handles.tf_freq2plot,'String',num2str(round(y)));

enames=get(handles.tf_pickelectrode,'String');
title(handles.tf_axes,[ 'Electrode ' enames{get(handles.tf_pickelectrode,'Value')} ' (#' num2str(get(handles.tf_pickelectrode,'Value')) ')' ]);

% plot topoplot
handles=update_topomap(handles); guidata(handles.figure1,handles); 

function tf_time2plot_Callback(hObject, eventdata, handles)
update_tfmap(handles);
guidata(handles.figure1,handles); 

% --- Executes during object creation, after setting all properties.
function tf_time2plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_time2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tf_freq2plot_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 
update_tfmap(handles)

% --- Executes during object creation, after setting all properties.
function tf_freq2plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_freq2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function topo_plotrad_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 

% --- Executes during object creation, after setting all properties.
function topo_plotrad_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function topo_clim_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 

% --- Executes during object creation, after setting all properties.
function topo_clim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on button press in topo_showselectedelectrode.
function topo_showselectedelectrode_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 


%%
function handles=update_topomap(handles)
% plot topoplot
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
set(handles.figure1,'CurrentAxes',handles.topo_axes)
if get(handles.topo_showselectedelectrode,'Value')==1, showelectrodes='on'; else showelectrodes='off'; end
toposhading='interp'; if get(handles.topo_plottype,'Value')==2, toposhading='interp'; elseif get(handles.topo_plottype,'Value')==3 toposhading='flat'; end

etype='off';
switch get(handles.topo_electrodetype,'value')
    case 1, etype = 'off';
    case 2, etype = 'on';
    case 3, etype = 'labels';
    case 4, etype = 'numbers';
end

[junk,fpoint] = min(abs(f-str2num(get(handles.tf_freq2plot,'string'))));
[junk,tpoint] = min(abs(t-str2num(get(handles.tf_time2plot,'string'))));
[handles.topoobject,handles.pltchans,handles.epos] = topoplotX(squeeze(ctf(:,fpoint,tpoint)),chanlocs,'plotrad',str2double(get(handles.topo_plotrad,'string')),'numcontour',str2double(get(handles.topo_numcontour,'string')),'electrodes',etype,'shading',toposhading);
% set color limits
topoclim = sscanf(get(handles.topo_clim,'string'),'%g');
if numel(topoclim)==1, topoclim=[-abs(topoclim) abs(topoclim)]; end
set(handles.topo_axes,'clim',topoclim);
% plot selected electrode
e2plot = handles.pltchans==get(handles.tf_pickelectrode,'value');
if any(e2plot)
    handles.selected_electrode = plot3(handles.epos(2,e2plot),handles.epos(1,e2plot),1,'ko','linew',2,'markerface','m','hittest','off');
    if get(handles.topo_showselectedelectrode,'value')==1
        set(handles.selected_electrode,'visible','on');
    else set(handles.selected_electrode,'visible','off');
    end
end

set(handles.topoobject,'buttondownfcn',{@topo_axes_ButtonDownFcn,handles});

title(handles.topo_axes,[ get(handles.tf_time2plot,'string') ' ms, ' get(handles.tf_freq2plot,'string') ' Hz' ]);

% update TFS point
set(handles.datapointvalue,'string',num2str(ctf(get(handles.tf_pickelectrode,'value'),fpoint,tpoint)))

% Update handles structure
guidata(handles.figure1,handles);

%%
function update_tfmap(handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
% get TF point
time2plot = str2double(get(handles.tf_time2plot,'String'));
freq2plot = str2double(get(handles.tf_freq2plot,'String'));
title(handles.tf_axes,[ 'Electrode ' chanlocs(get(handles.tf_pickelectrode,'Value')).labels ' (#' num2str(get(handles.tf_pickelectrode,'Value')) ')' ]);

if time2plot>t(end)
    time2plot=t(end);
    set(handles.tf_time2plot,'string',num2str(time2plot));
end
if freq2plot>f(end)
    freq2plot=f(end);
    set(handles.tf_freq2plot,'string',num2str(freq2plot));
end

% fix out-of-bounds values
if abs(freq2plot-f(end))<.001; freq2plot=freq2plot*.975;  end
if abs(freq2plot-f(1))<.001;   freq2plot=freq2plot*1.05; end
if time2plot>t(end), time2plot=t(end); end
if time2plot<t(1),   time2plot=t(1);   end
if get(handles.tf_plottype,'Value')==2
    [junk,freq2plot]=min(abs(f-freq2plot));
end

% update crosshairs
set(handles.tflineV,'xdata',[time2plot time2plot],'ydata',get(handles.tf_axes,'ylim'),'hittest','off');
set(handles.tflineH,'xdata',get(handles.tf_axes,'xlim'),'ydata',[freq2plot freq2plot],'hittest','off');
tf_showcrosshairs_Callback(7,7,handles);

clim = sscanf(get(handles.tf_clim,'String'),'%g');
if length(clim)==1, clim=[-abs(clim) abs(clim)];end
if isempty(clim), clim=[-3 3]; end
if length(clim)>2,  clim=clim(1:2); end
set(handles.tf_axes,'clim',clim);

yskip=5;
if get(handles.tf_freqscaling,'Value')==1
    if get(handles.tf_plottype,'Value')==2
        set(handles.tf_axes,'ydir','normal','ytick',1:yskip:length(f),'yticklabel',round(f(1:yskip:end)*10)/10);
    else
        set(handles.tf_axes,'yscale','log','ytick',round(logspace(log10(f(1)),log10(f(end)),length(f(1:yskip:end)))*10)/10);
    end
else
    if get(handles.tf_plottype,'Value')==3
        set(handles.tf_axes,'yscale','linear','ytick',round(linspace(f(1),f(end),length(f(1:yskip:end)))*10)/10);
    end
end

% update xlim
xlim = sscanf(get(handles.tf_xlim,'String'),'%g');
if numel(xlim)==3, xlim=xlim(1:2); end
if numel(xlim)==1, xlim(2)=t(end); end
if xlim(1)<t(1), xlim(1)=t(1); end
if xlim(end)>t(end), xlim(end)=t(end); end

set(handles.tf_axes,'xlim',xlim)

%% --- Executes on selection change in tf_pickelectrode.
function tf_pickelectrode_Callback(hObject, eventdata, handles)
tf_plottype_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles); 

% --- Executes during object creation, after setting all properties.
function tf_pickelectrode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function topo_axes_ButtonDownFcn(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
xy = get(handles.topo_axes,'Currentpoint');
x=xy(1,1); y=xy(1,2);

nearest_electrode = dsearchn(handles.epos',[y x]);
ctemp = chanlocs(handles.pltchans);
set(handles.tf_pickelectrode,'value',find(strcmpi(get(handles.tf_pickelectrode,'string'),ctemp(nearest_electrode).labels)));
handles=update_topomap(handles); guidata(handles.figure1,handles); 
tf_plottype_Callback(hObject, eventdata, handles)

% --- Executes on selection change in topo_electrodetype.
function topo_electrodetype_Callback(hObject, eventdata, handles)
handles=update_topomap(handles); guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function topo_electrodetype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% topoplotX() - hack of eeglab's topoplot to work with tfviewerx.
function [handle,pltchans,epos] = topoplotX(Values,chanlocs,varargin)

%% Set defaults

headrad = 0.5;          % actual head radius - Don't change this!
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
HEADCOLOR = [0 0 0];    % default head color (black)
HLINEWIDTH = 1.7;         % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring
HEADRINGWIDTH    = .007;% width of the cartoon head ring

nargs = nargin;
if nargs > 2
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'plotrad'
                plotrad = Value;
            case 'shading'
                SHADING = lower(Value);
                if ~any(strcmp(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
        end
    end
end

Values = Values(:); % make Values a column vector

%% Read channel location
labels={chanlocs.labels};
Th=[chanlocs.theta];
Rd=[chanlocs.radius];

Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);
plotchans = 1:length(chanlocs);

%% remove infinite and NaN values

inds = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
for chani=1:length(chanlocs)
    if isempty(chanlocs(chani).X); inds = [inds chani]; end
end

plotchans   = setdiff(plotchans,inds);

[x,y]       = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans   = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
Th          = Th(plotchans);
Rd          = Rd(plotchans);
x           = x(plotchans);
y           = y(plotchans);
labels      = char(labels(plotchans)); % remove labels for electrodes without locations
Values      = Values(plotchans);
intrad      = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location

%% Find plotting channels
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square

%% Eliminate channels not plotted

allx  = x;
ally  = y;
allchansind = allchansind(pltchans);
intTh = Th(intchans);           % eliminate channels outside the interpolation area
intRd = Rd(intchans);
intx  = x(intchans);
inty  = y(intchans);
Th    = Th(pltchans);              % eliminate channels outside the plotting area
Rd    = Rd(pltchans);
x     = x(pltchans);
y     = y(pltchans);

intValues = Values(intchans);
Values = Values(pltchans);

labels= labels(pltchans,:);

%% Squeeze channel locations to <= headrad
squeezefac = headrad/plotrad;
intRd = intRd*squeezefac; % squeeze electrode arc_lengths towards the vertex
Rd    = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
% to plot all inside the head cartoon
intx  = intx*squeezefac;
inty  = inty*squeezefac;
x     = x*squeezefac;
y     = y*squeezefac;
allx  = allx*squeezefac;
ally  = ally*squeezefac;

%% create grid
xmin = min(-headrad,min(intx)); xmax = max(headrad,max(intx));
ymin = min(-headrad,min(inty)); ymax = max(headrad,max(inty));
xi   = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
yi   = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

[Xi,Yi,Zi] = griddata(inty,intx,intValues,yi',xi,'v4'); % interpolate data

%% Mask out data outside the head
mask = (sqrt(Xi.^2 + Yi.^2) <= headrad); % mask outside the plotting circle
Zi(mask == 0) = NaN;                  % mask non-plotting voxels with NaNs
grid = plotrad;                       % unless 'noplot', then 3rd output arg is plotrad
delta = xi(2)-xi(1); % length of grid entry

%% Scale the axes and make the plot
cla  % clear current axis
hold on
h = gca; % uses current axes
AXHEADFAC = 1.05;     % do not leave room for external ears if head cartoon
set(gca,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC);
unsh = (GRID_SCALE+1)/GRID_SCALE; % un-shrink the effects of 'interp' SHADING

if strcmp(SHADING,'interp')
    handle = surface(Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
else
    handle = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
end
contour(Xi,Yi,Zi,CONTOURNUM,'k','hittest','off');

%% Plot filled ring to mask jagged grid boundary
hwidth = HEADRINGWIDTH;                   % width of head ring
hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius

if strcmp(SHADING,'interp')
    rwidth = BLANKINGRINGWIDTH*1.3;             % width of blanking outer ring
else
    rwidth = BLANKINGRINGWIDTH;         % width of blanking outer ring
end
rin    =  headrad*(1-rwidth/2);              % inner ring radius
if hin>rin
    rin = hin;                              % dont blank inside the head ring
end

circ = linspace(0,2*pi,CIRCGRID);
rx = sin(circ);
ry = cos(circ);
ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];
ringh = patch(ringx,ringy,0.01*ones(size(ringx)),get(gcf,'color'),'edgecolor','none','hittest','off'); hold on

%% Plot cartoon head, ears, nose
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
ringh = patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off'); hold on

% Plot ears and nose
base  = headrad-.0046;
basex = 0.18*headrad;                   % nose width
tip   = 1.15*headrad;
tiphw = .04*headrad;                    % nose tip half width
tipr  = .01*headrad;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % headrad = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = headrad/plotrad;                                          % squeeze the model ears and nose
% by this factor
plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,2*ones(size([basex;tiphw;0;-tiphw;-basex])),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');                 % plot nose
plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')    % plot left ear
plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')   % plot right ear

%% Mark electrode locations
if strcmp(ELECTRODES,'on')   % plot electrodes as spots
    hp2 = plot3(y,x,ones(size(x)),'.','Color',[0 0 0],'markersize',5,'linewidth',.5,'hittest','off');
elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
    text(double(y),double(x),{chanlocs(pltchans).labels},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off')
elseif strcmp(ELECTRODES,'numbers')
    for i = 1:size(labels,1)
        h = text(double(y(i)),double(x(i)),1,int2str(allchansind(i)),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0]);%,'hittest','off')
        set(h,'hittest','off');
    end
end

epos=[x; y];
axis off
axis equal


% --- Executes on slider movement.
function tf_timeslider_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
td=mean(diff(t)); 
if get(handles.tf_timeslider,'value')==1
    set(handles.tf_time2plot,'string',num2str(round(str2double(get(handles.tf_time2plot,'string'))+td)));
else
    set(handles.tf_time2plot,'string',num2str(round(str2double(get(handles.tf_time2plot,'string'))-td)));
end
update_tfmap(handles);
handles=update_topomap(handles); guidata(handles.figure1,handles);

% --- Executes on button press in tf_timeP.
function tf_timeP_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
td=mean(diff(t)); 
set(handles.tf_time2plot,'string',num2str(round(str2double(get(handles.tf_time2plot,'string'))+td)));
update_tfmap(handles);
handles=update_topomap(handles); guidata(handles.figure1,handles);

function tf_timeM_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
td=mean(diff(t)); 
set(handles.tf_time2plot,'string',num2str(round(str2double(get(handles.tf_time2plot,'string'))-td)));
update_tfmap(handles);
handles=update_topomap(handles); guidata(handles.figure1,handles);

% --- Executes on button press in tf_freqP.
function tf_freqP_Callback(hObject,eventdata,handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
[junk,currf]=min(abs(f-str2double(get(handles.tf_freq2plot,'string'))));
currf=min(length(f),currf+1);
set(handles.tf_freq2plot,'string',num2str(round(f(currf)*10)/10));
update_tfmap(handles);
handles=update_topomap(handles); guidata(handles.figure1,handles);

% --- Executes on button press in tf_freqM.
function tf_freqM_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
[junk,currf]=min(abs(f-str2double(get(handles.tf_freq2plot,'string'))));
currf=max(1,currf-1);
set(handles.tf_freq2plot,'string',num2str(round(f(currf)*10)/10));
update_tfmap(handles);
handles=update_topomap(handles); guidata(handles.figure1,handles);



function tf_xlim_Callback(hObject, eventdata, handles)
update_tfmap(handles);

% --- Executes during object creation, after setting all properties.
function tf_xlim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GoToGlobalMax.
function GoToGlobalMax_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;

% find max point and convert to TFE indices
[junk,maxpoint]  = max(ctf(:));
[maxE,maxF,maxT] = ind2sub(size(ctf),maxpoint);

% update
set(handles.tf_freq2plot,'string',num2str(round(f(maxF)*10)/10));
set(handles.tf_time2plot,'string',num2str(round(t(maxT))));
set(handles.tf_pickelectrode,'value',maxE);
handles=update_topomap(handles); guidata(handles.figure1,handles);
tf_plottype_Callback(hObject, eventdata, handles)

% --- Executes on button press in GoToGlobalMin.
function GoToGlobalMin_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;

% find max point and convert to TFE indices
[junk,minpoint]  = min(ctf(:));
[minE,minF,minT] = ind2sub(size(ctf),minpoint);

% update
set(handles.tf_freq2plot,'string',num2str(round(f(minF)*10)/10));
set(handles.tf_time2plot,'string',num2str(round(t(minT))));
set(handles.tf_pickelectrode,'value',minE);
handles=update_topomap(handles); guidata(handles.figure1,handles);
tf_plottype_Callback(hObject, eventdata, handles)

% --- Executes on button press in synchTo.
function synchTo_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
% write out time-frequency-electrode points
tfviewerx_synch_info = [ str2double(get(handles.tf_freq2plot,'string')) str2double(get(handles.tf_time2plot,'string')) get(handles.tf_pickelectrode,'value') ];
assignin('base','tfviewerx_synch_info',tfviewerx_synch_info);

% --- Executes on button press in synchFrom.
function synchFrom_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;
try
    tfviewerx_synch_info = evalin('base','tfviewerx_synch_info');
    % update
    set(handles.tf_freq2plot,'string',num2str(tfviewerx_synch_info(1)));
    set(handles.tf_time2plot,'string',num2str(tfviewerx_synch_info(2)));
    set(handles.tf_pickelectrode,'value',tfviewerx_synch_info(3));
    handles=update_topomap(handles); guidata(handles.figure1,handles);
    tf_plottype_Callback(hObject, eventdata, handles)
catch me;
end


% --- Executes on button press in GoToSensorMax.
function GoToSensorMax_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;

% find max point and convert to TFE indices
ctf=squeeze(ctf(get(handles.tf_pickelectrode,'Value'),:,:));
[junk,maxpoint]  = max(ctf(:));
[maxF,maxT] = ind2sub(size(ctf),maxpoint);

% update
set(handles.tf_freq2plot,'string',num2str(round(f(maxF)*10)/10));
set(handles.tf_time2plot,'string',num2str(round(t(maxT))));
handles=update_topomap(handles); guidata(handles.figure1,handles);
tf_plottype_Callback(hObject, eventdata, handles)


% --- Executes on button press in GoToSensorMin.
function GoToSensorMin_Callback(hObject, eventdata, handles)
data=get(handles.figure1,'UserData'); t=data.t; f=data.f; ctf=data.ctf; chanlocs=data.chanlocs;

% find max point and convert to TFE indices
ctf=squeeze(ctf(get(handles.tf_pickelectrode,'Value'),:,:));
[junk,minpoint]  = min(ctf(:));
[minF,minT] = ind2sub(size(ctf),minpoint);

% update
set(handles.tf_freq2plot,'string',num2str(round(f(minF)*10)/10));
set(handles.tf_time2plot,'string',num2str(round(t(minT))));
handles=update_topomap(handles); guidata(handles.figure1,handles);
tf_plottype_Callback(hObject, eventdata, handles)

%% 



function save_basefilename_Callback(hObject, eventdata, handles)
% hObject    handle to save_basefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_basefilename as text
%        str2double(get(hObject,'String')) returns contents of save_basefilename as a double


% --- Executes during object creation, after setting all properties.
function save_basefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_basefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_selectdir.
function save_selectdir_Callback(hObject, eventdata, handles)

set(handles.save_filepath,'string',uigetdir);


% --- Executes on button press in save_colorbar.
function save_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to save_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_colorbar


% --- Executes on selection change in save_fileformat.
function save_fileformat_Callback(hObject, eventdata, handles)
% hObject    handle to save_fileformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns save_fileformat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from save_fileformat


% --- Executes during object creation, after setting all properties.
function save_fileformat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_fileformat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function save_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to save_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_filepath as text
%        str2double(get(hObject,'String')) returns contents of save_filepath as a double


% --- Executes during object creation, after setting all properties.
function save_filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
