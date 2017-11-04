function varargout = dis_selchan(varargin)
% DIS_SELCHAN M-file for dis_selchan.fig
%      DIS_SELCHAN, by itself, creates a new DIS_SELCHAN or raises the existing
%      singleton*.
%
%      H = DIS_SELCHAN returns the handle to a new DIS_SELCHAN or the handle to
%      the existing singleton*.
%
%      DIS_SELCHAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIS_SELCHAN.M with the given input arguments.
%
%      DIS_SELCHAN('Property','Value',...) creates a new DIS_SELCHAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dis_selchan_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dis_selchan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: dis_selchan.m 417 2012-01-25 10:25:12Z jessica $

% Edit the above text to modify the response to help dis_selchan

% Last Modified by GUIDE v2.5 13-May-2010 23:11:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dis_selchan_OpeningFcn, ...
    'gui_OutputFcn',  @dis_selchan_OutputFcn, ...
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


% --- Executes just before dis_selchan is made visible.
function dis_selchan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dis_selchan (see VARARGIN)

% Choose default command line output for dis_selchan
set(0,'CurrentFigure',handles.figure1);
handles.output = hObject;

load CRC_electrodes.mat;
handles.names     = names;
handles.pos       = pos';
handles.crc_types = crc_types;

if isempty(varargin) || ~isfield(varargin{1},'file')
    % Filter for vhdr, mat and edf files
    try and(isfield(varargin{1},'multcomp'), varargin{1}.multcomp)
        %if the multiple comparison option is checked
        prefile = spm_select(Inf, 'any', 'Select imported EEG file','' ...
            ,pwd,'\.[mMvVeErR][dDhHaA][fFDdTtwW]');
        set(handles.Selectall,'enable','off','visible','off');
        set(handles.selhalf,'enable','off','visible','off');
        set(handles.otherhalf,'enable','off','visible','off');
        set(handles.selquarter,'enable','off','visible','off');
        set(handles.Quarter2,'enable','off','visible','off');
        set(handles.quarter3,'enable','off','visible','off');
        set(handles.quarter4,'enable','off','visible','off');
        set(handles.desall,'enable','off','visible','off');
        set(handles.checkcICA,'enable','off','visible','off');
        set(handles.Score,'enable','off','visible','off');
        set(handles.text4,'visible','off');
        set(handles.text5,'visible','off');
        handles.multcomp=1;
    catch
        prefile = spm_select(1, 'any', 'Select imported EEG file','' ...
            ,pwd,'\.[mMvVeErR][dDhHaA][fFDdTtwW]');
        handles.multcomp=0;
    end
    for i=1:size(prefile,1)
        D{i} = crc_eeg_load(deblank(prefile(i,:)));
        file = fullfile(D{i}.path,D{i}.fname);
        handles.file{i} = file;
        handles.chan{i} = upper(chanlabels(D{i}));
        if isfield(D{i}, 'info')
            try
                D{i}.info.date;
            catch
                D{i}.info.date = [1 1 1];
            end
            try
                D{i}.info.hour;
            catch
                D{i}.info.hour = [0 0 0];
            end
        else
            D{i}.info = struct('date',[1 1 1],'hour',[0 0 0]);
        end
        handles.Dmeg{i} = D{i};
        handles.date{i} = zeros(1,2);
        handles.date{i}(1) = datenum([D{i}.info.date D{i}.info.hour]);
        handles.date{i}(2) = handles.date{i}(1) + ...
                        datenum([ 0 0 0 crc_time_converts(nsamples(D{i})/ ...
                                                            fsample(D{i}))] );
        handles.dates(i,:) = handles.date{i}(:);
    end

else
    handles.file = varargin{1}.file;
    prefile = deblank(handles.file);
    index = varargin{1}.index;
    for i=1:size(varargin{1}.Dmeg,2)
%         handles.Dmeg{i} = crc_eeg_load([path(varargin{1}.Dmeg{i}),filesep,fname(varargin{1}.Dmeg{i})]);
        handles.Dmeg{i} = varargin{1}.Dmeg{i};
    end
    if isempty(index)
        index=1:nchannels(handles.Dmeg{1});
    end
    set(handles.Select,'String',upper(chanlabels(handles.Dmeg{1},varargin{1}.index)));
    diff = setdiff(upper(chanlabels(handles.Dmeg{1})),upper(chanlabels(handles.Dmeg{1},varargin{1}.index)));
    set(handles.Deselect,'String',diff);

    [dumb1,dumb2,index2]=intersect(upper(chanlabels(handles.Dmeg{1},varargin{1}.index)),upper(handles.names));

    idxred=index(find(handles.crc_types(index2)<-1));
    idxblue=index(find(handles.crc_types(index2)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])
    handles.chan{1}=get(handles.Select,'String');
end

if ~isempty(varargin) && isfield(varargin{1},'delmap')
    handles.delmap=varargin{1}.delmap;
end

if (~isempty(varargin) && isfield(varargin{1},'multcomp') && varargin{1}.multcomp) || ...
        isempty(varargin)
    chanset=handles.chan{1};
    for i=1:size(prefile,1)
        chanset=intersect(chanset,handles.chan{i});
    end
    set(handles.Deselect,'String',chanset);
    handles.chan=chanset;
end

try
    handles.Dmeg{1}.CRC.cICA;
catch
    set(handles.checkcICA,'Visible','off')
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dis_selchan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dis_selchan_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Select.
function Select_Callback(hObject, eventdata, handles)
% hObject    handle to Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Select
contents = get(hObject,'String');
if isempty(contents)
else
    %Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.Select,'String',temp);

    [dumb1,dumb2,index2]=intersect(upper(temp),upper(handles.names));

    idxred=index2(find(handles.crc_types(index2)<-1));
    idxblue=index2(find(handles.crc_types(index2)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])

    set(handles.Localizer,'XTick',[]);
    set(handles.Localizer,'YTick',[]);

    %Add the "activated" in the list "Selected Channels"
    if length(get(handles.Deselect,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.Deselect,'String')];
    end
    set(handles.Deselect,'String',temp);

    %Prevent crashing if the first/last item of the list is selected.
    set(handles.Select,'Value',max(index-1,1));
    set(handles.Deselect,'Value',1);

    %if multiple comparison, no more than one channel
    contents=get(handles.Select,'String');
    if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))
        set(handles.PLOT,'enable','off');
        set(handles.PLOT,'ForegroundColor',[1 0 0]);
        beep
        disp('Select only one channel for multiple files comparison')
    elseif isfield(handles,'multcomp') && (handles.multcomp && length(contents)==1)
        set(handles.PLOT,'enable','on');
        set(handles.PLOT,'ForegroundColor',[0 0 0]);
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Deselect.
function Deselect_Callback(hObject, eventdata, handles)
% hObject    handle to Deselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
if length(contents)==0
else
    %Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.Deselect,'String',temp);

    %Add the "activated" in the list "Selected Channels"
    if length(get(handles.Select,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.Select,'String')];
    end
    set(handles.Select,'String',temp);

    %Prevent crashing if the first/last item of the list is selected.
    set(handles.Deselect,'Value',max(index-1,1));
    set(handles.Select,'Value',1);

    [dumb1,dumb2,index]=intersect(upper(temp),upper(handles.names));

    idxred=index(find(handles.crc_types(index)<-1));
    idxblue=index(find(handles.crc_types(index)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])

    set(handles.Localizer,'XTick',[]);
    set(handles.Localizer,'YTick',[]);

    %if multiple comparison, no more than one channel
    contents=get(handles.Select,'String');
    if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))
        set(handles.PLOT,'enable','off');
        set(handles.PLOT,'ForegroundColor',[1 0 0]);
        beep
        disp('Select only one channel for multiple files comparison')
    elseif isfield(handles,'multcomp') && (handles.multcomp && length(contents)==1)
        set(handles.PLOT,'enable','on');
        set(handles.PLOT,'ForegroundColor',[0 0 0]);
    end

end
% Update handles structure
guidata(hObject, handles);

% Hints: contents = get(hObject,'String') returns Deselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Deselect

% --- Executes during object creation, after setting all properties.
function Deselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Deselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Selectall.
function Selectall_Callback(hObject, eventdata, handles)
% hObject    handle to Selectall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Select,'String',upper(chanlabels(handles.Dmeg{1})));
set(handles.Deselect,'String',cell(0));

[dumb1,dumb2,index]=intersect(upper(chanlabels(handles.Dmeg{1})),upper(handles.names));

idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in selhalf.
function selhalf_Callback(hObject, eventdata, handles)
% hObject    handle to selhalf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
total=1:1:length(chanlabels(handles.Dmeg{1}));
half=1:1:length(chanlabels(handles.Dmeg{1}))/2;
otherpart=setdiff(total,half);
half=upper(chanlabels(handles.Dmeg{1},half));
otherpart=upper(chanlabels(handles.Dmeg{1},otherpart));

set(handles.Deselect,'String',otherpart);
set(handles.Select,'String',half);

[dumb1,dumb2,index]=intersect(upper(half),upper(handles.names));
idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in selquarter.
function selquarter_Callback(hObject, eventdata, handles)
% hObject    handle to selquarter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
total=1:1:length(chanlabels(handles.Dmeg{1}));
quarter=1:1:length(chanlabels(handles.Dmeg{1}))/4;
otherpart=setdiff(total,quarter);
quarter=upper(chanlabels(handles.Dmeg{1},quarter));
otherpart=upper(chanlabels(handles.Dmeg{1},otherpart));

set(handles.Deselect,'String',otherpart);
set(handles.Select,'String',quarter);

[dumb1,dumb2,index]=intersect(upper(quarter),upper(handles.names));
idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in desall.
function desall_Callback(hObject, eventdata, handles)
% hObject    handle to desall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Deselect,'String',upper(chanlabels(handles.Dmeg{1})));
set(handles.Select,'String',cell(0));

cleargraph(handles)
xlim([0 1])
ylim([0 1])
set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);

% --- Executes on button press in PLOT
function PLOT_Callback(hObject, eventdata, handles)
% hObject    handle to PLOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents=get(handles.Select,'String');
if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))%no more than one channel if multiple file comparison
    beep
    disp('Select only one channel for multiple files comparison')
    return
elseif ~isfield(handles,'multcomp') || ~handles.multcomp
    [dumb1,index]=intersect(upper(chanlabels(handles.Dmeg{1})),upper(get(handles.Select,'String')));
    try
        flags.index=sortch(handles.Dmeg{1},index);
    catch
        flags.index=fliplr(sort(index));
    end
else
    flags.dates=handles.dates;
    flags.chanset=handles.chan;
    [dumb1,index]=intersect(upper(handles.chan),upper(get(handles.Select,'String')));
    flags.index=index;
    flags.multcomp=1;
end
flags.Dmeg=handles.Dmeg;
flags.file=handles.file;
if isfield(handles,'delmap')
    flags.delmap=handles.delmap;
end
crc_dis_main(flags);
delete(handles.figure1)

% --- Executes on button press in Quarter2.
function Quarter2_Callback(hObject, eventdata, handles)
% hObject    handle to Quarter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
total=1:1:length(chanlabels(handles.Dmeg{1}));
quarter=round(length(chanlabels(handles.Dmeg{1}))/4+1):1:round(2*length(chanlabels(handles.Dmeg{1}))/4);
otherpart=setdiff(total,quarter);
quarter=upper(chanlabels(handles.Dmeg{1},quarter));
otherpart=upper(chanlabels(handles.Dmeg{1},otherpart));

set(handles.Deselect,'String',otherpart);
set(handles.Select,'String',quarter);

[dumb1,dumb2,index]=intersect(upper(quarter),upper(handles.names));
idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in quarter3.
function quarter3_Callback(hObject, eventdata, handles)
% hObject    handle to quarter3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
total=1:1:length(chanlabels(handles.Dmeg{1}));
quarter=round(2*length(chanlabels(handles.Dmeg{1}))/4+1):1:round(3*length(chanlabels(handles.Dmeg{1}))/4);
otherpart=setdiff(total,quarter);
quarter=upper(chanlabels(handles.Dmeg{1},quarter));
otherpart=upper(chanlabels(handles.Dmeg{1},otherpart));

set(handles.Deselect,'String',otherpart);
set(handles.Select,'String',quarter);

[dumb1,dumb2,index]=intersect(upper(quarter),upper(handles.names));
idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in quarter4.
function quarter4_Callback(hObject, eventdata, handles)
% hObject    handle to quarter4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
total=1:1:length(chanlabels(handles.Dmeg{1}));
quarter=round(3*length(chanlabels(handles.Dmeg{1}))/4)+1:1:length(chanlabels(handles.Dmeg{1}));
otherpart=setdiff(total,quarter);
quarter=upper(chanlabels(handles.Dmeg{1},quarter));
otherpart=upper(chanlabels(handles.Dmeg{1},otherpart));

set(handles.Deselect,'String',otherpart);
set(handles.Select,'String',quarter);

[dumb1,dumb2,index]=intersect(upper(quarter),upper(handles.names));
idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in otherhalf.
function otherhalf_Callback(hObject, eventdata, handles)
% hObject    handle to otherhalf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

total=1:1:length(chanlabels(handles.Dmeg{1}));
half=round(length(chanlabels(handles.Dmeg{1}))/2+1):1:length(chanlabels(handles.Dmeg{1}));
otherpart=setdiff(total,half);
half=upper(chanlabels(handles.Dmeg{1},half));
otherpart=upper(chanlabels(handles.Dmeg{1},otherpart));

set(handles.Deselect,'String',otherpart);
set(handles.Select,'String',half);

[dumb1,dumb2,index]=intersect(upper(half),upper(handles.names));
idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(find(handles.crc_types(index)>-2));

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.Localizer,'XTick',[]);
set(handles.Localizer,'YTick',[]);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Score.
function Score_Callback(hObject, eventdata, handles)
% hObject    handle to Score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[dumb1,index,dumb2]=intersect(upper(chanlabels(handles.Dmeg{1})),get(handles.Select,'String'));
try
    flags.index=sortch(handles.Dmeg{1},index);
catch
    flags.index=fliplr(sort(index));
end

flags.Dmeg=handles.Dmeg;
flags.file=handles.file;
flags.scoresleep=1;
if isfield(handles,'delmap')
    flags.delmap=handles.delmap;
end
crc_dis_main(flags);
delete(handles.figure1)

% --- Executes on button press in checkcICA.
function checkcICA_Callback(hObject, eventdata, handles)
% hObject    handle to checkcICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dumb1,index,dumb2] = intersect(upper(chanlabels(handles.Dmeg{1})), ...
                                get(handles.Select,'String'));
try
    flags.index=sortch(handles.Dmeg{1},index);
catch
    flags.index=fliplr(sort(index))  ;
end
flags.Dmeg = handles.Dmeg;
flags.file = handles.file;
dis_cICAcheck(flags);
delete(handles.figure1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SORT CHANNEL FX
function idx=sortch(d,index)

% Find the EOG channels

A=strfind(upper(chanlabels(d,index)),'EOG');
EOGchan=[];
for i=1:length(A);
    if A{i}>0
        EOGchan=[EOGchan index(i)];
    end
end
% Find the ECG channels

A=strfind(upper(chanlabels(d,index)),'ECG');
ECGchan=[];
for i=1:length(A);
    if A{i}>0
        ECGchan=[ECGchan index(i)];
    end
end

A=strfind(upper(chanlabels(d,index)),'EMG');
EMGchan=[];
for i=1:length(A);
    if A{i}>0
        EMGchan=[EMGchan index(i)];
    end
end

A=strfind(upper(chantype(d)),'EEG');
EEGchan=[];
for i=1:length(A);
    if A{i}>0
        EEGchan=[EEGchan index(i)];
    end
end

allbad=[EOGchan ECGchan EMGchan EEGchan];
other=setdiff(index,allbad);

otherknown=[];
othernotknown=[];
chanstr=channels(d);
for ff = other
    if chanstr.order(ff)==0
        othernotknown = [othernotknown ff];
    else
        otherknown = [otherknown ff];
    end
end
other = [otherknown othernotknown];


allbad=[EOGchan ECGchan EMGchan other];
eeg=setdiff(index,allbad);

AFrontal = intersect(find(strncmp(chanlabels(d),'AF',2) ==1),eeg);
Frontal = intersect(find(strncmp(chanlabels(d),'F',1) ==1),eeg);
Coronal = intersect(find(strncmp(chanlabels(d),'C',1) ==1),eeg);
Temporal = intersect(find(strncmp(chanlabels(d),'T',1) ==1),eeg);
Parietal= intersect(find(strncmp(chanlabels(d),'P',1) ==1),eeg);
Occipital= intersect(find(strncmp(chanlabels(d),'O',1) ==1),eeg);

neweeg = [Occipital Parietal Temporal Coronal Frontal AFrontal];
eeg2=setdiff(eeg,neweeg);

eeg = [eeg2 neweeg];

idx=[otherknown othernotknown ECGchan EMGchan EOGchan eeg];
%%%%%%%%%%%%%%%
function cleargraph(handles)

A=get(handles.figure1,'Children');

idx=find(strcmp(get(A,'Type'),'axes')==1);

delete(get(A(idx),'Children'))

