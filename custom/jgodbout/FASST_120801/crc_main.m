function varargout = crc_main(varargin)
%__________________________________________________________________________
%     ________  _____________ ___ 
%    / __ _/  |/ ___/_/_  __/|__ \
%   / /_ / /| |\__ \ \ / /   __/ /
%  / __// ___ |__/ / // / _ / __/ 
% /_/  /_/  |_|___/_//_/ (_)____/ 
%
% fMRI Artefact removal and Sleep Scoring Toolbox, FASST.2
% http://www.montefiore.ulg.ac.be/~phillips/FASST.html
%__________________________________________________________________________
%
% CRC_MAIN M-file for crc_main.fig
%      CRC_MAIN, by itself, creates a new CRC_MAIN or raises the existing
%      singleton.
%
%      H = CRC_MAIN returns the handle to a new CRC_MAIN or the handle to
%      the existing singleton*.
%
%      CRC_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRC_MAIN.M with the given input arguments.
%
%      CRC_MAIN('Property','Value',...) creates a new CRC_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crc_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_main.m 462 2012-07-19 10:30:44Z christophe $

% Edit the above text to modify the response to help crc_main

% Last Modified by GUIDE v2.5 15-Oct-2009 11:43:00

% Display ASCII Welcome
disp('                                                            ');
disp('     ________  _____________ ___                            ');
disp('    / ____/  |/ ___/_/_  __/|__ \                           ');
disp('   / /_ / /| |\__ \ \ / /   __/ /                           ');
disp('  / __// ___ |__/ / // / _ / __/                            ');
disp(' /_/  /_/  |_|___/_//_/ (_)____/                            ');
disp('                                                            ');
disp(' fMRI Artefact removal and Sleep Scoring Toolbox, FASST.2   ');
disp(' http://www.montefiore.ulg.ac.be/~phillips/FASST.html       ');
disp(' An SPM8-compatible toolbox.                                ');
fprintf('\n');

% Check if SPM is available, and maybe more one day...
ok = check_installation;
if ~ok
    beep
    fprintf('INSTALLATION PROBLEM!');
    return
end

% Check for Signal Processing Toolbox
persistent flag_TBX
if isempty(flag_TBX)
    flag_TBX = license('checkout','signal_toolbox');
    if ~flag_TBX
        pth = fullfile(spm_str_manip(mfilename('fullpath'),'h'),'SPTfunctions');
        addpath(pth)
        disp(['warning: using freely distributed equivalent to filtering functions ', ...
          'as Signal Processing Toolbox is not available.']);
    end
end


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @crc_main_OpeningFcn, ...
    'gui_OutputFcn',  @crc_main_OutputFcn, ...
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

% --- Executes just before crc_main is made visible.
function crc_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crc_main (see VARARGIN)


[A] = imread('LOGO_Simple.png','BackgroundColor',0.94*[1 1 1]);
image(A)
axis off

% Choose default command line output for crc_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crc_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = crc_main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in push_dis_main.
function push_dis_main_Callback(hObject, eventdata, handles)
% hObject    handle to push_dis_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%delete(handles.figure1)
dis_selchan;

% --- Executes on button press in push_dis_cmp.
function push_dis_cmp_Callback(hObject, eventdata, handles)
% hObject    handle to push_dis_cmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% flags.multcomp=1;
% setappdata(hObject,'multcomp',1);
flags.multcomp=1;
dis_selchan(flags);

% --- Executes on button press in push_credits.
function push_credits_Callback(hObject, eventdata, handles)
% hObject    handle to push_credits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = figure;
set(fig,'Position',[73   145   398   503])
set(fig,'NumberTitle','off')
set(fig,'Name','License & Copyright')
h = axes('Position',[0 0 1 1],'Visible','off');

str{1}= 'Please refer to this version as "FASST, fMRI Artefact Removal and';
str{2}= 'Sleep scoring Toolbox" in papers and communication';
str{3}= ' ';
str{4}= '_____________________________________________________________________';
str{5}= ' ';
str{6}= '';
str{7}= 'For bug reports, contact directly';
str{8}= '';
str{9}= 'Jessica Schrouff, jschrouff@doct.ulg.ac.be';
str{10}='Christophe Phillips, c.phillips@ulg.ac.be';
str{11}=' ';
str{12}='Feel free to submit add-ons to this toolbox';
str{13}='but please note that support will not be provided for';
str{14}='"home-made" changes of the distributed code.';
str{15}=' ';
str{16}='More details should (soon) be available on';
str{17}='             http://www.montefiore.ulg.ac.be/~phillips/FASST.html';
str{18}='___________________________________________________________________';
str{19}=' ';
str{20}='FASST is developed by the Cyclotron Research Centre, part of';
str{21}='the University of Liege (ULg), BE. This work is supported by the ';
str{22}='FRS-FNRS, the Queen Elizabeth''s funding and the University of Liege.';
str{23}=' ';
str{24}='___________________________________________________________________';
str{25}=' ';
str{26}='FASST (being the collection of files given in the manifest in the';
str{27}='Contents.m file) is free but copyright software, distributed under the';
str{28}='terms of the GNU General Public Licence as published by the Free Software';
str{29}='Fundation (either version 2, as given in file CRC_LICENCE.man or at your';
str{30}='option, any later version). Further details on "copyleft" can be found at';
str{31}='http://www.gnu.org/copyleft/';
str{32}=' ';
str{33}='___________________________________________________________________';
str{34}='Copyright (C) 2010 Cyclotron Research Centre, University of Liege';
str{35}=' ';
text(.025,.5,str,'FontSize',8)

% --- Executes on button press in push_concatenate.
function push_concatenate_Callback(hObject, eventdata, handles)
% hObject    handle to push_concatenate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_concatenate;

% --- Executes on button press in push_disfrqcomp.
function push_disfrqcomp_Callback(hObject, eventdata, handles)
% hObject    handle to push_disfrqcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frqcomp;

% --- Executes on button press in push_freqplot.
function push_freqplot_Callback(hObject, eventdata, handles)
% hObject    handle to push_freqplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frq;

% --- Executes on button press in push_crcgar.
function push_crcgar_Callback(hObject, eventdata, handles)
% hObject    handle to push_crcgar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_gar;

% --- Executes on button press in push_crcpar.
function push_crcpar_Callback(hObject, eventdata, handles)
% hObject    handle to push_crcpar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_par;

% --- Executes on button press in push_score.
function push_score_Callback(hObject, eventdata, handles)
% hObject    handle to push_score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_selchan;

% --- Executes on button press in push_freqplotstat.
function push_freqplotstat_Callback(hObject, eventdata, handles)
% hObject    handle to push_freqplotstat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frq;

% --- Executes on button press in push_chunk.
function push_chunk_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_chunks;

% --- Executes on button press in sws.
function sws_Callback(hObject, eventdata, handles)
% hObject    handle to sws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_wave_detection;

function push_batch_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_batch;

%% SUBFUNCTION

function ok = check_installation
% function to check installation state of toolbox,
% particullarly the SPM path setup

ok=1;

if exist('spm.m','file')
    if ~strcmpi(spm('ver'),'spm8')
        beep
        fprintf('\nERROR:\n')
        fprintf('\tSPM8 should be installed on your computer, and\n')
        fprintf('\tbe available on MATLABPATH!\n\n')
        ok = 0;
    end
else
    beep
    fprintf('\nERROR:\n')
    fprintf('\tSPM8 should be installed on your computer, and\n')
    fprintf ('\tbe available on MATLABPATH!\n\n')
    ok = 0;
end

return


