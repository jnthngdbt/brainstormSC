% This code was copied from <figure_connect.m> (line 2438) (20130830). It
% corresponds to method <UpdateColormap>. It has been modified this way:
% 
% 1. Optional input "CLim" that specifies wanted colormap limit
% 2. Added a condition in the definition of "CLim" (use the input if
% present, otherwise do what BST was doing
% 3. Adapt calling of methods to be external
% 
% UNFORTUNATELY, due to active development of <figure_connect> function,
% future Brainstorm updates may result in bug in this present function (we
% recall that we copy-pasted code here). See the last line of the function
% below for an example of bug that appeared (OGL.requestRedraw() was
% renamed to OGL.repaint()).


% %% ===== UPDATE COLORMAP =====
% function UpdateColormap(hFig,CLim)    %%%%%
function sc_bst_fig_connect_clim(hFig,CLim)    %%%%%
% If the connectivity plot that you have to edit is the current figure,
% send "gcf" as HFIG. CLIM is a 2 element vector ([MIN, MAX], the value 
% limits)

    % Get selected frequencies and rows
    TfInfo = getappdata(hFig, 'Timefreq');
    if isempty(TfInfo)
        return
    end
    % Get data description
    iDS = bst_memory('GetDataSetTimefreq', TfInfo.FileName);
    if isempty(iDS)
        return
    end
    % Get colormap
    sColormap = bst_colormaps('GetColormap', hFig);
    % Get DataPair
    [DataPair, DataMask] = figure_connect('GetPairs',hFig);    %%%%%
    if sColormap.isAbsoluteValues
        DataPair(:,3) = abs(DataPair(:,3));
    end
    % Get figure method
    Method = getappdata(hFig, 'Method');
    % Get maximum values
    DataMinMax = getappdata(hFig, 'DataMinMax');
    % Get threshold min/max values
    ThresholdMinMax = getappdata(hFig, 'ThresholdMinMax');
    % === COLORMAP LIMITS ===
    if ~exist('CLim','var') %%%%%%%
        % Get units type
        if ismember(Method, {'granger', 'plv', 'plvt'})
            UnitsType = 'timefreq';
            CLim = [DataMinMax(1) DataMinMax(2)];
        elseif ismember(Method, {'corr'})
            UnitsType = 'connect';
            if sColormap.isNormalized
                CLim = ThresholdMinMax;
                if sColormap.isAbsoluteValues
                    CLim = abs(CLim);            
                end
            else
                if sColormap.isAbsoluteValues
                    CLim = [0, 1];
                else
                    CLim = [-1, 1];
                end
            end
        elseif ismember(Method, {'cohere'})
            UnitsType = 'connect';
            CLim = [0, 1];
        end
    else %%%%%%%%%%%
        UnitsType = 'other'; %%%%%%%%%%%%%
    end %%%%%%%%%%%
    setappdata(hFig, 'CLim', CLim);
    
    % === SET COLORMAP ===
    % Update colorbar font size
    hColorbar = findobj(hFig, '-depth', 1, 'Tag', 'Colorbar');
    if ~isempty(hColorbar)
        set(hColorbar, 'FontSize', bst_get('FigFont'), 'FontUnits', 'points');
    end
    % Get figure colormap
    ColormapInfo = getappdata(hFig, 'Colormap');
    sColormap = bst_colormaps('GetColormap', ColormapInfo.Type);
    % Set figure colormap
    set(hFig, 'Colormap', sColormap.CMap);
    % Create/Delete colorbar
    bst_colormaps('SetColorbarVisible', hFig, sColormap.DisplayColorbar);
    % Display only one colorbar (preferentially the results colorbar)
    bst_colormaps('ConfigureColorbar', hFig, ColormapInfo.Type, UnitsType);
    
    % === UPDATE DISPLAY ===
    CMap = sColormap.CMap;
    OGL = getappdata(hFig, 'OpenGLDisplay');
    is3DDisplay = getappdata(hFig, 'is3DDisplay');
    
    if (sum(DataMask) > 0)
        % Normalize DataPair for Offset
        Max = max(DataPair(:,3));
        Min = min(abs(DataPair(:,3)));
        Diff = (Max - Min);
        if (Diff == 0)
            Offset = DataPair(DataMask,3);
        else
            Offset = (abs(DataPair(DataMask,3)) - Min) ./ (Max - Min);
        end
        % Interpolate
        [StartColor, EndColor] = figure_connect('InterpolateColorMap',hFig, DataPair(DataMask,:), CMap, CLim);  %%%%%%%
        % Update color
        OGL.setMeasureLinkColorGradient( ...
            find(DataMask) - 1, ...
            StartColor(:,1), StartColor(:,2), StartColor(:,3), ...
            EndColor(:,1), EndColor(:,2), EndColor(:,3));
        if (~is3DDisplay)
            % Offset is always in absolute
            OGL.setMeasureLinkOffset(find(DataMask) - 1, Offset(:).^2 * 2);
        end
    end
    
    [RegionDataPair, RegionDataMask] = figure_connect('GetRegionPairs',hFig); %%%%%%%
    if (sum(RegionDataMask) > 0)
        % Normalize DataPair for Offset
        Max = max(RegionDataPair(:,3));
        Min = min(RegionDataPair(:,3));
        Diff = (Max - Min);
        if (Diff == 0)
            Offset = RegionDataPair(RegionDataMask,3);
        else
            Offset = (abs(RegionDataPair(RegionDataMask,3)) - Min) ./ (Max - Min);
        end
        % Normalize within the colormap range 
        [StartColor, EndColor] = figure_connect('InterpolateColorMap',hFig, RegionDataPair(RegionDataMask,:), CMap, CLim); %%%%%%
        % Update display
        OGL.setRegionLinkColorGradient( ...
            find(RegionDataMask) - 1, ...
            StartColor(:,1), StartColor(:,2), StartColor(:,3), ...
            EndColor(:,1), EndColor(:,2), EndColor(:,3));
        if (~is3DDisplay)
            % Offset is always in absolute
            OGL.setRegionLinkOffset(find(RegionDataMask) - 1, Offset(:).^2 * 2);
        end
    end
    
%     OGL.requestRedraw(); %%%%%% THIS FUNCTION DOES NOT EXIST ANYMORE
    OGL.repaint(); %%%%% Replaced by this one
end
