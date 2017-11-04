function [iChannels, ChannelNames, WrnMsg, ErrMsg] = sc_bst_channel_find(Channel, target)

% Initialization
ErrMsg = [];
WrnMsg = [];
iChannels = [];
ChannelNames = [];

% Get channel selection; make it a cell array of strings
if iscell(target)
    channelSelect = target;
else
    if isempty(target)
        channelSelect = {''};
    else
        channelSelect = textscan(target,'%s','delimiter',',;');
        channelSelect = channelSelect{1};
    end
end

% Channel names in current file
fileChanNames = {Channel.Name};
% Channel types
fileChanTypes = {Channel.Type};
iType = find(ismember(lower(channelSelect),lower(fileChanTypes)));
if isempty(iType)
    hasTypes = false;
else
    % Don't want to flag a type if channel type is also channel name
    iSameNameType = find(strcmpi(fileChanTypes,fileChanNames));
    hasTypes = any(~ismember(lower(channelSelect(iType)),lower(fileChanNames(iSameNameType))));
end
% Get selected channel names and indexes
if numel(channelSelect)<=1 || hasTypes
    % Use Brainstorm: sorted, unique, only not-bad and available
    iChannels = channel_find(Channel, target); % Cautious with CHANNEL_FIND; it returns a sorted list (may not be the order of how it is specified)
    channelNames = fileChanNames(iChannels);
    if isempty(iChannels)
        ErrMsg = sprintf('[%s] are either all marked as bad or unavailable',target);
        return;
    end
else
    % Index of channels specified not present in file
    iNotFound = find(~ismember(lower(channelSelect),lower(fileChanNames)));
    if ~isempty(iNotFound)
        ErrMsg = 'Following channels are unavailable: ';
        ErrMsg = [ErrMsg,sprintf('\n%s',channelSelect{iNotFound})];
        return;
    end
    channelNames = channelSelect;
    % Set indexes of selected channels
    iChannels = zeros(1,numel(channelNames));
    for i=1:numel(channelNames)
        iChannels(i) = find(ismember(lower(fileChanNames),lower(channelNames(i))));
    end
end

end