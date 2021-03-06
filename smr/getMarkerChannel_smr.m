function[data,h]=getMarkerChannel_smr(fID, chan, varargin)
% getMarkerChannel_smrCHANNEL reads a marker channel from a file.
%
% [data{, h}]=getMarkerChannel_smr(fID, CHAN)
% fID is the MATLAB file handle and CHAN is the channel number (1 to Max)
% DATA is a structure containing:
%   DATA.TIMINGS: a length n vector with the marker timestamps
%   DATA.MARKERS: an n x 4 array of uint8 type, containing the marker
%   values
%
% When present, OPTIONS must be the last input argument. Valid options
% are:
% 'ticks', 'microseconds', 'milliseconds' and 'seconds' cause times to
%    be scaled to the appropriate unit (seconds by default)in HEADER
% 'scale' - no effect
% 'progress' - causes a progress bar to be displayed during the read.

Info=getInfo_smr(fID,chan);
if isempty (Info)
    data=[];
    h=[];
    return;
end;
if(Info.kind ~= 5)
    warning('getMarkerChannel_smrChannel: Channel #%d No data or not a marker channel', chan);
    data=[];
    h=[];
    return;
end;


FileH=getHeader_smr(fID);
SizeOfHeader=20;                                            % Block header is 20 bytes long
header=getBlockHeaders_smr(fID,chan);

if isempty(header)
    warning('getRealMarkerChannel_smr: Channel number #%d has blank header',chan);
    data=[];
    h=[];
    return;
else
    
    ShowProgress=0;
    arguments=nargin;
    for i=1:length(varargin)
        if ischar(varargin{i})
            arguments=arguments-1;
            if strcmpi(varargin{i},'progress') && Info.blocks>10
                ShowProgress=1;
                progbar=progressbar(0,sprintf('Analyzing %d blocks on channel %d',Info.blocks,chan),...
                    'Name',sprintf('%s',fopen(fID)));
            end;
        end;
    end;
    
    switch arguments
        case {2}
            startBlock=1;
            endBlock=Info.blocks;
        case {3}
            startBlock=varargin{1};
            endBlock=varargin{1};
        otherwise
            startBlock=varargin{1};
            endBlock=min(Info.blocks,varargin{2});
    end;
    
    NumberOfMarkers=sum(header(5,startBlock:endBlock)); % Sum of samples in required blocks
    
    data.timings=zeros(NumberOfMarkers,1);
    data.markers=uint8(zeros(NumberOfMarkers,4));
    
    count=1;
    for block=startBlock:endBlock
        fseek(fID, header(1, block)+SizeOfHeader, 'bof');                     % Start of block
        for i=1:header(5,block)                                              % loop for each marker
            data.timings(count)=fread(fID,1,'int32');                    % Time
            data.markers(count,:)=fread(fID,4,'uint8=>uint8');                    % 4x marker bytes
            count=count+1;
            if ShowProgress==1
                done=(i-startBlock)/max(1,endBlock-startBlock);
                progressbar(done, progbar,sprintf('Reading Channel %d....     %d%% Done',chan,(int16(done*100)/5)*5));
            end;
        end;
    end
    
    
    if(nargout>1)
        h.FileName=Info.FileName;                                   % Set up the header information to return
        h.system=['Matoff' num2str(FileH.systemID)];                   % if it's been requested
        h.FileChannel=chan;
        h.phyChan=Info.phyChan;
        h.kind=Info.kind;
        h.npoints=NumberOfMarkers;
        h.comment=Info.comment;
        h.title=Info.title;
    end;
    
    [data.timings,h.TimeUnits]=convertSamplesToTime_smr(fID,data.timings, varargin{:});                % Convert time
    h.Epochs={startBlock endBlock 'of' Info.blocks 'blocks'};
    if ShowProgress==1
        close(progbar);
        drawnow;
    end;
end
