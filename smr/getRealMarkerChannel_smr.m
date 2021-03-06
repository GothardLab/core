function[data,h]=getRealMarkerChannel_smr(fID, chan, varargin)
% getRealMarkerChannel_smr reads an RealMark channel from a file.
%
% [DATA {, HEADER}]=getRealMarkerChannel_smr(fID, CHAN{, START{, STOP{, OPTIONS}}})
% fID is the matlab file handle, CHAN is the channel number (1=max)
%
% [DATA, HEADER]=getRealMarkerChannel_smr(fID, 1{, OPTIONS})
%       reads all the data on channel 1
% [DATA, HEADER]=getRealMarkerChannel_smr(fID, 1, 10{, OPTIONS})
%       reads disc block 10 for continuous data or epoch 10 for triggered
%       data
% [DATA, HEADER]=getRealMarkerChannel_smr(fID, 1, 10, 20{, OPTIONS})
%       reads disc blocks 10-20
%
% DATA is a structure with 3 fields.
%       DATA.TIMINGS contains timestamps
%       DATA.MARKERS contains 4 uint8 marker values for each event
%       DATA.ADC contains the realwave data associated with each timestamp
%
% When present, OPTIONS must be the last input argument. Valid options
% are:
% 'ticks', 'microseconds', 'milliseconds' and 'seconds' cause times to
%    be scaled to the appropriate unit (seconds by default)in HEADER
% 'scale' - no effects
% 'progress' - causes a progress bar to be displayed during the read.

Info=getInfo_smr(fID,chan);

if isempty (Info)
    data=[];
    h=[];
    return;
end;

if(Info.kind~=7)
    warning('getRealMarkerChannel_smr: Channel %d No data or wrong channel type',chan);
    return;
end;


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


FileH=getHeader_smr(fID);
SizeOfHeader=20;                                            % Block header is 20 bytes long
header=getBlockHeaders_smr(fID,chan);
if isempty(header)
     warning('getRealMarkerChannel_smr: Channel number #%d has blank header',chan);
     data=[];
    h=[];
     return;
else
    NumberOfMarkers=sum(header(5,startBlock:endBlock)); % Sum of samples in required blocks
    
    nValues=Info.nExtra/4;                                      % Each value has 4 bytes (single precision)
    data.timings=zeros(NumberOfMarkers,1);
    data.markers=char(zeros(NumberOfMarkers,4));
    data.real=single(zeros(NumberOfMarkers,nValues));
    
    count=1;
    for block=startBlock:endBlock
        fseek(fID, header(1, block)+SizeOfHeader, 'bof');                  % Start of block
        for i=1:header(5,block)                                            % loop for each marker
            data.timings(count)=fread(fID,1,'int32');            % Time
            data.markers(count,:)=fread(fID,4,'uint8=>uint8');             % 4x marker bytes
            data.real(count,:)=fread(fID,nValues,'single=>single');
            count=count+1;
        end;
        if ShowProgress==1
            done=(block-startBlock)/max(1,endBlock-startBlock);
            progressbar(done, progbar,sprintf('Reading Channel %d....     %d%% Done',chan,(int16(done*100)/5)*5));
        end
    end
    
    
    if(nargout>1)
        h.FileName=Info.FileName;                                   % Set up the header information to return
        h.system=['Matoff' num2str(FileH.systemID)];
        h.FileChannel=chan;
        h.phyChan=Info.phyChan;
        h.kind=Info.kind;
        h.npoints=NumberOfMarkers;
        h.values=Info.nExtra/4;
        h.preTrig=Info.preTrig;
        h.comment=Info.comment;
        h.title=Info.title;
        h.sampleinterval=getSampleInterval_smr(fID,chan);
        h.min=Info.min;
        h.max=Info.max;
        h.units=Info.units;
        if(isfield(Info,'interleave'))
            h.interleave=Info.interleave;
        end
    end
    
    [data.timings,h.TimeUnits]=convertSamplesToTime_smr(fID,data.timings,varargin{:});                % Convert time
    h.Epochs={startBlock endBlock 'of' Info.blocks 'blocks'};
    if ShowProgress==1
        close(progbar);
        drawnow;
    end
end