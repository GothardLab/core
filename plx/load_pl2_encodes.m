function [codes, ts] = load_pl2_encodes(path, method)

% Load the pl2 file index
pl2 = PL2GetFileIndex(path);


% Use default if no inputs specified
switch nargin
    case 2
        method = upper(method);
    case 1
        method = 'DEFAULT';
    otherwise
        error('Not enough inputs');
end


%Read in the first strobed words
[strobed] = PL2EventTs(path, 'Strobed');

% Seperate times and strobe values
vals = strobed.Strobed; %Sometimes off by +256 or -32512
times = strobed.Ts;

if strcmp(method, 'DEFAULT')
    
    codes = vals;
    ts = times;
    
elseif strcmp(method, 'HILO')
    
    tooClose = 0.002; %Closest two paired encodes can be, usually 0.002 seconds
    tooFar =0.005; %Farthest two paired encodes can be, usually 0.006 seconds
    
    timeDiffs = diff(times);
    
    c = 0;
    
    for i = 1:(size(vals,1)-1)
        if (timeDiffs(i)>tooClose && timeDiffs(i)<tooFar)
            c=c+1;
            lob(c)=vals(i);
            hib(c)=vals(i+1);
            ts(c)=times(i);
            
        end
    end
    
    codes = nan(1,c);
    
    for i=1:length(lob);
        codes(i)=double(lob(i))+(double(hib(i))*256);
    end
    
end



