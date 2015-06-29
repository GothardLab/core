function d = setdefaults(varargin)
%d = var2field([],variable name, variable value, rinse and repeat...)
% Set Defaults
%This is a very useful way to put existing variables into structure fields
%first argument must be the structure variable (empty or otherwise)

%Nathan Killian 110215

lastvar = nargin-1;
d = varargin{1};
for k = 2:2:lastvar
    try
        if ~isfield(d,varargin{k})
            d = setfield(d,varargin{k},varargin{k+1});
        end
    catch
        %         %         names = fieldnames(varargin{k});
        %         inputname(k);
        %         period = strfind(inputname(k),'.');
        %         %         d = setfield(d,inputname(k)(period+1:end),varargin{k});
    end
end