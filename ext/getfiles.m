function files = getfiles(datadir,datafiletype,o)
% function files = getfiles(datadir,datafiletype,o)
% get a list of files, 'advanced'
% options fields are leavein, pullout, keep, remove, strind, useOnlyValid, fullpath
%   EXAMPLE USAGE:
%     [Num Txt xlsRaw] = xlsread(['D:\Dropbox\DB\ArrayInfo.xls']);
%     o.valid = xlsRaw(find(strcmp(structure,xlsRaw(2:end,8))&[xlsRaw{2:end,11}]')+1,1);
%     o.useOnlyValid = 1;o.strind = 1:8;
%     o.leavein = 0;    o.keep    = {'100614'};
%     o.pullout = 1;    o.remove  = {'070502' '060803'};
%     nexfiles = getfiles(nexdir,nexfiletype,o)
%     
% Nathan Killian 100101


if ~(strcmp(datadir(end),'/')||strcmp(datadir(end),'\'))
    datadir = [datadir '/'];
end
if nargin<3, o = [];end
o = setdefaults(o,'leavein',0,'keep',{''},'pullout',0,'remove',{''},'strind',[3:8],'fullpath',0,'useOnlyValid',0);

datstr = strcat(datadir,datafiletype);
disp(['finding ' datstr])
d = dir(datstr);

files = [];
count = 1;
for k=1:length(d)
    if o.pullout
        if any(strcmp(d(k).name(o.strind),o.remove))
            continue
        end
    end
    if o.leavein
        if any(strcmp(d(k).name(o.strind),o.keep))
            files{count,1}=d(k).name;
            count = count + 1;
        end
    elseif o.useOnlyValid
        if any(strcmp(d(k).name(o.strind),o.valid))
            files{count,1}=d(k).name;
            count = count + 1;
        end
    else
        files{count,1}=d(k).name;
        count = count + 1;
    end
end

if o.fullpath
    for k = 1:length(files)
        files{k,1} = [datadir files{k,1}];
    end
end