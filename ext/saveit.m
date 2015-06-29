function saveit(titleLabel,newDir,figOptions)
%SAVEIT Saves current figure as a PNG and/or EPS
% saveit(titleLabel,newDir,figOptions)
% figOptions.saveEPS or .savePNG, save PNG only by default
% figOptions.close, close by default
% Requires setdefaults,export_fig and the software to use export_fig 
% JA Solyst made dis 140316
figOptions = setdefaults(figOptions,'close',1,'saveEPS',0,'savePNG',1,'save',1);

if figOptions.save
    if ~exist(newDir,'dir');
        mkdir(newDir);
    end
    
    
    saveName=[newDir,titleLabel];
    
    if figOptions.saveEPS
        if figOptions.savePNG
            export_fig([saveName,'.eps'],'-eps','-transparent')
            title(titleLabel)
            export_fig([saveName,'.png'],'-png')
        end
    elseif figOptions.savePNG
        title(titleLabel)
        export_fig([saveName,'.png'],'-png')
    else
    end
end
    if figOptions.close
        close(gcf)
    end

end

