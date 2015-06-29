function [h] = barsig2(D,options,raw)

% PLOT BAR GRAPH OF MEANS +/- SEM FOR 2 OR MORE VARIABLES FOR 2 GROUPS:
% plot significance values from a t-test or Wilcoxon rank-sum test on bar graphs
% x: x-location of bars on the chart
% D: data as cell array of variables, trial values within each cell
% a: critical alpha level for rank-sum test
% options.stars: use lines and options.stars only if signif diff
% options.siglvl: use 1 to 3 options.stars for signif. at a, a/10, a/100 levels
% options.color: RGB color codes, one color or a matrix of [RGB;RGB;...;RGB]
% options.checktesttype: use 1 to check that selected test is appropriate
% (equal sample sizes for ttest or normally distributed data for ttest2),
% use 0 to ignore
% options.group: only plot comparisons between groups (Dim 2)
% options.mean: display means
% labels: x-axis labels for the bars, e.g. {'Var 1' 'Var 2' 'Var 3'}
% test_type: 'ttest2', 'ttest', 'ranksum'
% errtype: 1 SEM; 2 std. dev.; 3 95% conf. int.
% trls: Plot text of # of samples and % of samples that are real values
% samps: # of samples that percentage is taken from, one for each percentage value in same size as D
% Nathan Killian 130115

if nargin<3,raw=[];end
if nargin<2,options = [];end

nvars = size(D,1);
nGroups = size(D,2);

testprop=0;
if ~iscell(D)
    testprop=1;
%     D =  mat2cell(D,size(D,1),ones(1,size(D,2)));
    D =  mat2cell(D,ones(1,nvars),ones(1,nGroups));
    % Total count
    DT = mat2cell(raw{2},ones(1,nvars),ones(1,nGroups));
    % Raw count
    DC = mat2cell(raw{1},ones(1,nvars),ones(1,nGroups));
end

options = setdefaults(options,'diff',0,'mean',0,'trls',0,'group',0,'checktesttype',1,'labels',[],'errtype',1,'test_type','ranksum','color',[0.5 0.5 0.5],'plotsig',1,'siglvl',1,'stars',0,'a',0.05,'x',[1:nvars],'cmap','nicejet','dotplot',1,'dotmsize',2,'dotcolor','k');

if options.diff
    rEnd=4;
else
    rEnd=3;
end

for d=1:2
    for k = 1:nvars
        nsamp(k,d) = length(D{k,d});
    end
end

if options.checktesttype
    if length(unique(nsamp))>1
        if strcmp(options.test_type,'ttest')
            fprintf('Different sample sizes, using ttest2 instead of ttest\n')
            options.test_type = 'ttest2';
        end
    end
    
    % TEST FOR NORMALITY
    if strcmp(options.test_type,'ttest2') || strcmp(options.test_type,'ttest')
        fprintf('ttest selected, testing for normality\n')
        useRS = zeros(nvars,1);
        for d=1:2
            for k = 1:nvars
                [h,p,ksstat,cv] = kstest(D{k,d});
                if h,
                    fprintf('The data in variable %g are not normally distributed, using nonparametric rank-sum test on medians\n',k);
                    useRS(k,d) = 1;
                end
            end
        end
        if any(useRS), options.test_type = 'ranksum';else fprintf('normality was not violated.\n');end
    end
    if testprop
        options.test_type = 'chisquare';
    end
end

% do the bar graph
for d=1:2;
    for k = 1:nvars
        mu(k,d) = nanmean(D{k,d});
        sg(k,d) = nanstd(D{k,d});
        trls(k,d)=sum(~isnan(D{k,d}));
        valz(k,d)=sum(D{k,d}~=0);
%         prop(k,d)=round((trls(k,d)/size(D{k,d},1))*100);
        if strcmp(options.test_type,'ttest')
            prop(k,d)=round((valz(k,d)/size(D{k,d},1))*100);
        else
            prop(k,d)=round((trls(k,d)/size(D{k,d},1))*100);
        end
        sems(k,d) = sg(k,d)./sqrt(trls(k,d));% SEM
    end
end
if options.errtype == 2
    sems = sg;%use the standard deviation instead
elseif options.errtype == 3
    sems = sems*1.96;
end

figure;h=bar(mu,'hist');

if nvars>1
    x=cell2mat(get(h,'Xdata'));
else
    x=get(h,'Xdata');
end
% x=cell2mat(get(h,'Xdata'));
xcenter = 0.5*(x(2:4:end,:)+x(3:4:end,:));
options.x=reshape(xcenter,1,numel(xcenter));
close(gcf)

% For some reason this isn't working 10.30.13
figure;barwitherr(sems,mu);
colormap(options.color);
yR=get(gca,'YLim');
% 
% [maxVal maxInd]=max(mu);
% [maxVal2 maxInd2]=max(maxVal);
% yR(2)=(maxVal2+sems(maxInd(maxInd2),maxInd2))*1.1;

scaleFactor=[.04,.1,.16,.22];
startVal=.1;

% Add text below Xticks, text starts at xTick
% options.xG=[((options.x(1)+options.x(2))/2),((options.x(3)+options.x(4))/2)]
% text(options.xG,[-8,-8],'hello')

if options.trls
    if options.mean
        c=1;
        for k=1:nvars
            for d=1:nGroups
                r=1;
                plotStr(r,c)={[num2str(prop(k,d)),'%']};
                textY(r,c)=yR(2)*scaleFactor(r);            
                
                r=2;
                plotStr(r,c)={num2str(trls(k,d))};
                textY(r,c)=yR(2)*scaleFactor(r); 
                
                
                r=3;
                plotStr(r,c)={num2str(round(mu(k,d)))};
                textY(r,c)=yR(2)*scaleFactor(r);
                if d==2
                    r=4;
                    gDiff=round(mu(k,1))-round(mu(k,2));
                    l=1;
                    m=1;
                    if gDiff<0
                        m=2;
                    else gDiff>0
                        l=2;
                    end
                    plotStr(r,c)={num2str(abs(gDiff))};
                    textY(r,c)=yR(2)*scaleFactor(r);
                end
                c=c+1;
            end
        end
        for r=1:rEnd
            t1=text(options.x-startVal,textY(r,:),plotStr(r,:));
%             if r==4
%                 if gDiff<0
%                     set(t1,'Color','r')
%                 elseif gDiff>0
%                     set(t1,'Color','g')
%                 else
%                     set(t1,'Color','w')
%                 end
%             else
                set(t1,'Color','g')
%             end
        end
    else
        c=1;
        for k=1:nvars
            for d=1:nGroups
                r=1;
                plotStr(r,c)={num2str(trls(k,d))};
                textY(r,c)=yR(2)*scaleFactor(r);
                
                r=2;
                plotStr(r,c)={[num2str(prop(k,d)),'%']};
                textY(r,c)=yR(2)*scaleFactor(r);
               
                c=c+1;
            end
        end
        for r=1:3
            t1=text(options.x-startVal,textY(r,:),plotStr(r,:));
            set(t1,'Color','g')            
        end
    end
elseif options.mean
    c=1;
    for k=1:nvars
        for d=1:nGroups            
            r=1;
            plotStr(r,c)={num2str(round(mu(k,d)))};
            textY(r,c)=yR(2)*scaleFactor(r);
            
            c=c+1;
        end
    end 
    for r=1:3
        t1=text(options.x-startVal,textY(r,:),plotStr(r,:));
        set(t1,'Color','g')
    end
end



% 
% 
%    
% if options.trls && ~options.mean        
%     c=1;
%     for k=1:nvars
%         for d=1:nGroups
%             trlStr(c)={num2str(trls(k,d))};
%             textY(1,c)=yR(2)*.1;
%             propStr(c)={[num2str(prop(k,d)),'%']};
%             textY(2,c)=yR(2)*.04;
%             c=c+1;
%         end
%     end   
% elseif ~options.trls && options.mean
%     c=1;
%     for k=1:nvars
%         for d=1:nGroups
%             plotStr(1,c)={num2str(trls(k,d))};
%             textY(1,c)=yR(2)*.1;
%             propStr(2,c)={[num2str(prop(k,d)),'%']};
%             textY(2,c)=yR(2)*.04;
%             c=c+1;
%         end
%     end
% elseif options.trls && options.mean
%     c=1;
%     for k=1:nvars
%         for d=1:nGroups
%             trlStr(1,c)={num2str(trls(k,d))};
%             textY(1,c)=yR(2)*.1;
%             propStr(2,c)={[num2str(prop(k,d)),'%']};
%             textY(2,c)=yR(2)*.04;
%             c=c+1;
%         end
%     end
% end
% 
% if c>1
%     startVal=.09;
%     t1=text(options.x-startVal,textY(1,:),trlStr);
%     t2=text(options.x-startVal,textY(2,:),propStr);
%     set(t1,'Color','g')
%     set(t2,'Color','g')
% end

if options.dotplot
    for k = 1:length(options.x)
        size(D{k})
        ph = plotSpread(makecol(D{k}),'distributionColors',options.dotcolor,'xvalues',options.x(k));
        set(ph{1},'markersize',options.dotmsize)
    end
end

set(gca,'xtick',1:nvars)
if ~isempty(options.labels),set(gca,'xticklabel',options.labels);end


% plot some significance info

if options.stars
    if options.group
%         disp('plotting with ''stars''')
        sigi=1;
        for k=1:nvars
            if (trls(k,1)~=0 && trls(k,2)~=0 && valz(k,1)>1 && valz(k,2)>1) || testprop
                switch options.test_type
                    case 'ttest'
                        [h p] = ttest(D{k,1},D{k,2});
                    case 'ttest2'
                        [h p] = ttest2(D{k,1},D{k,2});
                    case 'ranksum'
                        [p h] = ranksum(D{k,1},D{k,2});
                    case 'chisquare'
                        [h,p,chi2stat,df]=prop_test([DC{k,1} DC{k,2}],[DT{k,1} DT{k,2}],0);
%                         [h,p,chi2stat,df]=prop_test([((D{k,1}/100)*DT{k,1}) ((D{k,2}/100)*D{k,2})],[DT{k,1} D{k,2}],0);
                end
                top = max([nanstd(D{k,1},1)/sqrt(size(D{k,1},1)*.5)+nanmean(D{k,1},1) nanstd(D{k,2},1)/sqrt(size(D{k,2},1)*.5)+nanmean(D{k,2},1)]);
                xr = [options.x((k*2)-1) options.x(k*2)];
                xmid = diff(xr)/2+xr(1);
                if p<options.a & options.plotsig
                    if options.siglvl, if p<options.a/100, lvl = 3;elseif p<options.a/10, lvl = 2;else lvl = 1;end;else lvl = 1;end
                    yshift = sigi*top/20;top2 = top+yshift;
                    line(xr,[top2 top2],'color',[0 0 0],'linewidth',2)
                    text(xmid,top2+top2/60,repmat(['*'],1,lvl),'fontsize',14,'fontweight','bold','color',[0 0 0],'horizontalalignment','center')
                    sigi = sigi+1;
                end
            end
        end
    elseif options.grouplvl
        disp('plotting with ''stars''')
        sigi=1;
        for k=1:nvars
            switch options.test_type
                case 'ttest'
                    [h p] = ttest(D{k,1},D{k,2});
                case 'ttest2'
                    [h p] = ttest2(D{k,1},D{k,2});
                case 'ranksum'
                    [p h] = ranksum(D{k,1},D{k,2});
            end
            top = max([nanstd(D{k,1},1)/sqrt(size(D{k,1},1)*.5)+nanmean(D{k,1},1) nanstd(D{k,2},1)/sqrt(size(D{k,2},1)*.5)+nanmean(D{k,2},1)]);
            xr = [options.x((k*2)-1) options.x(k*2)];
            xmid = diff(xr)/2+xr(1);
            if p<options.a & options.plotsig
                if options.siglvl, if p<options.a/100, lvl = 3;elseif p<options.a/10, lvl = 2;else lvl = 1;end;else lvl = 1;end
                yshift = sigi*top/20;top2 = top+yshift;
                line(xr,[top2 top2],'color',[0 0 0],'linewidth',2)
                text(xmid,top2+top2/60,repmat(['*'],1,lvl),'fontsize',14,'fontweight','bold','color',[0 0 0],'horizontalalignment','center')
                sigi = sigi+1;
            end
        end        
    else
        disp('plotting with ''stars''')
        pi = 1;maini = 1;
        % Convert matrix into vector
        muV=reshape(mu',1,numel(mu));
        semsV=reshape(sems',1,numel(sems));
        [dum heighti] = sort(muV,'descend');
        kset = heighti(1:end-1);%tallest first, last one already eval'd
        kkset = setdiff([1:length(options.x)],kset(1));% all others
        [dum newi] = sort(abs([kset(1)-kkset]),'descend');kkset = kkset(newi);% need to keep longest line on bottom
        for k = kset
            sigi = 0;
            for kk = kkset(maini:end);
                switch options.test_type
                    case 'ttest'
                        [h p] = ttest(D{k},D{kk});
                    case 'ttest2'
                        [h p] = ttest2(D{k},D{kk});
                    case 'ranksum'
                        [p h] = ranksum(D{k},D{kk});
                end
                top = max([nanstd(D{k},1)/sqrt(size(D{k},1)*.5)+nanmean(D{k},1) nanstd(D{kk},1)/sqrt(size(D{kk},1)*.5)+nanmean(D{kk},1)]);
                xr = [options.x(k) options.x(kk)];
                xmid = diff(xr)/2+xr(1);
                if p<options.a & options.plotsig
                    if options.siglvl, if p<options.a/100, lvl = 3;elseif p<options.a/10, lvl = 2;else lvl = 1;end;else lvl = 1;end
                    yshift = sigi*top/20;top2 = top+yshift;
                    line(xr,[top2 top2],'color',[0 0 0],'linewidth',2)
                    text(xmid,top2+top2/60,repmat(['*'],1,lvl),'fontsize',14,'fontweight','bold','color',[0 0 0],'horizontalalignment','center')
                    sigi = sigi+1;
                end
                pi = pi+1;
            end
            maini = maini + 1;
        end
    end
else
    if nchoosek(length(options.x),2)>1
        %         clrs = gencs(nchoosek(length(x),2));
        clrs = getmap(nchoosek(length(options.x),2),options.cmap);
    else
        clrs = [0 0 1];
    end
    pi = 1;
    for k = 1:length(options.x)-1
        for kk = k+1:length(options.x)
            switch options.test_type
                case 'ttest'
                    [h p] = ttest(D{k},D{kk});
                case 'ttest2'
                    [h p] = ttest2(D{k},D{kk});
                case 'ranksum'
                    [p h] = ranksum(D{k},D{kk});
            end
            top = nanmean([nanmean(D{k},1) nanmean(D{kk},1)]);
            xr = [options.x(k) options.x(kk)];
            xmid = diff(xr)/3+xr(1);
            if p<options.a & options.plotsig
                line(xr,[top top],'color',clrs(pi,:),'linewidth',3)
                text(xmid,top+top/10,['p=' num2str(rsig(p,4))],'fontsize',14,'fontweight','bold','color',clrs(pi,:),'horizontalalignment','center')
            elseif options.plotsig
                line(xr,[top top],'color',clrs(pi,:),'linewidth',1)
                text(xmid,top+top/10,['p=' num2str(rsig(p,4))],'fontsize',12,'color',clrs(pi,:),'horizontalalignment','center')
            end
            pi = pi+1;
        end
    end
end


%% sub functions
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
    end
end
return;

function [cset] = getmap(n,map)
%get a colorset based on a set number of variables
% and a specified map
% Nathan Killian 110827
if nargin < 2
    map = 'jet';
end
colormap('default');
map = flipud(colormap(map));
cset = map(round(linspace(1,size(map,1),n)),:);
cset = cset/max(cset(:));
% cset = map(round(linspace(1,size(map,1),n+2)),:);
% cset(1,:) = [];cset(2,:) = [];
return;

function map = nicejet(numvals)
if nargin<1, numvals = 0;end
jet = [0,0,0.500000000000000;0,0,0.563492063492064;0,0,0.626984126984127;0,0,0.690476190476191;0,0,0.753968253968254;0,0,0.817460317460317;0,0,0.880952380952381;0,0,0.944444444444444;0,0.00793650793650791,1;0,0.0714285714285714,1;0,0.134920634920635,1;0,0.198412698412698,1;0,0.261904761904762,1;0,0.325396825396825,1;0,0.388888888888889,1;0,0.452380952380952,1;0,0.515873015873016,1;0,0.579365079365079,1;0,0.642857142857143,1;0,0.706349206349206,1;0,0.769841269841270,1;0,0.833333333333333,1;0,0.896825396825397,1;0,0.960317460317460,1;0.0238095238095237,1,0.976190476190476;0.0873015873015872,1,0.912698412698413;0.150793650793651,1,0.849206349206349;0.214285714285714,1,0.785714285714286;0.277777777777778,1,0.722222222222222;0.341269841269841,1,0.658730158730159;0.404761904761905,1,0.595238095238095;0.468253968253968,1,0.531746031746032;0.531746031746032,1,0.468253968253968;0.595238095238095,1,0.404761904761905;0.658730158730159,1,0.341269841269841;0.722222222222222,1,0.277777777777778;0.785714285714286,1,0.214285714285714;0.849206349206349,1,0.150793650793651;0.912698412698413,1,0.0873015873015874;0.976190476190476,1,0.0238095238095237;1,0.960317460317461,0;1,0.896825396825397,0;1,0.833333333333334,0;1,0.769841269841270,0;1,0.706349206349207,0;1,0.642857142857143,0;1,0.579365079365080,0;1,0.515873015873016,0;1,0.452380952380953,0;1,0.388888888888889,0;1,0.325396825396826,0;1,0.261904761904762,0;1,0.198412698412699,0;1,0.134920634920635,0;1,0.0714285714285716,0;1,0.00793650793650791,0;0.944444444444445,0,0;0.880952380952381,0,0;0.817460317460318,0,0;0.753968253968254,0,0;0.690476190476191,0,0;0.626984126984127,0,0;0.563492063492064,0,0;0.500000000000000,0,0;];
%
map = jet(9:end-8,:);
if numvals
    samps = round(linspace(1,size(map,1),numvals));
    map =  map(samps,:);
end
return;

function [colorset] = gencs(N)
%generate a options.colorset with N options.colors, use linear spacing to avoid sameness
%and double randomization of RGB to avoid sequential sameness
%this adds a bit more fun to things than just linearly interpolating a options.colormap
%Nathan Killian 2/21/09

Rset = linspace(0,1,N);
Gset = linspace(0,1,N);
Bset = linspace(0,1,N);

colorset0 = [randsample(Rset,N,false)',randsample(Gset,N,false)',...
    randsample(Bset,N,false)'];
rearrange1 = randperm(N);rearrange2 = randperm(N);rearrange3 = randperm(N);

colorset = [colorset0(rearrange1,1), colorset0(rearrange2,2),...
    colorset0(rearrange3,3)];
return;