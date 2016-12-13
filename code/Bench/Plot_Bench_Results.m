function Plot_Bench_Results(benchRes,benchName,figures_path)
if ~exist('benchName','var')
    benchName = 'default';
end
if ~exist('figures_path','var')
    code_path = fileparts(mfilename('fullpath'));
    root_path = fileparts(code_path);
    data_path = fullfile(root_path,'data');
    figures_path = fullfile(data_path,'figures',benchName);
	MakeDir(figures_path)
    clear code_path root_path data_path
end

%% params

transformationClasses = GetShrec2010ClassNames();
nClasses = numel(transformationClasses);
nStrengh = 5;
nOverlapVals = 30;
overlapVals = linspace(0,1,nOverlapVals);
%plotFormats = {'ks-','bo-','gd-','r+-','m*-','kx-','ch-','g-.','r-.','k-.'};

strenghValsToPlot = [3 5];
%%
nPassTh   = @(x,th) nnz(x > th);
pctPassTh = @(x,th) nPassTh(x,th) * 100 / numel(x);

% plots to do:
template_ = struct(...
    'repeatability',    zeros(nOverlapVals,1),...
    'nCorrespondences', zeros(nOverlapVals,1),...
    'matchingScore',    zeros(nOverlapVals,1),...
    'nCorrectMatches',    zeros(nOverlapVals,1),...
    'matchingScore_SI',    zeros(nOverlapVals,1),...
    'nCorrectMatches_SI',  zeros(nOverlapVals,1));

strenghAvg = repmat(template_,nStrengh,1);
for maxStrengh = 1:nStrengh
    
    classesVals = repmat(template_,nClasses,1);
    
    for iClass = 1:nClasses
        currClass = transformationClasses{iClass};
        classRes  = [benchRes.(currClass)];
%         classRes  = reshape(classRes,nStrengh,3);
        classResXform  = [classRes(1:maxStrengh,:).xform];
        classResMats  = [classRes(1:maxStrengh,:).mats];

        classesVals(iClass).nCorrespondences = ComputeClassVals(nPassTh,  {classResXform.maxIntersection},overlapVals);
        classesVals(iClass).repeatability    = ComputeClassVals(pctPassTh,{classResXform.maxIntersection},overlapVals);
        
        matchIntersection = [classResXform.matchIntersection];
        classesVals(iClass).nCorrectMatches     = ComputeClassVals(nPassTh,  {matchIntersection.vs_12},overlapVals);
        classesVals(iClass).matchingScore       = ComputeClassVals(pctPassTh,{matchIntersection.vs_12},overlapVals);
        classesVals(iClass).nCorrectMatches_SI  = ComputeClassVals(nPassTh,  {matchIntersection.siDesc},overlapVals);
        classesVals(iClass).matchingScore_SI    = ComputeClassVals(pctPassTh,{matchIntersection.siDesc},overlapVals);
    end
    
    fN = fieldnames(classesVals);
    for iF = 1:numel(fN)
        strenghAvg(maxStrengh).(fN{iF}) = mean([classesVals.(fN{iF})],2);
    end
    
    classesVals(nClasses + 1) = strenghAvg(maxStrengh);
    if ismember(maxStrengh,strenghValsToPlot)
        FormatAndSave(...
            overlapVals,...
            classesVals,...
            fullfile(figures_path ,['PerClass_le' int2str(maxStrengh)]),...
            [transformationClasses,{'AVERAGE'}],...
            benchName)
    end
    
end

avgLegend = num2cell([repmat('strengh \leq ',nStrengh,1) int2str((1:nStrengh)')],2);
FormatAndSave(...
    overlapVals,...
    strenghAvg,...
    fullfile(figures_path ,'PerStrengh' ),...
    avgLegend,...
    benchName)



function FormatAndSave(xVals,dataStruct,figure_path,legendVal,benchName)
markerFormats = {'s','o','d','+','*','x','p','<','>'};
plotTypes = fieldnames(dataStruct);
for iPlot = 1:numel(plotTypes)
    fHandle = figure;
    lineHs = plot(xVals,[dataStruct.(plotTypes{iPlot})]');
    for iLine = 1:numel(lineHs)
        if strcmp(legendVal{iLine},'AVERAGE')
            set(lineHs(iLine),'LineWidth',3,'color','k','LineStyle','-.');
        else set(lineHs(iLine),'Marker',markerFormats{iLine}); 
        end
    end
    grid on
    xlabel('overlap')
    SetYaxisProps(plotTypes{iPlot})
    %if iPlot==1,legend(legendVal,'Location','Best');end
    
    currFigPath = [figure_path '_' plotTypes{iPlot}];
    saveas(fHandle,[currFigPath,'.fig'])
    title(benchName)
    saveas(fHandle,[currFigPath,'.jpg'])
    close(fHandle)
end


function resVals = ComputeClassVals(fHandle,valsCell,thVals)
nVals = numel(thVals);
resVals = zeros(nVals,1);
for iVal = 1:nVals
    currThFun = @(x) fHandle(x,thVals(iVal));
    
    transformedVals = cellfun(currThFun,valsCell);
    transformedMean = mean(transformedVals);
    
    resVals(iVal) = transformedMean;    
end

function SetYaxisProps(plotType)
plotType = strrep(plotType,'_SI','');
switch plotType
    case 'nCorrespondences'
        ylabel('# of correspondences')
    case 'repeatability'
        ylabel('repeatability (%)')
    case 'nCorrectMatches'
        ylabel('# of correct matches')
    case 'matchingScore'
        ylabel('matching score (%)')
end
if ~strncmp(plotType,'n',1)
    ylim([0 100])
end
