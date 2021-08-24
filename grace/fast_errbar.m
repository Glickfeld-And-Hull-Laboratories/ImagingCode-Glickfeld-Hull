function fast_errbar(xvals,data,trialDim,varargin)

p = inputParser;
p.addParamValue('shaded', false, @islogical);
p.addParamValue('continuous',true, @islogical);
p.addParamValue('color', [0 0 0], @isnumeric);
p.addParamValue('marker_size',6, @isnumeric);
p.addParamValue('norm',false,@islogical);
p.addParamValue('SD',false,@islogical);
p.addParamValue('stats',false,@islogical);
p.addParamValue('cells_as_x',true,@islogical);

% parse inputs
p.parse(varargin{:});
params = p.Results;

if ~isempty(data)
if params.shaded
    if isempty(xvals)
            xvals = 1:numel(nanmean(data,trialDim));
    end

    mean_val = smooth(nanmean(data,trialDim));
    if params.SD
        std_val = nanstd(data,[],trialDim);
    else
        std_val = nanstd(data,[],trialDim)/sqrt(size(data,trialDim));
    end

    h = shadedErrorBar(xvals,mean_val,std_val);
    h.mainLine.Color = params.color;
    h.patch.FaceColor = params.color;
else
    if iscell(data)
        hold on;
        for cell_i = 1:numel(data)
            if params.SD
                    var = nanstd(data{cell_i},[],trialDim);
            else
                var = nanstd(data{cell_i},[],trialDim)/sqrt(numel(data{cell_i}));
            end
            if isempty(xvals)
                xvals = 1:numel(nanmean(data{cell_i},trialDim));
            end
            if params.cells_as_x
                    x_idx = xvals(cell_i);
            else
                x_idx = xvals;
            end
            if params.continuous
                if size(params.color,1) > 1
                    errorbar(x_idx,nanmean(data{cell_i},trialDim),var,'o-','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color(cell_i,:),'Color',params.color(cell_i,:))
                else 
                    errorbar(x_idx,nanmean(data{cell_i},trialDim),var,'o-','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color,'Color',params.color);
                end
            else
                if size(params.color,1) > 1
                    errorbar(x_idx,nanmean(data{cell_i},trialDim),var,'o','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color(cell_i,:),'Color',params.color(cell_i,:));
                else 
                    errorbar(x_idx,nanmean(data{cell_i},trialDim),var,'o','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color,'Color',params.color);
                end
            end
                
        end
        if params.stats
            nDist = numel(data);
            distIdx = num2cell(1:nDist,1);
            nCellsByDist = cellfun(@numel,data,'un',0);
            
            temp = cellfun(@(x,y) repmat(x,y,1),distIdx,nCellsByDist,'un',0);
            dist_ids = vertcat(temp{:});
            try
                [p,~,stats_result] = anova1(cell2mat(data'),dist_ids,'off');
            catch
                [p,~,stats_result] = anova1(cell2mat(data),dist_ids,'off');
            end
            disp(p);
            comparison = multcompare(stats_result,'display','off');
            disp(comparison);
        end
    else
        if isempty(xvals)
            xvals = 1:size(nanmean(data,trialDim));
        end
        if params.norm
%             tempy = nanmean(data,trialDim)./max(nanmean(data,trialDim));
%             tempstd = nanstd(data,[],trialDim)/(max(nanmean(data,trialDim))*sqrt(size(data,trialDim)));
            tempy = nanmean(data,trialDim);
            if params.SD
                tempstd = nanstd(data,[],trialDim)/(tempy(1));
            else
                tempstd = nanstd(data,[],trialDim)/(tempy(1)*sqrt(size(data,trialDim)));
            end
            tempy = tempy/tempy(1);
        else
            tempy = nanmean(data,trialDim);
            if params.SD
                tempstd = nanstd(data,[],trialDim); 
            else
                tempstd = nanstd(data,[],trialDim)/(sqrt(size(data,trialDim))); 
            end
        end
        if params.continuous
            errorbar(xvals,tempy,tempstd,'o-','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color,'Color',params.color);
        else
            errorbar(xvals,tempy,tempstd,'o','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color,'Color',params.color);
        end
        if params.stats
%             disp('no stats support right now for matrix');
              if trialDim == 2
                  data = data';
              end
            [p,~,stats_result] = anova1(data,[],'off');
            disp(p);
            comparison = multcompare(stats_result,'display','off');
            disp(comparison);
              
        end
    end
end
    
end