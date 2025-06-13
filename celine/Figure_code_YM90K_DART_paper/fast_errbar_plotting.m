function fast_errbar_plotting(xvals,data,trialDim,varargin)

p = inputParser;
p.addParamValue('color', [0 0 0], @isnumeric);
p.addParamValue('marker_size',6, @isnumeric);
p.addParamValue('cap_size',10,@isnumeric);

% parse inputs
p.parse(varargin{:});
params = p.Results;

    
if isempty(xvals)
    xvals = 1:size(nanmean(data,trialDim));
end
        
tempy = mean(data,trialDim);
tempstd = std(data,[],trialDim)/(sqrt(size(data,trialDim))); 

h = errorbar(xvals,tempy,tempstd,'o-','LineWidth',1.5,'MarkerSize',params.marker_size,'MarkerFaceColor',params.color,'Color',params.color);

        
h.CapSize = params.cap_size;