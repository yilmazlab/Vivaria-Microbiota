% Function: Print annotation statistics
function varargout = fAnnotationStats(tableRxns, varargin)

if nargin == 1
    nrxn = height(tableRxns);
else
    nrxn = varargin{1};
end

stats = [sum(tableRxns.nEmpty>0), ...
    sum(tableRxns.nMultiple>0), ...
    sum(~(tableRxns.rxnMetaNetXID=="")), ...
    sum(~(tableRxns.rxnSEEDID=="")), ...
    sum(~(tableRxns.rxnKEGGID==""))];
rstats = 100 * stats / nrxn;

fprintf('[%s] %5i rxns: %3.1f%% empty, %3.1f%% multiple, %3.1f%% MetaNetX, %3.1f%% SEED, %3.1f%% KEGG.\n', ...
    datestr(now), nrxn, rstats);

varargout{1} = stats;

return