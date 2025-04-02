%% Function: Consolidate annotations for AGORA2 models

function [tableRxns,varargout] = fConsolidateModels(tableRxns,tableSubsExceptions, ...
    tableRxnExceptions,options,varargin)

ner = numel(unique(tableRxnExceptions.rxnidx));
nes = numel(unique(tableSubsExceptions.rxnidx));
nr  = height(tableRxns);
nre = height(tableRxnExceptions);

fprintf('[%s] Unique reactions = %i, reaction exceptions = %4.1f%%, subsystem exceptions = %4.1f%%.\n', ...
    datestr(now), nr, 100*ner/nr, 100*nes/nr);

% Compile reaction identifier exceptions

fprintf('[%s] Consolidating reaction exceptions: %i\n', datestr(now),nre);
tableRxnExceptions = sortrows(tableRxnExceptions,[{'rxnidx'},options.fnIDs(2:end),{'modelidx'}]);
[TR, ~, idx] = unique(tableRxnExceptions(:, 2:end), 'rows', 'stable');
n    = height(TR);

nModels = nan(n,1);
for z = 1:n
    nModels(z) = sum(idx == z);
end
TR            = addvars(TR, nModels, 'Before', 1);
TR.fracModels = TR.nModels ./ tableRxns.nModels(TR.rxnidx);

% Transfer reaction identifiers according to threshold
fprintf('[%s] Transferring reaction identifiers (%i): ', datestr(now),height(TR));
nids    = length(options.fnIDs);
[nt,nd] = deal(0);

for z = 1:nids
    fn = options.fnIDs{z};
    if any(strcmp(TR.Properties.VariableNames, fn))
        idxf = ~(TR.(fn) == "");
        
        % transfer entries
        idx = find(idxf & TR.fracModels > options.TransferThreshold);
        for zr = 1:length(idx)
            tableRxns.(fn)(TR.rxnidx(idx(zr))) = TR.(fn)(idx(zr));
            nt = nt + 1;
        end
        
        % delete inconsistent entries
        idx = find(idxf & TR.fracModels <= options.TransferThreshold & ...
            TR.fracModels > (1 - options.TransferThreshold));
        for zr = 1:length(idx)
            tableRxns.(fn)(TR.rxnidx(idx(zr))) = "";
            nd = nd + 1;
        end
    end
end
fprintf('transferred %i, deleted %i.\n',nt,nd);

varargout{1} = tableRxnExceptions;
varargout{2} = TR;

return
