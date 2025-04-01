%% Function: Annotate KEGG pathways based on KEGG reaction identifiers

function [tableRxns] = fAnnotateKEGGSubs(tableRxns,tableKO,tableKEGGSubs,options)

ncat             = length(options.fnKEGG);
tableKO.ID       = string(tableKO.ID);

ns    = 0;
nk    = height(tableKEGGSubs);
nsubs = cellfun(@numel, tableKEGGSubs.ptrKO) == 1; % unique annotation

fprintf('[%s] Transferring KEGG subsystems (%i): ',datestr(now),nk);
nprog = floor(nk/20);

% add / overwrite subsystem columns if necessary
for zc = 1:ncat
    if ~any(strcmp(tableRxns.Properties.VariableNames,options.fnKEGG{zc})) || ...
            options.ReplaceKEGGSubs
        tableRxns.(options.fnKEGG{zc}) = strings(height(tableRxns),1);
    end
end
tableRxns.subsKEGG = strings(height(tableRxns),1);

for z = 1:nk
    rxn = tableKEGGSubs.rxnKEGGID(z);
    ptr = tableKEGGSubs.ptrKO{z};
    
    if length(ptr) > 1  % filter for metabolism
        subM = tableKO.(options.fnKO{1})(ptr);
        ptr  = ptr(contains(subM, 'Metabolism'));
    end
    
    if ~isempty(ptr)
        idxr = strcmp(tableRxns.rxnKEGGID, rxn);                
        if any(idxr)
            tableRxns.subsKEGG(idxr) = tableKEGGSubs.subsKEGG(z);
            if length(ptr) == 1
                for zc = 1:ncat
                    tableRxns.(options.fnKEGG{zc})(idxr) = tableKO.(options.fnKO{zc})(ptr);
                end
                ns = ns + 1;
            end
        end
    end
    
    if ~mod(z,nprog)
        fprintf('.');
    end
end

fprintf('\n');
fprintf('[%s] Transferred KEGG subsystems: %i / %i\n',datestr(now),ns,nk);

return