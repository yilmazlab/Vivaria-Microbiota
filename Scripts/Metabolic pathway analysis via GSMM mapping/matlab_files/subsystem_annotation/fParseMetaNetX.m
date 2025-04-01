%% Function: Parse MetaNetX cross-references

function [tableRxns, varargout] = fParseMetaNetX(tableRxns,options)

fprintf('[%s] Loading MetaNetX xref.\n',datestr(now));
tableMetaNetX = readtable(fullfile(options.dirData,options.fnMetanetxXRef), ...
    'Delimiter','\t','FileType','text');

% Select relevant xrefs
flkegg = contains(tableMetaNetX.source,'kegg.reaction');
flseed = contains(tableMetaNetX.source,'seed.reaction');
flbigg = contains(tableMetaNetX.source,'bigg.reaction');
fl     = flkegg | flbigg | flseed;
flkegg = flkegg(fl);
flseed = flseed(fl);
flbigg = flbigg(fl);
tableMetaNetX = tableMetaNetX(fl,:);

n     = height(tableRxns);
nprog = floor(n/20);
fprintf('[%s] Transferring MetaNetX annotations (%i): ',datestr(now),n);
[ns,nk,nm,ncorr] = deal(0);

for z = 1:n
    rx = tableRxns.rxnMetaNetXID(z);
    
    % search by BIGG identifier
    s = tableRxns.rxn(z);
    fl = strcmp("bigg.reaction:" + s, tableMetaNetX.source);
    if sum(fl & flbigg) == 1
        rxnew = string(tableMetaNetX.ID(fl));
        if ~(rx == "") && ~strcmp(rx, rxnew)
            if options.Verbose
                warning('Correcting MetaNetX identifier for <%s>: %s -> %s\n', ...
                    s, rx, rxnew);
            end
            if options.CorrectIdentifiers
                tableRxns.orgMetaNetXID(z) = rx;
                rx    = rxnew;
                ncorr = ncorr + 1;
            end
        elseif (rx == "")
            rx = rxnew;
            nm = nm+1;
        end
    elseif ~(rx == "") && options.DeleteIdentifiers
        idx = strcmp(tableMetaNetX.ID, rx);
        if ~any(idx) % remove identifier
            tableRxns.orgMetaNetXID(z) = rx;
            rx    = "";
            tableRxns.rxnMetaNetXID(z) = rx;
            ncorr = ncorr + 1;
        end
    end

    % search by SEED identifier
    if rx == ""
        fl = strcmp("seed.reaction:" + s,tableMetaNetX.source);
        if sum(fl & flseed) == 1
            rx = string(tableMetaNetX.ID(fl));
            nm = nm+1;
        end
    end
    
    % gather data based on MetaNetX identifier
    if ~(rx == "")
        tableRxns.rxnMetaNetXID(z) = rx;
        
        fl = strcmp(rx,tableMetaNetX.ID);
        if any(fl & flkegg)
            s = tableMetaNetX.source(fl & flkegg);
            s = regexp(s,'\w\d{5}','match');
            if tableRxns.rxnKEGGID(z) == ""
                tableRxns.rxnKEGGID(z) = string(s{1});
                nk = nk+1;
            end
        end
        if any(fl & flseed)                     
            s = tableMetaNetX.source(fl & flseed);
            s = regexp(s,'rxn\d{5}','match');
            
            % only unique mappings
            if sum(fl & flseed) == 1 
                if tableRxns.rxnSEEDID(z) == ""
                    tableRxns.rxnSEEDID(z) = string(s{1});
                    ns = ns+1;
                end
                
                % extract KEGG Id from SEED description
                sraw = tableMetaNetX.description(fl & flseed);
                s = regexp(sraw,'R\d{5}','match');
                if ~isempty(s{1}) && tableRxns.rxnKEGGID(z) == ""
                    tableRxns.rxnKEGGID(z) = string(s{1});
                    nk = nk+1;
                end
            end
        end        
    end
    
    if ~mod(z,nprog)
        fprintf('.');
    end
end
fprintf('\n[%s] New annotations: MetaNetX = %i (%i corrected), KEGG = %i, SEED = %i\n', ...
    datestr(now),nm,ncorr,nk,ns);

varargout{1} = tableMetaNetX;

return