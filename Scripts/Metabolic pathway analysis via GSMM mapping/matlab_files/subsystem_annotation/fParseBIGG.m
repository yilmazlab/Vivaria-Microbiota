%% Function: Parse BIGG database files
function [tableRxns, varargout] = fParseBIGG(tableRxns,options)

% Load database
fprintf('[%s] Loading BIGG database.\n',datestr(now));
tableBIGG              = readtable(fullfile(options.dirData,options.fnBIGGReactions));
tableBIGG.bigg_id      = string(tableBIGG.bigg_id);
tableBIGG.old_bigg_ids = string(tableBIGG.old_bigg_ids);

n = height(tableBIGG);

f = tableBIGG.Properties.VariableNames;
for z = 1:length(options.fnIDs)
    if ~strcmp(f,options.fnIDs{z})
        tableBIGG.(options.fnIDs{z}) = strings(n,1);
    end
end

% Process identifiers
fprintf('[%s] Transferring BIGG annotations (%i): ',datestr(now),n);
nprog = floor(n/20);
[ns,nm,nk] = deal(0);

for z = 1:n
    s = tableBIGG.database_links{z};
    
    % MetaNetX
    sm = regexp(s,'MNXR\d+','match');
    if numel(sm) == 1
        tableBIGG.rxnMetaNetXID(z) = sm{1};
    end
    
    % SEED
    sm = regexp(s,'rxn\d{5}','match');
    if numel(sm) == 1
        tableBIGG.rxnSEEDID(z) = sm{1};
    end
    
    % KEGG
    sm = regexp(s,'kegg.reaction/R\d{5}','match');
    if numel(sm) == 1
        tableBIGG.rxnKEGGID(z) = strrep(sm{1},'kegg.reaction/','');
    end
    
    % Transfer data to reaction list
    s = tableBIGG.bigg_id(z) + "; " + tableBIGG.old_bigg_ids(z);
    s = unique(strtrim(regexp(s,';','split')));
    s = s(~(s==""));
    s = strrep(s, '_DASH_', '-');
    s = strrep(s, '_LPAREN_', '(');
    s = strrep(s, '_RPAREN_', ')');
    s = strrep(s, '__', '_');
           
    for zs = 1:length(z)
        idx = strcmp(tableRxns.rxn, s(zs));    
        if sum(idx) == 1
            tableRxns.rxnBiGGID(idx) = s(zs);
            if ~(tableBIGG.rxnMetaNetXID(z)=="")  ...
                    && (tableRxns.rxnMetaNetXID(idx)=="")
                tableRxns.rxnMetaNetXID(idx) = tableBIGG.rxnMetaNetXID(z);
                nm = nm+1;
            end
            if ~(tableBIGG.rxnSEEDID(z)=="")  ...
                    && (tableRxns.rxnSEEDID(idx)=="")
                tableRxns.rxnSEEDID(idx) = tableBIGG.rxnSEEDID(z);
                ns = ns+1;
            end
            if ~(tableBIGG.rxnKEGGID(z)=="")  ...
                    && (tableRxns.rxnKEGGID(idx)=="")
                tableRxns.rxnKEGGID(idx) = tableBIGG.rxnKEGGID(z);
                nk = nk+1;
            end
        end    
    end
    
    if ~mod(z,nprog)
        fprintf('.');
    end
end
fprintf('\n[%s] New annotations: MetaNetX = %i, KEGG = %i, SEED = %i\n',datestr(now),nm,nk,ns);

varargout{1} = tableBIGG;

return