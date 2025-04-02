%% Function: Processing of SEED reaction file annotations / addition to reactionAnnotation
function [tableRxns, varargout] = fParseSEED(tableRxns,tableKO,options,varargin)

if nargin == 4
    tableSEED = varargin{1};
    
else % Compile SEED reactions
    
    fprintf('[%s] Loading SEED reaction annotations.\n',datestr(now));
    tableSEED = readtable(fullfile(options.dirData,options.fnSEEDReactions), ...
        'Delimiter','\t','FileType','text');
    
    % Prepare data structures
    
    n = height(tableSEED);
    [tableSEED.MetaCycID,tableSEED.MetaCycName,tableSEED.MetaCycLevel] = deal(cell(n,1));
    fn = tableRxns.Properties.VariableNames;
    
    ncat = length(options.fnKEGG);
    for z = 1:ncat
        tableSEED.(options.fnKEGG{z}) = cell(n,1);
        if ~any(strcmp(fn,options.fnKEGG{z}))
            tableRxns.(options.fnKEGG{z}) = strings(height(tableRxns),1);
        end
    end
    tableSEED.KEGGPtr   = cell(n,1);
    tableSEED.rxnKEGGID = strings(n,1);
    tableSEED.rxnBiGGID = strings(n,1);
    tableSEED.aliases   = string(tableSEED.aliases);
    tableSEED.pathways  = string(tableSEED.pathways);
    
    % Process SEED annotations
    
    fprintf('[%s] Processing SEED annotations (%i): ',datestr(now),n);
    nprog = floor(n/20);
    
    UKEGGPtr = [];
    
    for z = 1:n
        % extract KEGG and BiGG identifiers
        s = tableSEED.aliases(z);
        sm = string(regexp(s,'KEGG: (R\d{5})','tokens'));
        if ~isempty(sm)
          tableSEED.rxnKEGGID(z) = sm;
        end
        sm = string(regexp(s,'BiGG: ([^\|]+)','tokens'));
        if ~isempty(sm)
          tableSEED.rxnBiGGID(z) = sm;
        end
                        
        % extract pathway annotations
        s = tableSEED.pathways(z) + ";";
        c = split(s,'|');
        for zc = 1:length(c)
            pw = regexp(c{zc},'(\s[^;]+;)','tokens');
            pw(cellfun(@isempty,pw)) = {[""]};
            pw =  strtrim(strrep(string(pw'),';',''));
            
            pwname = regexp(pw,'\((.+)\)','tokens');
            pwname(cellfun(@isempty,pwname)) = {[""]};
            pwname = string(pwname);
            
            pwid   = regexp(pw,'([^\s]+)\s\(','tokens','once');
            pwid(cellfun(@isempty,pwid)) = {[""]};
            pwid = string(pwid);
            
            if contains(c{zc},'MetaCyc')
                tableSEED.MetaCycID{z}   = pwid;
                tableSEED.MetaCycName{z} = pwname;
                
            elseif contains(c{zc},'KEGG')
                np    = length(pwid);
                pwptr = nan(np,1);
                
                for zd = 1:np
                    idx = find(strcmp(tableKO.Name,pwname{zd}));
                    if length(idx) == 1
                        pwptr(zd) = idx;
                        UKEGGPtr  = union(UKEGGPtr,idx);
                    end
                end
                
                fl = ~isnan(pwptr); % remove undefined / higher-level categories
                tableSEED.KEGGPtr{z}  = pwptr(fl);
                tableSEED.KEGGID{z}   = pwid(fl);
                tableSEED.(options.fnKEGG{3}){z} = pwname(fl);
                tableSEED.(options.fnKEGG{2}){z} = unique(string(tableKO.Subclass(pwptr(fl))));
                tableSEED.(options.fnKEGG{1}){z} = unique(string(tableKO.Class(pwptr(fl))));
            end
        end
        
        if ~mod(z,nprog)
            fprintf('.');
        end
    end
    
    fprintf('\n');
end

%% Transfer annotations to model table

tableSEED.id  = string(tableSEED.id);
n = height(tableRxns);
fprintf('[%s] Transferring SEED annotations (%i): ',datestr(now),n);
nprog = floor(n/20);

[ns,nm,nk] = deal(0);
for z = 1:n
    idx = strcmp(tableRxns.rxn(z), tableSEED.rxnBiGGID) | ...
        (~(tableRxns.rxnSEEDID(z) == "") & strcmp(tableRxns.rxnSEEDID(z), tableSEED.id));
    
    if sum(idx) == 1
        if ~(tableSEED.id(idx)=="")  ...
                && (tableRxns.rxnSEEDID(z)=="")
            tableRxns.rxnSEEDID(z) = tableSEED.id(idx);
            ns = ns+1;
        end
        if ~(tableSEED.rxnKEGGID(idx)=="")  ...
                && (tableRxns.rxnKEGGID(z)=="")
            tableRxns.rxnKEGGID(z) = tableSEED.rxnKEGGID(idx);
            nk = nk+1;
        end
        
        % Transfer KEGG subsystem annotations
        for zs = 1:length(options.fnKEGG)
            if numel(tableSEED.(options.fnKEGG{zs}){idx}) == 1
                tableRxns.(options.fnKEGG{zs})(z) = string(tableSEED.(options.fnKEGG{zs}){idx});
            end
        end              
    end
    
    if ~mod(z,nprog)
        fprintf('.');
    end
    
    
end
fprintf('\n[%s] New reaction identifiers: MetaNetX = %i, KEGG = %i, SEED = %i\n',datestr(now),nm,nk,ns);

varargout{1} = tableSEED;

return