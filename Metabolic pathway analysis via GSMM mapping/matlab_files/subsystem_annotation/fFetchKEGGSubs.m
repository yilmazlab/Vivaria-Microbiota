%% Function: Fetch KEGG reaction data for subsystem annotations

function [tableKEGGSubsOrg,varargout] = fFetchKEGGSubs(tableRxns, options, varargin)

% process optional inputs
tableKO = [];
if nargin == 3
    tableKO = varargin{1};
end

ncat      = length(options.fnKEGG);

% get reaction identifiers

idx = ~(tableRxns.rxnKEGGID=="");
ids = unique(tableRxns.rxnKEGGID(idx));
n   = length(ids);
nfail = 0;

tableKEGGSubsOrg = [];
tableKEGGSubs = table(strings(n,1),strings(n,1),cell(n,1), ...
    'VariableNames',{'rxnKEGGID','subsKEGG','ptrKO'});

fprintf('[%s] Fetching KEGG reactions (%i): ',datestr(now), n);
nprog = floor(n/20);

for zid = 1:n
    id = ids(zid);
    
    % get reaction record via REST request
    s = strcat('https://rest.kegg.jp/get/',id);
    try
        data = webread(s);
        data = splitlines(data);
        
        % get pathway identifiers
        pw  = cellfun(@(x) regexp(x,'rn\d{5}','match'), data, 'UniformOutput', false);
        idx = find(~cellfun(@isempty,pw));
        pw  = string(pw(idx));
        
        % compile pathway names
        nr    = length(idx);
        pname = strings(nr,1);
        [plevel,pptr] = deal(nan(nr,1));
        
        for z = 1:nr
            s = data{idx(z)};
            c = regexp(s,'(rn\d{5}\s+)','split');
            if length(c) == 2
                pname(z) = string(c{2});
                
                % map to KEGG ontology
                for zc = 1:ncat
                    idxo = strcmp(tableKO.(options.fnKO{zc}),pname(z));
                    if sum(idxo) == 1
                        plevel(z) = zc;
                        pptr(z)   = find(idxo);
                    end
                end                
            end
        end
        tableKEGGSubsOrg = [tableKEGGSubsOrg; table(repmat(string(id),nr,1),pw,pname,plevel,pptr, ...
            'VariableNames',{'rxnKEGGID','rxnKEGGSubsID','subsKEGG','levelKEGG','ptrKO'})];
        
        % compile list for subsystems
        sidx = find(plevel==3);
        ns = numel(sidx);
        s  = '';
        for zs = 1:ns
            t = tableKO(pptr(sidx(zs)),:);
            s = sprintf('%s|%s;%s;%s',s,t.(options.fnKO{1}),t.(options.fnKO{2}),t.(options.fnKO{3}));
        end
        tableKEGGSubs.rxnKEGGID(zid) = string(id);
        tableKEGGSubs.subsKEGG(zid)  = string(s(2:end));        
        tableKEGGSubs.ptrKO{zid}     = pptr(sidx);
        nfail = 0;
    catch
        fprintf('-');
        nfail = nfail + 1;
    end
    
    if nfail > 10
        error('KEGG REST API not functional - aborting!');
    end
    
    if ~mod(zid,nprog)
        fprintf('.');
    end
    
end
fprintf('\n');

varargout{1} = tableKEGGSubs;

return