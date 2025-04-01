%% Function: Processing of KEGG ontology

function [tableKO] = fParseKEGGOntology(options)

fprintf('[%s] Processing KEGG pathway ontology.\n',datestr(now));

KOText = fileread(fullfile(options.dirData,options.fnKEGGOntology));
KOText = splitlines(KOText);

% remove gene entries (level D)
fl     = cellfun(@isempty,regexp(KOText,'D\s'));
KOText = KOText(fl);

idx    = find(strcmp(KOText,'!'));
KOText = KOText(idx(1)+1:idx(2)-1);

fl     = cellfun(@length,KOText) > 1;
KOText = KOText(fl);

% parse ontology levels
n      = length(KOText);
ccurr  = cell(1,3);
ckegg  = cell(n,3);
idkegg = cell(n,1);
level  = cell(n,1);
ptr    = 1;

for z = 1:n
    s  = KOText{z};
    sl = double(s(1))-double('A')+1;
    level{ptr} = sl;
    
    m = regexp(s,'(\d{5})','match');
    idkegg{ptr} = ['rn' m{1}];
    
    switch sl
        case 1
            m        = regexp(s,'A\d+\s([\w\s]+)','tokens');
            ccurr(1) = m{1};
            ccurr{2} = '';
            ccurr{3} = '';
        case 2
            m        = regexp(s,'B\s+\d+\s([\w\s]+)','tokens');
            ccurr(2) = m{1};
            ccurr{3} = '';
        case 3
            m        = regexp(s,'C\s+\d+\s([^\[]+)','tokens');
            ccurr(3) = strtrim(m{1});
    end
    ckegg(ptr,:) = ccurr;
    ptr = ptr+1;
end

fns = {'ID','Level','Class','Subclass','Name'};
tableKO = cell2table([idkegg(1:ptr-1) level(1:ptr-1) ckegg(1:ptr-1,:)],'VariableNames', ...
    fns);
for z = 1:length(fns)
    tableKO.(fns{z}) = string(tableKO.(fns{z}));
end

return