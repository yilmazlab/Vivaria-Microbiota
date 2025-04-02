%% Function: Load and parse AGORA2 models

function [tableRxns,mapRxns,tableSubsExceptions,tableRxnExceptions,modelDir,varargout] = ...
    fParseModels(options,varargin)

%% Input arguments

if nargin == 1
    modelDir = dir(fullfile(options.dirModels,'*.mat'));
elseif nargin == 2
    modelDir = varargin{1};
else
    error('Invalid number of function arguents!');
end
nm = length(modelDir);

%% Parameters
tblock = 1000; % block size for tables
nids   = length(options.fnIDs);

%% Data structures
% standard annotation table (empty template)
s0 = strings(tblock,1);
n0 = zeros(tblock,1);
T0 = table(s0, s0, s0, s0, s0, s0, 1+n0, n0, n0, ...
    'VariableNames',[{'rxn'},options.fnIDs, ...
    {'subs','nModels','nMultiple','nEmpty'}]);

% exceptions table for subsystems (empty template)
TS0 = table(n0, s0, n0, s0, n0, ...
    'VariableNames',{'modelidx','rxn','rxnidx','subs','subsidx'});

% exceptions table for reaction identifiers (empty template)
TR0 = table(n0, s0, n0, s0, s0, s0, s0, ...
    'VariableNames',[{'modelidx','rxn','rxnidx'},options.fnIDs]);
TR0row = TR0(1,:);

M0 = false(tblock,nm);

%% Parse models for unique reactions and annotations
fprintf('\n[%s] Constructing unique reaction table:\n',datestr(now));

tableRxns           = T0;
mapRxns             = M0;
tableSubsExceptions = TS0;
tableRxnExceptions  = TR0;
[ptrT,ptrS,ptrR]    = deal(1);

for z = 1:nm
    
    if ~mod(z,20)
        nrxn = ptrT - 1;
        fprintf('\t%4i / %4i: %5i : %3.1f%% empty, %3.1f%% multiple, %3.1f%% MetaNetX, %3.1f%% SEED, %3.1f%% KEGG.\n', ...
            z, nm, nrxn, ...
            100*sum(tableRxns.nEmpty>0)/nrxn, 100*sum(tableRxns.nMultiple>0)/nrxn, ...
            100*sum(~(tableRxns.rxnMetaNetXID==""))/nrxn, ...
            100*sum(~(tableRxns.rxnSEEDID==""))/nrxn, ...
            100*sum(~(tableRxns.rxnKEGGID==""))/nrxn);
    end
    
    load(fullfile(options.dirModels,modelDir(z).name),'model');
    
    % flags for existence of identifier fields
    flfield = false(1,nids);
    for zf = 1:nids
        flfield(zf) = isfield(model,options.fnIDs{zf});
    end
    
    rxn  = string(model.rxns);
    subs = string(model.subSystems);
    fle  = double(cellfun(@isempty,model.subSystems));
       
    % handle matching reactions
    [~,idxM,idxM0]          = intersect(tableRxns.rxn,rxn);
    tableRxns.nModels(idxM) = tableRxns.nModels(idxM) + 1;
    
    % temp variables for efficiency
    trxns = cell(1,nids);
    for zf = 1:nids
        trxns{zf} = string(tableRxns.(options.fnIDs{zf}));
    end
    
    for zr = 1:length(idxM)
        idx  = idxM(zr);
        idx0 = idxM0(zr);
        mapRxns(idx,z) = true;
        
        % mismatches in subsystem annotations
        if ~strcmp(subs(idx0),tableRxns.subs(idx))
            tableRxns.nMultiple(idx) = tableRxns.nMultiple(idx) + 1;
            if fle(idx0)
                tableRxns.nEmpty(idx) = tableRxns.nEmpty(idx) + 1;
            end
            
            % collect exceptions
            if ptrS > height(tableSubsExceptions)
                tableSubsExceptions = vertcat(tableSubsExceptions,TS0);
            end
            tableSubsExceptions.modelidx(ptrS) = z;
            tableSubsExceptions.rxn(ptrS)      = rxn(idx0);
            tableSubsExceptions.rxnidx(ptrS)   = idx;
            tableSubsExceptions.subs(ptrS)     = subs(idx0);
            ptrS = ptrS+1;
        end
        
        % reaction identifiers
        flnew  = false;
        trtmp  = TR0row;
        
        for zf = 1:nids
            f = options.fnIDs{zf};
             
            if flfield(zf)
                rt = trxns{zf}(idx); %string(tableRxns.(f)(idx));
                rm = string(model.(f){idx0});
                
                if isempty(rt) || isempty(rm)
                elseif rt == "" && ~(rm == "") % new
                    tableRxns.(f)(idx) = rm;
                    
                elseif ~(rt == "") && ~(rm == "") && ~strcmp(rt,rm) % mismatch
                    % collect exceptions
                    if ~flnew
                        trtmp.modelidx = z;
                        trtmp.rxn      = rxn(idx0);
                        trtmp.rxnidx   = idx;
                    end
                    trtmp.(f) = rm;
                    flnew     = true;
                end
            end
        end
        if flnew
            if ptrR > height(tableRxnExceptions)
                tableRxnExceptions = vertcat(tableRxnExceptions,TR0);
            end
            tableRxnExceptions(ptrR,:) = trtmp;
            ptrR = ptrR + 1;
        end
    end
    
    % identify and handle new reactions
    [~,idxDiff] = setdiff(rxn,tableRxns.rxn);
    
    for idx0 = 1:length(idxDiff)
        if ptrT > height(tableRxns)
            tableRxns = vertcat(tableRxns,T0);
            mapRxns   = vertcat(mapRxns,M0);
        end
        tableRxns.rxn(ptrT)       = string(rxn(idxDiff(idx0)));
        tableRxns.subs(ptrT)      = string(subs(idxDiff(idx0)));
        tableRxns.nMultiple(ptrT) = 0;
        tableRxns.nEmpty(ptrT)    = fle(idx0);
        
        % additional identifiers
        for zf = 1:nids
            f = options.fnIDs{zf};
            if flfield(zf) && ~isempty(model.(f){idx0})
                tableRxns.(f)(ptrT) = string(model.(f){idx0});
            end
        end
        mapRxns(ptrT,z) = true;
        ptrT = ptrT+1;
    end
 
end

tableRxns           = tableRxns(1:ptrT-1,:);
tableSubsExceptions = tableSubsExceptions(1:ptrS-1,:);
tableRxnExceptions  = tableRxnExceptions(1:ptrR-1,:);

fprintf('\n');

nr = height(tableRxns);

varargout{1} = nr;

return