%% Annotate reactions for AGORA2 models via KEGG / SEED / BiGG / MetaNetX

%% Definitions

options = [];

options.dirData      = '../../data/subsystem_assignation/'; % path to external data directory
options.dirModels    = '../../data/AGORA2_models/mat/'; % path to AGORA2 models

% SEED reactions including annotations
options.urlSEEDReactions = 'https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/reactions.tsv';
options.fnSEEDReactions  = 'SEED_reactions.tsv';

% KEGG pathway ontology and reaction list
options.urlKEGGOntology  = 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=';
options.fnKEGGOntology   = 'kegg-ontology.txt';
options.urlKEGGAllRxns   = 'https://rest.kegg.jp/list/reaction';
options.fnKEGGAllRxns    = 'kegg-reaction.txt';

% BIGG reactions
options.urlBIGGReactions = 'http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt';
options.fnBIGGReactions  = 'bigg_models_reactions.txt';

% MetaNetX xref
options.urlMetanetxXRef  = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv';
options.fnMetanetxXRef   = 'metanetx-rxns-xref.tsv';
options.procMetanetxXRef = true;

% Variable names for column identifiers
options.fnIDs  = {'rxnBiGGID','rxnMetaNetXID','rxnSEEDID','rxnKEGGID'}; % reaction identifiers
options.fnKEGG = {'subsKEGG1','subsKEGG2','subsKEGG3'}; % KEGG subsystem hierarchy (output)
options.fnKO   = {'Class','Subclass','Name'}; % KEGG subsystem hierarchy (original)

% Annotation transfer from AGORA2 models
options.TransferThreshold  = 0.9; % threshold for consistency on annotations for transfer from models
options.UniqueKEGGSubs     = true; % annotate KEGG subsystems only when unique
options.ReplaceKEGGSubs    = true; % overwrite SEED-derived subsystem annotations
options.CorrectIdentifiers = true; % correct reaction identifiers transferred from models
options.DeleteIdentifiers  = false; % delete reaction identifiers transferred from models

% Variables to save
options.saveVars = {'options','tableRxns','mapRxns','tableSubsExceptions','tableRxnExceptions', ...
    'modelDir','nr','tableKEGGSubs'};
options.fnSave   = 'tmp-annotation-agora2.mat';
options.flagSave = true;
options.Verbose  = false; % verbose output

%% Main script

% Set up environment
[exitflag] = fSetEnvironment(options);
if ~exitflag
    error('Data sources are incomplete!');
end
    
% Load and parse models
if ~exist('tableRxns','var')
    [tableRxns,mapRxns,tableSubsExceptions,tableRxnExceptions,modelDir,nr] = fParseModels(options); 
end

% Prepare data structures (identifier columns)
n = height(tableRxns);
f = tableRxns.Properties.VariableNames;
fn = [options.fnIDs strrep(options.fnIDs,'rxn','org')];
for z = 1:length(fn)
    if ~strcmp(f,fn{z})
        tableRxns.(fn{z}) = strings(n,1);
    end
end

% Consolidate model annotations
[tableRxns,~,TR] = fConsolidateModels(tableRxns,tableSubsExceptions,tableRxnExceptions,options);

% Get KEGG pathway ontology
[tableKO] = fParseKEGGOntology(options);

% Map reaction identifiers via MetaNetX / BIGG / SEED cross-refs
fAnnotationStats(tableRxns);
[tableRxns,tableBIGG]     = fParseBIGG(tableRxns,options);
fAnnotationStats(tableRxns);
[tableRxns,tableMetaNetX] = fParseMetaNetX(tableRxns,options);
fAnnotationStats(tableRxns);
[tableRxns,tableSEED]     = fParseSEED(tableRxns,tableKO,options);
fAnnotationStats(tableRxns);

% Annotate with KEGG subsystems
if ~exist('tableKEGGSubs', 'var')
    [~,tableKEGGSubs] = fFetchKEGGSubs(tableRxns,options,tableKO);
end
tableRxns = fAnnotateKEGGSubs(tableRxns,tableKO,tableKEGGSubs,options);

% Save analysis and annotation results
if options.flagSave
    save(options.fnSave,options.saveVars{:});
    writetable(tableRxns, "../../data/subsystem_assignation/subsystem_mapping.csv")
end
