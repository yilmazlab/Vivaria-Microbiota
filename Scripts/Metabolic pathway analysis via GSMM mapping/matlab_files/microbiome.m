clc; clear; close all;

% This code is adapted from Swiss IBD Cohort Investigators *et al.* , 2019. DOI: 10.1038/s41591-018-0308-z
%% Choose cohort and load data

fprintf('Reading otu data \n');

load('..\data\processed_files\otumat.mat')
% To avoid rewriting the code
otumat = struct2table(otumat);
otumat = table2array(otumat);

fprintf('Otu data read \n');

fprintf('Reading metadata \n');

opts = detectImportOptions('..\data\processed_files\taxonomy.csv');
tax = readtable('..\data\processed_files\taxonomy.csv');

model_class_family=readtable('../data/AGORA2_models/AGORA2_TableS1.xlsx');

%% Compute the relative abundances

otumat = otumat ./ sum(otumat);

%% Load models

% Load the same models as in the previous collaboration
model_dir = '../data/AGORA2_models/mat/'; 
model_files = dir(strcat(model_dir,'*.mat'));
models = struct;
n = length(model_files);
h = waitbar(0,'Loading models...');
for i=1:n
    waitbar(i/n,h);
    % Load and unpack model
    model = load(strcat(model_dir,model_files(i).name),'model');
    f = fieldnames(model);
    model = model.(f{1});
    % Save to struct
    models(i).rxns = model.rxns;
    models(i).rxnNames = model.rxnNames;
    models(i).subSystems = model.subSystems;
    models(i).name = model_files(i).name(1:end-4);
end
h.delete;


%% Check unique OTUs

taxcat=strcat(tax.Kingdom,tax.Phylum,tax.Class,tax.Order,tax.Family,tax.Genus,tax.Species);
fprintf('Unique OTUs: %d\n',length(unique(taxcat))) % The OTUs are obtained by concatenating 7 levels of the classification

%% Get models and reactions at species, genus, and family level

% Create rxn set and associated weigths for each subspecie
rxnSet_specie=arrayfun(@(s) s.rxns, models,'uni',0);
weights_specie=cellfun(@(c) ones(length(c),1),rxnSet_specie,'uni',0);

% Extract the model classification 
model_files = dir(strcat(model_dir,'*.mat')); % Repeat to correct for mistake
model_names=extractfield(model_files,'name'); 
model_names=strrep(model_names,'_','-');
model_class_sp=regexp(model_names,'\<(\w+)-(\w+)([a-zA-Z_0-9-]+).mat','tokens');
model_class_sp=[model_class_sp{:}];

% model_class is a model classification with the first column being the
% genus and second column the specie
model_class1=cellfun(@(c) c(1),model_class_sp,'uni',0)';
model_class2=cellfun(@(c) c(2),model_class_sp,'uni',0)';
model_class_sp=[[model_class1{:}]',[model_class2{:}]'];

model_class_cat=strcat(model_class_sp(:,1),'_',model_class_sp(:,2));


% add another column family in model_class_sp


model_class_family.organism=strrep(model_class_family.MicrobeID,' ','_');
model_class_family.organism=strrep(model_class_family.MicrobeID,'-','_');
model_class_family.organism=strrep(model_class_family.MicrobeID,'sp.','sp');
model_class_family.organism=strrep(model_class_family.MicrobeID,'str.','str');
model_class_family.organism=strrep(model_class_family.MicrobeID,'subsp.','subsp');
model_class_family.organism=strrep(model_class_family.MicrobeID,'.','_');
model_class_family.organism=strrep(model_class_family.MicrobeID,'/','_');
model_class_family.organism=strrep(model_class_family.MicrobeID,'__','_');

model_class_family.organism=strrep(model_class_family.MicrobeID,',','');
model_names2=strrep(model_names,'-','_');
model_names2=strrep(model_names2,'.mat','');
[~,inClass]=ismember(model_names2,model_class_family.MicrobeID);

%0 means that it is not a member
model_class_sp(:,3)=model_class_family.Family(inClass);

%un_genus is the unique set of genus represented in the model collection

[un_genus,iungenus]=unique(model_class_sp(:,1));
[rxnSet_genus,w_rxnSet_genus]=create_super_rxn_set(rxnSet_specie,model_class_sp(:,1));

%un_family is the unique set of families represented in the model collection

un_family=unique(model_class_sp(:,3));
[rxnSet_family,w_rxnSet_family]=create_super_rxn_set(rxnSet_specie,model_class_sp(:,3)); 
rxnSet=[rxnSet_specie';rxnSet_genus;rxnSet_family];
w_rxnSet_specie=cellfun(@(c) ones(length(c),1),rxnSet_specie,'un',0);
w=[w_rxnSet_specie';w_rxnSet_genus;w_rxnSet_family];

nsp=length(models);
model_class_all=model_class_sp;
model_class_all(nsp+1:nsp+length(un_genus),1)=un_genus;
model_class_all(nsp+1:nsp+length(un_genus),3)=model_class_sp(iungenus,3); % family of corresponding genus 
model_class_all(end+1:end+length(un_family),3)=un_family;


%% Remove rows not mapped to family level or below from OTU matrix

sp_barcode=[tax.Genus,tax.Species,tax.Family]; % Rank6, Rank 7, Rank5
% format the species in the data 
sp_barcode(:,1)=strrep(sp_barcode(:,1),'g__','');
sp_barcode(:,2)=strrep(sp_barcode(:,2),'s__','');
sp_barcode(:,1)=strrep(sp_barcode(:,1),'NA','');
sp_barcode(:,2)=strrep(sp_barcode(:,2),'NA','');
sp_barcode(:,3)=strrep(sp_barcode(:,3),'NA','');
sp_barcode(:,3)=strrep(sp_barcode(:,3),'f__','');

% find the unique genus specie combination
sp_barcode_cat=strcat(sp_barcode(:,3),'-',sp_barcode(:,1),'-',sp_barcode(:,2));
[unSp_data,unSp_data_id,~]=unique(sp_barcode_cat);
[~,spid]=ismember(sp_barcode_cat,unSp_data);
otumat_red=cell2mat(arrayfun(@(x) accumarray(spid,...
 otumat(:,x), [length(unSp_data) 1], @sum,NaN),1:size(otumat,2),'un',0));
% Sum abundances for OTUs with the same taxonomy


not_mapped_counts = otumat_red(strcmp(unSp_data,'--'),:)./sum(otumat_red);


%% Remove rows not mapped to model collection from OTU matrix

% 0s correspond to taxonomies that are not in any of the models
model_class_cat=strcat(model_class_all(:,3),'-',model_class_all(:,1),'-',model_class_all(:,2));
[~, posinmodcoc]=ismember(unSp_data,model_class_cat); % posinmodcoc constrains the locations of membership

%visualize mapped species
mapped_species=[unSp_data, num2cell(posinmodcoc)];

% percentage of counts not mapped to model collection
not_mapped_counts=sum(otumat_red(posinmodcoc==0,:))./sum(otumat_red);


%%
% otumat reduced to mapped species in model collection
otumat_red_mapped=otumat_red(posinmodcoc~=0,:);

%% Build reaction matrix

rxnSet_red=rxnSet(posinmodcoc(posinmodcoc~=0)); % corresponding rxnset
un_rxnSet_red=unique(vertcat(rxnSet_red{:}));
w_red=w(posinmodcoc(posinmodcoc~=0)); % corresponding weights


n = size(otumat_red,2);
abund_rxns=zeros(length(un_rxnSet_red),n);
h = waitbar(0,'Building reaction matrix...');

for i=1:n
    waitbar(i/n,h);
    abund=cellfun(@(w,f) w*f,w_red(otumat_red_mapped(:,i)~=0),num2cell(otumat_red_mapped(otumat_red_mapped(:,i)~=0,i)),'un',0);
    [~, weights]=weighted_union(rxnSet_red(otumat_red_mapped(:,i)~=0), ...
    abund,un_rxnSet_red);
    abund_rxns(:,i)=weights;
end
h.delete;

% removing exchange and biomass reactions
exc=startsWith(un_rxnSet_red,'EX_');
bio=startsWith(un_rxnSet_red,'biomass');
demand=startsWith(un_rxnSet_red,'DM_');
sink=startsWith(un_rxnSet_red,'sink_');
otherdemand=contains(un_rxnSet_red,{'dreplication','pbiosynthesis','rtranscription'});
un_rxnSet_red=un_rxnSet_red((exc+bio+demand+sink+otherdemand)==0);
abund_rxns=abund_rxns((exc+bio+demand+sink+otherdemand)==0,:);

%% Get reaction names and subsystems

rxnnames = cell(length(un_rxnSet_red),1);
rxn_fullnames = cell(length(un_rxnSet_red),1);
subsys = repmat({""},length(un_rxnSet_red),1);
unique_subsys = repmat({''},length(un_rxnSet_red),1);
n = length(models)
h = waitbar(0,'Getting reaction names and subsystems...');
% Count how many useful reactions appear in each model
model_count =zeros(n,1);
for i=1:n
    waitbar(i/n,h);
    model = models(i);
    for j=1:length(model.rxns)
        rxnidx = find(strcmp(un_rxnSet_red,model.rxns(j)));
        if ~isempty(rxnidx)
            rxnnames(rxnidx) = model.rxns(j); 
            rxn_fullnames(rxnidx) = model.rxnNames(j);
            subsys{rxnidx} = [subsys{rxnidx},convertCharsToStrings(model.subSystems{j})];
        end
    end
end

h.delete;
rxnnames= vertcat(rxnnames);

%%
subsys2 = repmat({""},length(un_rxnSet_red),1);
for i=1:length(subsys2)
    isubsys=subsys{i};
    [val,~,index] = unique(isubsys);
    if length(val)>2 % Empty string + 1 subsystem
        i
    end
    counts = accumarray(index,1);
    val = val.transpose();
    value_counts = [val, counts];
    [~,i3]=max(counts);
    subsys2{i}=val(i3,1);
end

subsys=subsys2;

%% Dfs of interest

% Save the name of the species we mapped to an OTU 
mapped_species_c = mapped_species(posinmodcoc ~= 0);
mapped_species_c = cell2table(mapped_species);


% Save the matrix of reaction abundances
abund_rxns_c = array2table(abund_rxns);
abund_rxns_c.full_names = rxn_fullnames;
abund_rxns_c.name = rxnnames;
abund_rxns_c.subsystems = subsys;

% Normalized reactions abundance:
writetable(abund_rxns_c, "..\data\processed_files\normalized_reaction_abundance.csv") 
% Species-model mapping used for computing the normalized reaction
% abundance. This mapping is also used to generate Extended Figure 6d (see
% python_files/exploratory_analyses/mapped_species_AGORA2)
writetable(mapped_species_c, "..\data\processed_files\mapped_species_matlab.csv") 