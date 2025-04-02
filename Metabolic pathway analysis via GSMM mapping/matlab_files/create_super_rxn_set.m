function [new_rxn_set,new_w]=create_super_rxn_set(old_rxn_set,identifiers)

    un_ident=unique(identifiers);
    [~,ident_id]=ismember(identifiers,un_ident);

    new_rxn_set=cell(length(un_ident),1);
    new_w=cell(length(un_ident),1);

    for i=1:length(un_ident)
        modi=find(ident_id==i);
        if length(modi)==1
            new_rxn_seti=old_rxn_set{modi};
            new_wi=ones(length(new_rxn_seti),1);
        else
            [new_rxn_seti, sortedCount] = count_unique(vertcat(old_rxn_set{modi}));
            new_wi=sortedCount/length(modi);
        end
        new_rxn_set{i}=new_rxn_seti;
        new_w{i}=new_wi;
    end
end