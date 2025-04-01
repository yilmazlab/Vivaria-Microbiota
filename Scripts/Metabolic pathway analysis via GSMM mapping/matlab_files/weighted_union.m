function [new_list, weights]=weighted_union(obj, weights,un_obj)
list_obj=vertcat(obj{:});
list_weights=vertcat(weights{:});
res=accumarray(ismember_v2(list_obj,un_obj),list_weights,[length(un_obj),1],@sum);
new_list=un_obj;
weights=res;
end

function locb=ismember_v2(A,B)
    [~,locb]=ismember(A,B);
end

