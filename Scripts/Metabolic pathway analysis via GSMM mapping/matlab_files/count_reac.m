for i=1:830
    t = cell2mat(w_rxnSet_specie(i));
    jmax = max(t);
    jmin = min(t);
    if jmax>1
        jmax
    end
    if jmin<1
       % jmin
    end
end

