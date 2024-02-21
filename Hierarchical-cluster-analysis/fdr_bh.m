function adj_p = fdr_bh(pvals, method)

s=size(pvals);
if (length(s)>2) || s(1)>1
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[~, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests

if strcmpi(method,'pdep')
    %BH procedure for independence or positive dependence
    wtd_p=m*p_sorted./(1:m);
    
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./(1:m));
    wtd_p=denom*p_sorted./[1:m];
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end


%compute adjusted p-values; This can be a bit computationally intensive
adj_p=zeros(1,m)*NaN;
[wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
nextfill = 1;
for k = 1 : m
    if wtd_p_sindex(k)>=nextfill
        adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
        nextfill = wtd_p_sindex(k)+1;
        if nextfill>m
            break;
        end
    end
end
adj_p=reshape(adj_p(unsort_ids),s);

end









