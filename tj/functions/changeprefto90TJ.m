function [new_preflist] = ...
    changeprefto90TJ(og_preflist)

new_preflist = double(og_preflist>90);
new_preflist(new_preflist>0) = 180;
new_preflist = abs(new_preflist-og_preflist);

end