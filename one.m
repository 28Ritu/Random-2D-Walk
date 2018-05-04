function u=one(m, L)
    if (m>=L)
        u=L;      % reflective boundary 
    else
        u=m+1;
    end
end