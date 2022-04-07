function finitefind(a, val, find)
    b = deepcopy(a)
    b[.!isfinite.(a)] .= val
    return find(b)
end

finitefindmin(a) = finitefind(a, Inf, findmin)
finitefindmax(a) = finitefind(a, -Inf, findmax)

