function Lu(du,u)
    #= Lu is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid. =#
    return du
end

function Lv(dv,v)
    #= Lv is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid. =#
    return dv
end

function LLu(du,u)
    #= LLu is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#
    return du
end

function LLv(dv,v)
    #= LLv is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#
    return dv
end
