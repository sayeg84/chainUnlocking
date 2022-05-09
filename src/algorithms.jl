include("intersections.jl")

using Distributions    # need for truncated gaussian distribution
using IterativeSolvers # julia's `\` operator is not good for solving this systems
using Random           # generating random permutations for multi mutation steps


function distToLine(Q::AbstractChain,P::AbstractChain=makeLine(Q))
    if length(P) == length(Q)
        return overlapedRmsd(Q,P)^2
    else
        error("Chains must be of same length")
    end
end

function distToFlat(Q::AbstractChain,P::AbstractChain=flatten(Q))
    if length(P) == length(Q)
        return overlapedRmsd(Q,P)^2
    else
        error("Chains must be of same length")
    end
end

function squaredMaxSpan(Q::AbstractChain)
    v = Q[1] - Q[end]
    return -dot(v,v)
end

function squaredMaxSpanNormed(Q::AbstractChain)
    v = Q[1] - Q[end]
    return -dot(v,v)/totalLength(Q)
end

function fourKnotSpan(Q::AbstractChain)
    v = Q[2] - Q[end-1]
    return -dot(v,v)
end

function demaineEnergySimple(Q::AbstractChain)
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in i+2:(n+1)
            sum += 1/pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
        end
    end
    return sum
end

function demaineEnergy(Q::AbstractChain)
    n = length(Q)
    sum = 0
    for i in 1:n
        #=
        All of these proceedures with two four loops grouped under a 
        single for loop using an `if o looping over vcat(1:(i-1),(i+2):(n+1))

        Both of these options are discarded as suboptimal for performance: 
        using the `if` makes a comparison over the index in each step of the loop,
        and using the `vcat` allocates memory instead of using a lazy range.
        =#
        for j in 1:(i-1)
            sum += 1/pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
        end
        for j in (i+2):(n+1)
            sum += 1/pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
        end
    end
    return sum
end

function tangentPointKernel(p::Point,q::Point,tang::Point,alpha::Real,beta::Real)
    dir = p-q
    return norm(cross(tang,dir))^alpha/norm(dir)^beta
end

function tangentPointDiscreteKernel(Q::AbstractChain,i::Integer,j::Integer,
                            alpha::Real,beta::Real)
    sum = 0.0
    Ti = unitVector(Q[i+1] - Q[i])
    for a in 0:1, b in 0:1
        sum += tangentPointKernel(Q[i+a],Q[j+b],Ti,alpha,beta)
    end
    return sum/4
end

function tangentEnergySimple(Q::AbstractChain;alpha::Real=3,beta::Real=6)
    n = length(Q)
    sum = 0
    for i in 1:n
        li = distance(Q[i+1],Q[i])
        for j in (i+2):n
            lj = distance(Q[j+1],Q[j])
            sum += tangentPointDiscreteKernel(Q,i,j,alpha,beta)*li*lj
        end
    end
    return sum
end

function tangentEnergy(Q::AbstractChain;alpha::Real=3,beta::Real=6)
    n = length(Q)
    sum = 0
    for i in 1:n
        li = distance(Q[i+1],Q[i])
        # see `demaineEnergy` comment
        for j in 1:(i-2)
            lj = distance(Q[j+1],Q[j])
            sum += tangentPointDiscreteKernel(Q,i,j,alpha,beta)*li*lj
        end
        for j in (i+2):n
            lj = distance(Q[j+1],Q[j])
            sum += tangentPointDiscreteKernel(Q,i,j,alpha,beta)*li*lj
        end
    end
    return sum
end

function tangentEnergyFrac(Q::AbstractChain)
    return tangentEnergy(Q,alpha=2,beta=4.5)
end

function tangentEnergyFracNormed(Q::AbstractChain)
    return tangentEnergyFrac(Q)/totalLength(Q)
end

function needleMixture(Q::AbstractChain)
    return tangentEnergyFrac(Q) + squaredMaxSpan(Q)
end

function fourKnotMixture(Q::AbstractChain)
    return tangentEnergyFrac(Q) + fourKnotSpan(Q)
end


function uniformAngle(angmin::Real,angmax::Real)
    return rand()*(angmax-angmin) + angmin
end

# bad practice as it has to create a new Distribution.TruncatedNormal object for each call
function gaussianAngle(angmin::Real,angmax::Real)
    rand(Distributions.TruncatedNormal(0,1,angmin,angmax))
end

# we can use a single integer `k` to encode changes in both dihedral and internal ang_vals
# if the integer is of size 1 <= k <= n-3, then it representes a change in dihedral angle
function generateAngleIndex(n::Integer,internal::Bool)::Integer
    return internal ? rand(1:(2*n-3)) : rand(1:(n-2))
end

function checkSingleMutation(P::AbstractChain,ang_idx::Integer,alpha::Real;debug::Bool=false)
    n = length(P)
    intersection = false
    if 1 <= ang_idx <= n-2
        newP = moveBeforeDihedral(P,ang_idx)
        intersection = checkRotationIntersection(newP,ang_idx,alpha,false,debug=debug)
        if !intersection; dihedralRotateFast!(newP,ang_idx,alpha); end;
    elseif n-1 <= ang_idx <= 2*n-3
        ang_idx = ang_idx - n + 2 # adjusting to be in range 
        newP = moveBeforeInternal(P,ang_idx)
        intersection = checkRotationIntersection(newP,ang_idx,alpha,true,debug=debug)
        if !intersection; internalRotateFast!(newP,ang_idx,alpha); end;
    else
        error("index $(ang_idx) is not supported")
    end
    return intersection,newP
end

function checkMultipleMutation(P::AbstractChain,mov_ang_idxs::Array{<:Integer,1},alphas::Array{<:Real,1};debug::Bool=false)
    n = length(P)
    intersection = false
    m1 = length(mov_ang_idxs)
    newP = P
    i = 1
    while i <= m1 && !intersection
        loc_inter,newP = checkSingleMutation(newP,mov_ang_idxs[i],alphas[i],debug=debug)
        intersection = intersection || loc_inter
        i += 1
    end
    return intersection, newP
end

function localRandomSearchStep(Q,newQ,inter_flag,c,minf_val,minf_newval,
                               fun_vals,ang_idx,ang_idxs,alpha,ang_vals)
    if !inter_flag && minf_newval < minf_val
        ang_idxs[c] = ang_idx
        ang_vals[c] = alpha
        fun_vals[c] = minf_newval
        return minf_newval, newQ
    else
        fun_vals[c] = minf_val
        return minf_val, Q
    end
end

function basicSMetaheuristic(Q::PolygonalChain,minFunc::Function,
                             tolerance::Real,phimax::Real,
                             phimin::Real,temp_init::Real,
                             max_iter::Integer,
                             advanceFunc!::Function)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end

    ang_idxs = zeros(Int16,max_iter)
    ang_vals = zeros(eltype(Q),max_iter)
    fun_vals = zeros(eltype(Q),max_iter)
    nq = length(Q)
    minf_val = minFunc(Q)
    c = 1
    while minf_val > tolerance && c <= max_iter
        phi = rand()*(phimax-phimin) + phimin
        ang_idx = rand(1:(nq-2))
        newQ = moveBeforeDihedral(Q,ang_idx)
        inter_flag = checkRotationIntersection(newQ,ang_idx,phi)
        newQ = dihedralRotate(newQ,ang_idx,phi)
        minf_newval = minFunc(newQ)
        minf_val,Q = advanceFunc!(Q,newQ,inter_flag,c,minf_val,minf_newval,
                                  fun_vals,ang_idx,ang_idxs,phi,ang_vals)
        c += 1
    end
    return Q,ang_vals[1:(c-1)],ang_idxs[1:(c-1)],fun_vals[1:(c-1)]
end



function localRandomSearch2(Q::PolygonalChain,minFunc::Function,tolerance::Real=1e-2,phimax::Real=pi/2,phimin::Real=-pi/2;temp_init=1,max_iter::Integer=1000)
    return basicSMetaheuristic(Q,minFunc,
    tolerance,phimax,
    phimin,temp_init,
    max_iter,localRandomSearchStep)
end

# `temp_init` argument is useless and only put for compatibility reasons
function localRandomSearch(Q::PolygonalChain,
    minFunc::F,
    tolerance::Real=1e-2,
    alphamax::Real=pi/2,
    alphamin::Real=-pi/2,
    internal::Bool=false;
    population::Integer=8,
    max_iter::Integer=1000,
    debug::Bool=false) where {F}
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    ang_idxs  = zeros(Int16,max_iter)
    ang_vals  = zeros(ftype(Q),max_iter)
    fun_vals =  zeros(ftype(Q),max_iter)
    nq = length(Q)
    fval = minFunc(Q)
    c = 1
    while fval > tolerance && c <= max_iter
        alpha = rand()*(alphamax-alphamin) + alphamin
        (internal) && (alpha = alpha/2)
        ang_idx = generateAngleIndex(nq,internal)
        inter_flag , newQ = checkSingleMutation(Q,ang_idx,alpha)
        if !inter_flag
            fvalnew = minFunc(newQ)
            if fvalnew < fval
                Q = newQ
                fval = fvalnew
                ang_idxs[c] = ang_idx
                ang_vals[c] = alpha
                fun_vals[c] = fvalnew
            end
        else
            fun_vals[c] = fval
        end
        c += 1
    end
    c = c > max_iter ? max_iter : c
    return Q,ang_vals[1:c],ang_idxs[1:c],fun_vals[1:c]
end

abstract type AbstractTempProgram end

struct LinearProgram <:AbstractTempProgram
    dec::Float64
    LinearProgram(x::Real=1e-2) = new(x) 
end

function kthTemp(lin::LinearProgram,init_temp::Real,k::Integer)
    return init_temp - k*lin.dec
end

function updateTemp(lin::LinearProgram,curr_temp::Real,k::Integer)
    return curr_temp - lin.dec
end

struct ExponentialProgram <:AbstractTempProgram
    expon::Float64
    ExponentialProgram(x::Real=0.99) = new(x) 
end

function kthTemp(ex::ExponentialProgram,init_temp::Real,k::Integer)
    return init_temp*ex.expon^k
end

function updateTemp(ex::ExponentialProgram,curr_temp::Real,k::Integer)
    return curr_temp*ex.expon
end

struct LogarithmicProgram <:AbstractTempProgram
    base::Float64
    LogarithmicProgram(x::Real=2.0) = new(x) 
end

function kthTemp(loga::LogarithmicProgram,init_temp::Real,k::Integer)
    return log(loga.base,init_temp)/log(loga.base,k+1)
end

function updateTemp(loga::LogarithmicProgram,curr_temp::Real,k::Integer)
    return curr_temp*log(loga.base,k+1)/log(loga.base,k+2)
end

function singleSimulatedAnnealing(Q::PolygonalChain,
    minFunc::F,
    tolerance::Real=1e-2,
    alphamax::Real=pi/2,
    alphamin::Real=-pi/2,
    internal::Bool=false;
    temp_init::Float64=1.0,
    temp_f::Float64=1e-4,
    iter_per_temp::Integer=20,
    population::Integer=8,
    tempProgram::AbstractTempProgram=ExponentialProgram(),
    max_iter::Integer=1000,
    debug::Bool=false) where {F}
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    elseif minFunc == distToLine
        auxQ = makeLine(Q)
        minFunc = Q -> distToLine(Q,auxQ)
    end
    dist = Distributions.TruncatedNormal(0,pi/20,alphamin,alphamax)
    ang_idxs = zeros(Int16,max_iter,population)
    ang_vals = zeros(ftype(Q),max_iter,population)
    fun_vals = zeros(ftype(Q),max_iter+1,population)
    Qs = [perturbe(Q) for _ in 1:population]
    initQs = copy(Qs)
    fun_vals[1,:]  = [minFunc(R) for R in Qs]
    nq = length(Q)
    fval = minFunc(Q)
    temp = temp_init*abs(fval)
    temp_f = temp_f*abs(fval) 
    c = 1
    c2 = 1
    while c <= max_iter && fval > tolerance && temp > temp_f
        for i in 1:iter_per_temp
            for j in 1:population
                fval = fun_vals[c,j]
                #alpha = uniformAngle(alphamin,alphamax)
                alpha = rand(dist)
                (internal) && (alpha = alpha/2)
                ang_idx = generateAngleIndex(nq,internal)
                if debug
                    println("c = $c")
                    println("c2 = $(c2)")
                    println("i = $i")
                    println("Q = $Q")
                    println("ang_idx = $(ang_idx)")
                    println("ang_val  = $(alpha)")
                    println("# testing intersection")
                    println("\n")
                end
                inter_flag,newQ = checkSingleMutation(Qs[j],ang_idx,alpha,debug=debug)
                debug && println("inter  = $(inter_flag)")
                debug && println("newQ = $(newQ[j])")
                if !inter_flag
                    
                    fvalnew = minFunc(newQ)
                    r = log(rand())
                    p = (-fvalnew + fval)/temp
                    if debug
                        println("# no inter")
                        println("fval = $fval")
                        println("fvalnew = $(fvalnew)")
                        println("p = $p")
                        println("r = $r")
                    end
                    if r < p
                        debug && println("# accepted")
                        Qs[j] = newQ
                        fval            = fvalnew
                        ang_idxs[c,j]   = ang_idx
                        ang_vals[c,j]   = alpha
                        fun_vals[c+1,j] = fvalnew
                    else
                        fun_vals[c+1,j] = fval
                    end
                else
                    fun_vals[c+1,j] = fval
                end
            end
            c += 1
            debug && println("\n\n")
        end
        c2 += 1
        temp = updateTemp(tempProgram,temp,c2)
        println(temp)
    end
    c = c > max_iter ? max_iter : c
    return initQs,Qs,fun_vals[1:c+1,:],ang_vals[1:c,:],ang_idxs[1:c,:]
end

function multipleSimulatedAnnealing(Q::PolygonalChain,
    minFunc::F,
    tolerance::Real=1e-2,
    alphamax::Real=pi/2,
    alphamin::Real=-pi/2,
    internal::Bool=false;
    temp_init::Float64=1.0,
    temp_f::Float64=1e-4,
    iter_per_temp::Integer=20,
    tempProgram::AbstractTempProgram=ExponentialProgram(),
    population::Integer=8,
    max_iter::Integer=1000,
    mut_k::Integer=3,
    debug::Bool=false) where {F}
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    elseif minFunc == distToLine
        auxQ = makeLine(Q)
        minFunc = Q -> distToLine(Q,auxQ)
    end
    dist = Distributions.TruncatedNormal(0,pi/20,alphamin,alphamax)
    ang_idxs = zeros(Int16,mut_k*max_iter,population)
    ang_vals = zeros(ftype(Q),mut_k*max_iter,population)
    fun_vals = zeros(ftype(Q),max_iter+1,population)
    Qs = [perturbe(Q) for _ in 1:population]
    initQs = copy(Qs)
    fun_vals[1,:]  = [minFunc(R) for R in Qs]
    nq = length(Q)
    fval = minFunc(Q)
    temp = temp_init*abs(fval)
    temp_f = temp_f*abs(fval) 
    c = 1
    c2 = 1
    while  c <= max_iter && fval > tolerance && temp > temp_f
        for i in 1:iter_per_temp
            for j in 1:population
                fval = fun_vals[c,j]
                mut_length = rand(1:mut_k)
                alphas = [rand(dist) for _ in 1:mut_length]
                (internal) && (alphas = 0.5*alphas)
                mov_ang_idxs = internal ? Random.randperm(2*nq-3)[1:mut_length] : Random.randperm(nq-2)[1:mut_length]
                if debug
                    println("c = $c")
                    println("c2 = $(c2)")
                    println("i = $i")
                    println("Q = $(Qs[j])")
                    println("mov_ang_idxs = $(mov_ang_idxs)")
                    println("ang_val  = $(alphas)")
                    println("# testing intersection")
                    println("\n")
                end
                inter_flag,newQ = checkMultipleMutation(Qs[j],mov_ang_idxs,alphas,debug=debug)
                debug && println("inter  = $(inter_flag)")
                debug && println("newQ = $(newQ)")
                if !inter_flag
                    fvalnew = minFunc(newQ)
                    r = log(rand())
                    p = (-fvalnew + fval)/temp
                    if debug
                        println("# no inter")
                        println("fval = $fval")
                        println("fvalnew = $(fvalnew)")
                        println("temp = $(temp)")
                        println("p = $p")
                        println("r = $r")
                    end
                    if r < p
                        debug && println("# accepted")
                        Qs[j] = newQ
                        fval = fvalnew
                        li = mut_k*(c-1)+1
                        len = length(alphas)
                        ang_idxs[li:li+len-1,j] = mov_ang_idxs
                        ang_vals[li:li+len-1,j] = alphas
                        fun_vals[c+1,j] = fvalnew
                    else
                        fun_vals[c+1,j] = fval
                    end
                else
                    fun_vals[c+1,j] = fval
                end
            end
            c += 1
            debug && println("\n\n")
        end
        c2 += 1
        temp = updateTemp(tempProgram,temp,c2)
    end
    c = c > max_iter ? max_iter : c-1
    return initQs,Qs,fun_vals[1:c+1,:],ang_idxs[1:(c*mut_k),:],ang_vals[1:(c*mut_k),:]
end

# implementation taken from Wheeler et al "Algorithms for Optimization"
# current implementation avoids to make 

abstract type SelectionMethod end

struct TruncationSelection <: SelectionMethod
    k::Int64 
    TruncationSelection(k::Integer=3) = new(k)
end

function select(t::TruncationSelection,fvals::Array{<:Real,1})
    p = sortperm(fvals)
    return [p[rand(1:t.k)] for _ in fvals]
end

struct TournamentSelection <: SelectionMethod
    k::Int64 
    TournamentSelection(k::Integer=3) = new(k)
end

function select(t::TournamentSelection,fvals::Array{<:Real,1})
    n = length(fvals)
    vals = zeros(Int64,n) 
    for i in 1:n
        p = randperm(n)
        vals[i] = p[argmin(fvals[p[1:t.k]])]
    end
    return vals
end

struct RouletteWheelSelection <: SelectionMethod
    # initializer needs parameter to share interface with other SelectionMethod objects 
    RouletteWheelSelection(k::Integer=3) = new()
end

#=
function normalize(fvals::Array{<:Real,1})
    fvals = [Float64(maximum(fvals) - val) for val in fvals]
    return 0
end
=#

function select(t::RouletteWheelSelection,fvals::Array{<:Real,1})
    fvals = [Float64(maximum(fvals) - val) for val in fvals]
    if !isapprox(sum(fvals),0.0;atol=1e-15)
        pvec = LinearAlgebra.normalize(fvals,1)
    else
        # case where function values are all the same 
        println("fvals = $fvals")
        println("fvals are all the same")
        pvec = [1/length(fvals) for _ in fvals]
    end
    cat = Distributions.Categorical(pvec)
    return rand(cat,length(fvals))
end

function genetic(Q::PolygonalChain,
    minFunc::F,
    tolerance::Real=1e-2,
    alphamax::Real=pi/2,
    alphamin::Real=-pi/2,
    internal::Bool=false;
    selection::SelectionMethod=RouletteWheelSelection(),
    population::Integer=8,
    max_iter::Integer=1000,
    mut_k::Integer=3,
    debug::Bool=false) where {F}
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    elseif minFunc == distToLine
        auxQ = makeLine(Q)
        minFunc = Q -> distToLine(Q,auxQ)
    end
    ang_idxs = zeros(Int16,max_iter*mut_k,population)
    ang_vals = zeros(ftype(Q),max_iter*mut_k,population)
    fun_vals = zeros(ftype(Q),max_iter+1,population)
    parents = zeros(Int8,max_iter,population)
    dist = Distributions.TruncatedNormal(0,pi/20,alphamin,alphamax)
    nq = length(Q)
    Qs = [perturbe(Q) for _ in 1:population]
    initQs = copy(Qs)
    fun_vals[1,:]  = [minFunc(R) for R in Qs]
    c = 1
    while c <= max_iter && maximum(fun_vals[c,:]) > tolerance 
        # mutation
        for j in 1:population
            mut_length = rand(1:mut_k)
            alphas = [rand(dist) for _ in 1:mut_length]
            (internal) && (alphas = 0.5*alphas)
            mov_ang_idxs = internal ? Random.randperm(2*nq-3)[1:mut_length] : Random.randperm(nq-2)[1:mut_length]
            if debug
                println("c = $c")
                println("j = $j")
                println("Q = $(Qs[j])")
                println("ang_idxs = $(mov_ang_idxs)")
                println("ang_vals  = $(alphas)")
                println("# testing intersection")
                println("\n")
            end
            inter_flag,newQ = checkMultipleMutation(Qs[j],mov_ang_idxs,alphas,debug=false)
            if debug
                println("inter  = $(inter_flag)")
                println("newQ = $(newQ)")
            end
            if !inter_flag
                Qs[j] = newQ
                li = mut_k*(c-1)+1
                len = length(alphas)
                ang_idxs[li:li+len-1,j] = mov_ang_idxs
                ang_vals[li:li+len-1,j] = alphas
                fun_vals[c+1,j]         = minFunc(newQ)
            else
                fun_vals[c+1,j]         = fun_vals[c,j]
            end
        end
        # selection
        parents[c,:] = select(selection,fun_vals[c+1,:])
        if sum(abs.(diff(parents[c,:])))==0
            println(c)
            println("equal parents")
        end
        if debug
            for Q in Qs
                println(Q)
            end
        end
        # changing indexes
        Qs = Qs[parents[c,:]]
        initQs = initQs[parents[c,:]]
        ang_idxs[1:(c*mut_k),:] = ang_idxs[1:(c*mut_k),parents[c,:]]
        ang_vals[1:(c*mut_k),:] = ang_vals[1:(c*mut_k),parents[c,:]]
        fun_vals[1:(c+1),:]     = fun_vals[1:(c+1),parents[c,:]]
        c += 1
        debug && println("\n\n")
    end
    c = c > max_iter ? max_iter : c-1
    return initQs,Qs,fun_vals[1:c+1,:],ang_idxs[1:(c*mut_k),:],ang_vals[1:(c*mut_k),:],parents[1:c,:]
end



#=

function specialSimulatedAnnealing(Q::PolygonalChain,minFunc::Function,tolerance::Real=1e-2,phimax::Real=pi/2,phimin::Real=-pi/2; temp_init = 1,max_iter::Integer=1000,debug::Bool=false)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    ang_idxs = zeros(Int16,max_iter)
    ang_vals = zeros(typeof(Q[1].x),max_iter)
    fun_vals = zeros(typeof(Q[1].x),max_iter)
    nq = length(Q)
    debug && println("bien aca")
    temp = temp_init*minFunc(Q)
    delta_temp = (temp - 1e-6)/max_iter
    debug && println("intermedio")
    d = minFunc(Q)
    debug && println("mal aca")
    c = 1
    while d > tolerance && c <= max_iter
        
        debug && println(c)
        debug && println()
        debug && println(i)
        debug && println(Q)
        debug && println(linkLengths(Q))
        inter_flag = false
        newQ = copy(Q)
        for i in 1:(nq-2)
            newQ = moveBeforeDihedral(newQ,i)
            phi = rand()*(phimax-phimin) + phimin
            inter_flag = inter_flag || checkRotationIntersection(newQ,i,phi)
            newQ = dihedralRotate(newQ,i,phi)
        end
        debug && println(newQ)
        debug && println(linkLengths(newQ))
        dnew = minFunc(newQ)
        debug && println(c)
        debug && println(inter_flag)
        if !inter_flag && dnew < d
            Q = newQ
            d = dnew
            fun_vals[c] = dnew
        elseif !inter_flag && dnew >= d
            r = log(rand())
            p = (-dnew + d)/temp
            if r < p
                Q = newQ
                d = dnew
                fun_vals[c] = dnew
            else
                fun_vals[c] = d
            end
        else
            fun_vals[c] = d
        end
        c += 1
        temp -= delta_temp
    end
    return Q,ang_vals[1:(c-1)],ang_idxs[1:(c-1)],fun_vals[1:(c-1)]
end

=#

function weight(P::PolygonalChain,i::Integer,j::Integer,sigma::Real=0.75)
    li = distance(P[i+1],P[i])
    lj = distance(P[j+1],P[j])
    weight = 0
    for a in 0:1, b in 0:1
        weight += 1/(distance(P[i+a],P[j+b])^(2*sigma+1)) 
    end
    weight  = weight*li*lj/4
    return weight
end

function weight(P::PolygonalChain,i::Integer,j::Integer,li::Real,lj::Real,sigma::Real=0.75)
    weight = 0
    for a in 0:1, b in 0:1
        weight += 1/(distance(P[i+a],P[j+b])^(2*sigma+1)) 
    end
    weight  = weight*li*lj/4
    return weight
end

function Bmatrix(P::PolygonalChain,sigma::Real=0.75;debug::Bool=false)
    n = length(P)
    T = ftype(P)
    B = zeros(T,n+1,n+1)
    for i in 1:n
        # see `demaineEnergy` comment
        for j in 1:(i-2)
            li = distance(P[i+1],P[i])
            lj = distance(P[j+1],P[j])
            Ti = (P[i+1]-P[i])/li
            Tj = (P[j+1]-P[j])/lj
            tij = dot(Ti,Tj)
            wij = weight(P,i,j,li,lj,sigma)
            if debug
                println("i = $i")
                println("j = $j")
                println("li = $li")
                println("lj = $lj")
                println("Ti = $Ti")
                println("Tj = $Tj")
                println("tij = $tij")
                println("wij = $wij")
            end
            denom = li*lj
            for a in 0:1, b in 0:1
                sign = (-1)^(a+b)
                B[i+a,i+b] += sign*wij/li^2
                B[j+a,j+b] += sign*wij/lj^2
                B[i+a,j+b] -= sign*wij*tij/denom
                B[j+a,i+b] -= sign*wij*tij/denom
            end
        end
        for j in i+2:n
            li = distance(P[i+1],P[i])
            lj = distance(P[j+1],P[j])
            Ti = (P[i+1]-P[i])/li
            Tj = (P[j+1]-P[j])/lj
            tij = dot(Ti,Tj)
            wij = weight(P,i,j,li,lj,sigma)
            if debug
                println("i = $i")
                println("j = $j")
                println("li = $li")
                println("lj = $lj")
                println("Ti = $Ti")
                println("Tj = $Tj")
                println("tij = $tij")
                println("wij = $wij")
            end
            denom = li*lj
            for a in 0:1, b in 0:1
                sign = (-1)^(a+b)
                B[i+a,i+b] += sign*wij/li^2
                B[j+a,j+b] += sign*wij/lj^2
                B[i+a,j+b] -= sign*wij*tij/denom
                B[j+a,i+b] -= sign*wij*tij/denom
            end
        end
    end
    return B
end

function weight0(P::PolygonalChain,i::Integer,j::Integer,sigma::Real=0.75)
    li = distance(P[i+1],P[i])
    lj = distance(P[j+1],P[j])
    Ti = (P[i+1]-P[i])/li
    weight = 0
    for a in 0:1, b in 0:1
        weight += tangentPointKernel(P[i+a],P[j+a],Ti,2,4)/(distance(P[i+a],P[j+b])^(2*sigma+1)) 
    end
    weight  = weight*li*lj/4
    return weight
end

function weight0(P::PolygonalChain,i::Integer,j::Integer,li::Real,lj::Real,Ti::Point,sigma::Real=0.75)
    weight = 0
    for a in 0:1, b in 0:1
        weight += tangentPointKernel(P[i+a],P[j+a],Ti,2,4)/(distance(P[i+a],P[j+b])^(2*sigma+1)) 
    end
    weight = weight*li*lj/4
    return weight
end

function B0matrix(P::PolygonalChain,sigma::Real=0.75)
    n = length(P)
    T = ftype(P)
    B0 = zeros(T,n+1,n+1)
    for i in 1:n
        for j in 1:(i-2)
            w0ij = weight0(P,i,j,sigma)/4
            for a in 0:1, b in 0:1
                B0[i+a,i+b] += w0ij
                B0[j+a,j+b] += w0ij
                B0[i+a,j+b] -= w0ij
                B0[j+a,i+b] -= w0ij
            end
        end
        for j in (i+2):n
            w0ij = weight0(P,i,j,sigma)/4
            for a in 0:1, b in 0:1
                B0[i+a,i+b] += w0ij
                B0[j+a,j+b] += w0ij
                B0[i+a,j+b] -= w0ij
                B0[j+a,i+b] -= w0ij
            end
        end
    end
    return B0
end

function Amatrix(P::PolygonalChain,sigma::Real=0.75)
    n = length(P)
    T = ftype(P)
    B  = zeros(T,n+1,n+1)
    B0 = zeros(T,n+1,n+1)
    for i in 1:n
        li = distance(P[i+1],P[i])
        Ti = (P[i+1]-P[i])/li
        for j in 1:(i-2)
            lj = distance(P[j+1],P[j])
            Tj = (P[j+1]-P[j])/lj
            tij = dot(Ti,Tj)
            wij = weight(P,i,j,li,lj,sigma)
            denom = li*lj
            w0ij = weight0(P,i,j,li,lj,Ti,sigma)/4
            for a in 0:1, b in 0:1
                sign = (-1)^(a+b)
                B[i+a,i+b] += sign*wij/li^2
                B[j+a,j+b] += sign*wij/lj^2
                B[i+a,j+b] -= sign*wij*tij/denom
                B[j+a,i+b] -= sign*wij*tij/denom
                B0[i+a,i+b] += w0ij
                B0[j+a,j+b] += w0ij
                B0[i+a,j+b] -= w0ij
                B0[j+a,i+b] -= w0ij
            end
        end
        for j in i+2:n
            lj = distance(P[j+1],P[j])
            Tj = (P[j+1]-P[j])/lj
            tij = dot(Ti,Tj)
            wij = weight(P,i,j,li,lj,sigma)
            denom = li*lj
            w0ij = weight0(P,i,j,li,lj,Ti,sigma)/4
            for a in 0:1, b in 0:1
                sign = (-1)^(a+b)
                B[i+a,i+b] += sign*wij/li^2
                B[j+a,j+b] += sign*wij/lj^2
                B[i+a,j+b] -= sign*wij*tij/denom
                B[j+a,i+b] -= sign*wij*tij/denom
                B0[i+a,i+b] += w0ij
                B0[j+a,j+b] += w0ij
                B0[i+a,j+b] -= w0ij
                B0[j+a,i+b] -= w0ij
            end
        end
    end
    return B + B0
end

function Aline(P::PolygonalChain,sigma::Real=0.75)
    n = length(P)
    A = Amatrix(P,sigma)
    T = ftype(P)
    
    #=
    zer = zeros(T,n+1,n+1)
    fin = vcat(
        hcat(A,zer,zer),
        hcat(zer,A,zer),
        hcat(zer,zer,A)
    )
    =#
    id = [1 0 0 ; 0 1 0 ; 0 0 1]
    fin = zeros(T,3*n+3,3*n+3)
    for i in 1:n+1
        for j in 1:n+1
            fin[3*i-2:3*i,3*j-2:3*j] = A[i,j]*id
        end
    end
    return fin
end

function constraints(P::PolygonalChain,ls,bas,das,internal::Bool)
    n = length(P)
    T = ftype(P)
    if internal
        res = zeros(T,n)
        for i in 1:n
            res[i] = (distance(P[i+1],P[i])-ls[i])^2
        end
        return res
    else
        res = zeros(T,2*n-1)
        for i in 1:n
            res[i] = distance(P[i+1],P[i])-ls[i]
        end
        for i in 1:n-1
            res[n+i] = bangle(P[i],P[i+1],P[i+2])-bas[i]
        end
        return res
    end
end

#=

function bangle(x::Array{<:Real,1})
    u1x = x[4]-x[1]
    u1y = x[5]-x[2]
    u1z = x[6]-x[3]
    u2x = x[7]-x[4]
    u2y = x[8]-x[5]
    u2z = x[9]-x[6]
    d = u1x*u2x +  u1y*u2y + u1z*u2z
    n1 = distance([u1x,u1y,u1z,0,0,0]) 
    n2 = distance([u2x,u2y,u2z,0,0,0]) 
    return pi-acos(d/(n1*n2))
end
=#

function distance(x::Array{<:Real,1})
    return sqrt(pow2(x[1]-x[4]) + pow2(x[2]-x[5]) + pow2(x[3]-x[6]))
end

# wrappers to existing type seem to work better than implementations from zero

function bangle(x::Array{<:Real,1})
    return bangle(Point(x[1:3]),Point(x[4:6]),Point(x[7:9]))
end

function dihedral(x::Array{<:Real,1})
    return dihedral(Point(x[1:3]),Point(x[4:6]),Point(x[7:9]),Point(x[10:12]))
end

function bangleDev(a,b,c)
    h = 1e-3
    ah = Point(Hyperdual(a.x,h),Hyperdual(a.y,0),Hyperdual(a.z,0))
    return bangle(ah,b,c).a1/h
end


function deriv(f,p::Point,h::Real=1e-3)
    xp = Point(Hyperdual(p.x,h),Hyperdual(p.y,0),Hyperdual(p.z,0))
    yp = Point(Hyperdual(p.x,0),Hyperdual(p.y,h),Hyperdual(p.z,0))
    zp = Point(Hyperdual(p.x,0),Hyperdual(p.y,0),Hyperdual(p.z,h))
    return [f(xp).a1/h,f(yp).a1/h,f(zp).a1/h]
end


function hej(j::Integer,h::Real)
    T = typeof(h)
    if     j==1
        return Point(Hyperdual{T}(0.0,h,h,0.0),Hyperdual{T}(0.0,0.0,0.0,0.0),Hyperdual{T}(0.0,0.0,0.0,0.0))
    elseif j==2
        return Point(Hyperdual{T}(0.0,0.0,0.0,0.0),Hyperdual{T}(0.0,h,h,0.0),Hyperdual{T}(0.0,0.0,0.0,0.0))
    elseif j==3
        return Point(Hyperdual{T}(0.0,0.0,0.0,0.0),Hyperdual{T}(0.0,0.0,0.0,0.0),Hyperdual{T}(0.0,h,h,0.0))
    else
        error("index $j out of range")
    end
end

function deriv(f,P::PolygonalChain,h::Real=1e-3)
    n = length(P)
    T = ftype(P)
    res = zeros(T,3*n+3)
    hyperP = PolygonalChain([Point{Hyperdual}(p) for p in P.vertices])
    for i in 1:n+1
        #println("i = $i")
        for j in 1:3
            newP = copy(hyperP)
            #println("j = $j")
            newP[i] += hej(j,h)
            #println(newP)
            fval = f(newP)
            #println("f = $fval")
            res[3*(i-1)+j] = fval.a1/h
        end
    end
    return res
end

function diffLocalChange(f,P::PolygonalChain,h::Real=1e-3)
    n = length(P)
    T = ftype(P)
    res = zeros(T,3*n+3)
    for i in 1:n+1
        #println("i = $i")
        for j in 1:3
            #println("j = $j")
            newP = copy(P)
            newP[i] += hej(j,h)
            #println(newP)
            fval = f(newP)
            #println("f = $fval")
            res[3*(i-1)+j] = fval.a1/h
        end
    end
    return res
end

using DelimitedFiles

function constraints(P::PolygonalChain,ls,bas,das,internal::Bool)
    n = length(P)
    T = ftype(P)
    if internal
        res = zeros(T,n)
        for i in 1:n
            res[i] = distance(P[i+1],P[i])-ls[i]
        end
        return res
    else
        res = zeros(T,2*n-1)
        for i in 1:n
            res[i] = distance(P[i+1],P[i])-ls[i]
        end
        for i in 1:n-1
            res[n+i] = bangle(P[i],P[i+1],P[i+2])-bas[i]
        end
        return res
    end
end

function constraints(arr::Array{<:Real,1},ls,bas,das,internal::Bool)
    n = length(ls)
    T = typeof(arr[1])
    if internal
        res = zeros(T,n)
        for i in 1:n
            res[i] = distance(arr[3*i-2:3*i+3])-ls[i]
        end
        return res
    else
        res = zeros(T,2*n-1)
        for i in 1:n
            res[i] = distance(arr[3*i-2:3*i+3])-ls[i]
        end
        for i in 1:n-1
            res[n+i] = bangle(arr[3*i-2:3*i+6])-bas[i]
        end
        return res
    end
end

function hardConstraints(P::PolygonalChain,ls,bas,das)
    n = length(P)
    T = ftype(P)
    res = zeros(T,3*n-3)
    for i in 1:n
        res[i] = distance(P[i+1],P[i])-ls[i]
    end
    for i in 1:n-1
        res[n+i] = bangle(P[i],P[i+1],P[i+2])-bas[i]
    end
    d = (2*n-1)
    for i in 1:n-2
        res[d+i] = dihedral(P[i],P[i+1],P[i+2],P[i+3])-das[i]
    end
    return res
end

function hardConstraints(arr::Array{<:Real,1},ls,bas,das)
    n = length(ls)
    T = typeof(arr[1])
    res = zeros(T,3*n-3)
    for i in 1:n
        res[i] = distance(arr[3*i-2:3*i+3])-ls[i]
    end
    for i in 1:n-1
        res[n+i] = bangle(arr[3*i-2:3*i+6])-bas[i]
    end
    d = (2*n-1)
    for i in 1:n-2
        res[d+i] = dihedral(arr[3*i-2:3*i+9])-das[i]
    end
    return res
end

function helixifiedConstraints(P::PolygonalChain,ndens::Array{<:Integer,1},ls,bas,das;debug::Bool=false)
    nhelices = length(ndens)
    T = ftype(P)
    resn = 0
    for i in 1:nhelices 
        n = ndens[i]+1 # length of helix
        resn += 3*n-3
    end
    debug && println("resn=$resn")
    res = zeros(T,resn)
    c  = 1 # conunting vertex positions
    cl = 1 # counting lengths positions
    cb = 1 # counting bangles positions
    cd = 1 # counting dangles positions
    r = 0
    for i in 1:nhelices
        n  = ndens[i]+1 # length of helix
        d  = c+n
        dl = cl + n-1 
        db = cb + n-2
        dd = cd + n-3
        if debug
            println("c = $c, d = $d")
            println("cl = $cl,  dl = $dl")
            println("cb = $cb,  dd = $db")
            println("cd = $cl,  dd = $dd")
        end
        res[r+1:r+3*n-3] = hardConstraints(PolygonalChain(P[c:d]),ls[cl:dl],bas[cb:db],das[cd:dd])
        c  = d
        cl = dl+1
        cb = db+2
        cd = dd+3
        r += 3*n-3
    end
    debug && println("r=$r")
    return res
end

function helixifiedConstraints(arr::Array{<:Real,1},ndens::Array{<:Integer,1},ls,bas,das;debug::Bool=false)
    nhelices = length(ndens)
    T = typeof(arr[1])
    resn = 0
    for i in 1:nhelices 
        n = ndens[i]+1 # length of helix
        resn += 3*n-3
    end
    debug && println("resn=$resn")
    res = zeros(T,resn)
    c  = 1 # conunting vertex positions
    cl = 1 # counting lengths positions
    cb = 1 # counting bangles positions
    cd = 1 # counting dangles positions
    r = 0
    for i in 1:nhelices
        n  = ndens[i]+1 # length of helix
        d  = c+n
        dl = cl + n-1 
        db = cb + n-2
        dd = cd + n-3
        if debug
            println("c = $c, d = $d")
            println("cl = $cl,  dl = $dl")
            println("cb = $cb,  dd = $db")
            println("cd = $cl,  dd = $dd")
        end
        mat = hardConstraints(arr[3*c-2:3*d],ls[cl:dl],bas[cb:db],das[cd:dd])
        res[r+1:r+3*n-3] = mat
        c  = d
        cl = dl+1
        cb = db+2
        cd = dd+3
        r += 3*n-3
    end
    debug && println("r=$r")
    return res
end

## most performant version

function constraintsJacobian(P::PolygonalChain,internal::Bool)
    n = length(P)
    T = ftype(P)
    if internal
        C = zeros(T,n,3*n+3)
        # length restrictions
        for i in 1:n
            l = distance(P[i+1],P[i])
            C[i,3*i-2] = (P[i]-P[i+1]).x/l
            C[i,3*i-1] = (P[i]-P[i+1]).y/l
            C[i,3*i]   = (P[i]-P[i+1]).z/l
            C[i,3*i+1] = (P[i+1]-P[i]).x/l
            C[i,3*i+2] = (P[i+1]-P[i]).y/l
            C[i,3*i+3] = (P[i+1]-P[i]).z/l
        end
    else
        C = zeros(T,2*n-1,3*n+3)
        # length restrictions
        for i in 1:n
            l = distance(P[i+1],P[i])
            C[i,3*i-2] = (P[i]-P[i+1]).x/l
            C[i,3*i-1] = (P[i]-P[i+1]).y/l
            C[i,3*i]   = (P[i]-P[i+1]).z/l
            C[i,3*i+1] = (P[i+1]-P[i]).x/l
            C[i,3*i+2] = (P[i+1]-P[i]).y/l
            C[i,3*i+3] = (P[i+1]-P[i]).z/l
        end
        # internal angle restrictions
        for i in 1:n-1
            grad1 = deriv(p->bangle(p,P[i+1],P[i+2]),P[i])
            grad2 = deriv(p->bangle(P[i],p,P[i+2]),P[i+1])
            grad3 = deriv(p->bangle(P[i],P[i+1],p),P[i+2])
            C[n+i,3*i-2:3*i]   = grad1
            C[n+i,3*i+1:3*i+3] = grad2
            C[n+i,3*i+4:3*i+6] = grad3
        end
    end
    return C
end



function hardConstraintsJacobian(P::PolygonalChain)
    n = length(P)
    T = ftype(P)
    C = zeros(T,3*n-3,3*n+3)
    # length restrictions
    for i in 1:n
        l = distance(P[i+1],P[i])
        C[i,3*i-2] = (P[i]-P[i+1]).x/l
        C[i,3*i-1] = (P[i]-P[i+1]).y/l
        C[i,3*i]   = (P[i]-P[i+1]).z/l
        C[i,3*i+1] = (P[i+1]-P[i]).x/l
        C[i,3*i+2] = (P[i+1]-P[i]).y/l
        C[i,3*i+3] = (P[i+1]-P[i]).z/l
    end
    # internal angle restrictions
    for i in 1:n-1
        grad1 = deriv(p->bangle(p,P[i+1],P[i+2]),P[i])
        grad2 = deriv(p->bangle(P[i],p,P[i+2]),P[i+1])
        grad3 = deriv(p->bangle(P[i],P[i+1],p),P[i+2])
        C[n+i,3*i-2:3*i]   = grad1
        C[n+i,3*i+1:3*i+3] = grad2
        C[n+i,3*i+4:3*i+6] = grad3
    end
    # dihedral angle restrictions
    d = (2*n-1)
    for i in 1:n-2
        grad1 = deriv(p->dihedral(p,P[i+1],P[i+2],P[i+3]),P[i])
        grad2 = deriv(p->dihedral(P[i],p,P[i+2],P[i+3]),P[i+1])
        grad3 = deriv(p->dihedral(P[i],P[i+1],p,P[i+3]),P[i+2])
        grad4 = deriv(p->dihedral(P[i],P[i+1],P[i+2],p),P[i+3])
        C[d+i,3*i-2:3*i]   = grad1
        C[d+i,3*i+1:3*i+3] = grad2
        C[d+i,3*i+4:3*i+6] = grad3
        C[d+i,3*i+7:3*i+9] = grad4
    end
    return C
end

function helixifiedConstraintsJacobian(P::PolygonalChain,ndens::Array{<:Integer,1};debug::Bool=false)
    ntot = length(P)
    nhelices = length(ndens)
    T = ftype(P)
    resn = 0
    for i in 1:nhelices 
        n = ndens[i]+1 # length of helix
        resn += 3*n-3
    end
    debug && println("resn=$resn")
    res = zeros(T,resn,3*ntot+3)
    c  = 1 # conunting vertex positions
    r = 0
    for i in 1:nhelices
        n  = ndens[i]+1 # length of helix
        d  = c+n
        if debug
            println("c = $c, d = $d")
            println("res[r+1:r+3*n-3,3*c-2:3*d]=")
            println(res[r+1:r+3*n-3,3*c-2:3*d])
            println("hardConstraintsJacobian")
            println(hardConstraintsJacobian(PolygonalChain(P[c:d])))
        end
        mat = hardConstraintsJacobian(PolygonalChain(P[c:d]))
        res[r+1:r+3*n-3,3*c-2:3*d] = mat
        c  = d
        r += 3*n-3
    end
    debug && println("r=$r")
    return res
end

function sobolevGradient(P::PolygonalChain,internal::Bool,alpha::Real=2.0,beta::Real=4.5)
    sigma = ((beta-1)/alpha)-1
    n = length(P)
    T = ftype(P)
    k = internal ? n : 2*n-1
    ener = Q -> tangentEnergy(Q,alpha=alpha,beta=beta)
    ener_grad = deriv(ener,P)
    #println(ener_grad)
    #ener_grad = reshape(ener_grad,(3,n+1))
    #ener_grad = transpose(ener_grad)
    #ener_grad = reshape(ener_grad,length(ener_grad))
    #println(ener_grad)
    Al = Aline(P,sigma)
    C = constraintsJacobian(P,internal)
    mat = vcat(
        hcat(Al,transpose(C)),
        hcat(C,zeros(T,k,k))
    )
    ener_grad_ext = vcat(ener_grad,zeros(T,k))
    return IterativeSolvers.minres(mat,ener_grad_ext)
end

function movement(P::PolygonalChain,internal::Bool,
    ls::Array{<:Real,1},bas::Array{<:Real,1},das::Array{<:Real,1},
    alpha::Real=2.0,beta::Real=4.5;
    tau::Real=1e-2,debug::Bool=false)
    sigma = ((beta-1)/alpha)-1
    n = length(P)
    T = ftype(P)
    k = internal ? n : 2*n-1
    ener = Q -> tangentEnergy(Q,alpha=alpha,beta=beta)
    ener_grad = deriv(ener,P)
    #println(ener_grad)
    #ener_grad = reshape(ener_grad,(3,n+1))
    #ener_grad = transpose(ener_grad)
    #ener_grad = reshape(ener_grad,length(ener_grad))
    #println(ener_grad)
    Al = Aline(P,sigma)
    C = constraintsJacobian(P,internal)
    mat = vcat(
        hcat(Al,transpose(C)),
        hcat(C,zeros(T,k,k))
    )
    ener_grad_ext = vcat(ener_grad,zeros(T,k))
    cons_proy_dir = IterativeSolvers.minres(mat,ener_grad_ext)
    if debug
        open("../extra/Al.csv","w+") do io
            writedlm(io,Al,',')
        end
        open("../extra/mat.csv","w+") do io
            writedlm(io,mat,',')
        end
        open("../extra/ener.csv","w+") do io
            writedlm(io,ener_grad,',')
        end
        println("\n Al \n\n")
        display(Al)
        println()
        println("\n rank(Al) \n\n")
        display(rank(Al))
        println()
        println("\n ener_grad \n\n")
        display(ener_grad)
        println()
        #=
        sol = Al \ ener_grad
        println("\n sol \n\n")
        display(sol)
        println()
        println("\n Al*sol - ener_grad \n\n")
        display(Al*sol - ener_grad)
        =#
        println()
        println("\n mat \n\n")
        display(mat)
        println()
        println("\n ener_grad_ext \n\n")
        display(ener_grad_ext)
        println()
        println("\n cons_proy_dir \n\n")
        display(cons_proy_dir)
        println()
        println("\n mat*cons_proy_dir - ener_grad_ext \n\n")
        display(mat*cons_proy_dir - ener_grad_ext)
        println()
        println()
        println("\n cons_proy_dir \n\n")
        display(cons_proy_dir)
        println()
        
    end
    new_chain = toArray(P)-tau*cons_proy_dir[1:(3*n+3)]
    if debug
        println()
        println("\n new_chain \n\n")
        display(new_chain)
        println()
    end
    cons   = constraints(new_chain,ls,bas,das,internal)
    newSol = vcat(zeros(T,3*n+3),-1*cons)
    chain_proy_dir = IterativeSolvers.minres(mat,newSol)
    new_chain = new_chain + chain_proy_dir[1:(3*n+3)]
    c = 0
    cons = constraints(new_chain,ls,bas,das,internal)
    con_norm = LinearAlgebra.norm(cons,2)
    while con_norm > 1e-4 && c < 10
        if debug
            println("\n")
            println(con_norm)
            println("\n")
            display(cons)
            println("\n")
        end
        c += 1
        newSol[3*n+4:end] = -1*cons
        IterativeSolvers.minres!(chain_proy_dir,mat,newSol)
        new_chain = new_chain + chain_proy_dir[1:(3*n+3)]
        cons = constraints(new_chain,ls,bas,das,internal)
        con_norm = LinearAlgebra.norm(cons,2)
    end
    if debug
        println("\n c \n\n")
        println(c)
        println()
        println("\n cons \n\n")
        display(cons)
        println()
        println("\n chain_proy_dir \n\n")
        display(chain_proy_dir)
        println()
        println("\n new_chain_final \n\n")
        display(new_chain)
        println()
    end
    return new_chain
end




function sobolevGradientDescent(P::AbstractChain,
    iter::Integer
    ,alpha::Real = 2.0, beta::Real = 4.5;
    tau::Real=1e-2,
    debug::Bool=false)
    Q = copy(P)
    ls,bas,das = internalCoordinates(P)
    n = length(P)
    T = ftype(P)
    res = zeros(T,iter+1,3*n+3)
    res[1,:] = toArray(Q)
    for i in 1:iter
        if i % 10 == 0
            println("iter = $i")
        end
        newQ = PolygonalChain(res[i,:])
        res[i+1,:] = movement(newQ,true,ls,bas,das,alpha,beta,tau=tau,debug=debug)
    end
    return res
end

# movement for helixifiedChains
function movement(P::PolygonalChain,ndens::Array{<:Integer,1},
    ls::Array{<:Real,1},bas::Array{<:Real,1},das::Array{<:Real,1},
    alpha::Real=2.0,beta::Real=4.5;
    tau::Real=1e-2,debug::Bool=false)
    sigma = ((beta-1)/alpha)-1
    n = length(P)
    T = ftype(P)
    # counting restrictions
    k = 0
    for nden in ndens 
        nloc = nden+1 # length of helix
        k += 3*nloc-3
    end
    ener = Q -> tangentEnergy(Q,alpha=alpha,beta=beta)
    ener_grad = deriv(ener,P)
    Al = Aline(P,sigma)
    C = helixifiedConstraintsJacobian(P,ndens)
    mat = vcat(
        hcat(Al,transpose(C)),
        hcat(C,zeros(T,k,k))
    )
    ener_grad_ext = vcat(ener_grad,zeros(T,k))
    # it is very important to give the solver enough iterations to converge
    cons_proy_dir = IterativeSolvers.minres(mat,ener_grad_ext,maxiter=length(mat),reltol=1e-10)
    if debug
        open("../extra/Al.csv","w+") do io
            writedlm(io,Al,',')
        end
        open("../extra/mat.csv","w+") do io
            writedlm(io,mat,',')
        end
        open("../extra/ener.csv","w+") do io
            writedlm(io,ener_grad,',')
        end
        println("\n Al \n\n")
        display(Al)
        println()
        println("\n rank(Al) \n\n")
        display(rank(Al))
        println()
        println("\n ener_grad \n\n")
        display(ener_grad)
        println()
        #=
        sol = Al \ ener_grad
        println("\n sol \n\n")
        display(sol)
        println()
        println("\n Al*sol - ener_grad \n\n")
        display(Al*sol - ener_grad)
        =#
        println()
        println("\n mat \n\n")
        display(mat)
        println()
        println("\n ener_grad_ext \n\n")
        display(ener_grad_ext)
        println()
        println("\n cons_proy_dir \n\n")
        display(cons_proy_dir)
        println()
        println("\n mat*cons_proy_dir - ener_grad_ext \n\n")
        display(mat*cons_proy_dir - ener_grad_ext)
        println()
        println()
        println("\n cons_proy_dir \n\n")
        display(cons_proy_dir)
        println()
        
    end
    new_chain = toArray(P)-tau*cons_proy_dir[1:(3*n+3)]
    if debug
        println()
        println("\n new_chain \n\n")
        display(new_chain)
        println()
    end
    cons   = helixifiedConstraints(new_chain,ndens,ls,bas,das)
    newSol = vcat(zeros(T,3*n+3),-1*cons)
    chain_proy_dir = IterativeSolvers.minres(mat,newSol,maxiter=length(mat),reltol=1e-10)
    new_chain = new_chain + chain_proy_dir[1:(3*n+3)]
    c = 0
    cons = helixifiedConstraints(new_chain,ndens,ls,bas,das)
    con_norm = LinearAlgebra.norm(cons,2)
    while con_norm > 1e-4 && c < 10
        if debug
            println("\n")
            println(con_norm)
            println("\n")
            display(cons)
            println("\n")
        end
        c += 1
        newSol[3*n+4:end] = -1*cons
        IterativeSolvers.minres!(chain_proy_dir,mat,newSol,maxiter=length(mat),reltol=1e-10)
        new_chain = new_chain + chain_proy_dir[1:(3*n+3)]
        cons = helixifiedConstraints(new_chain,ndens,ls,bas,das)
        con_norm = LinearAlgebra.norm(cons,2)
    end
    if debug
        println("\n c \n\n")
        println(c)
        println()
        println("\n cons \n\n")
        display(cons)
        println()
        println("\n chain_proy_dir \n\n")
        display(chain_proy_dir)
        println()
        println("\n new_chain_final \n\n")
        display(new_chain)
        println()
    end
    return new_chain
end

function sobolevGradientDescent(P::AbstractChain,
    ndens::Array{<:Integer,1},
    iter::Integer,
    alpha::Real = 2.0,beta::Real = 4.5;
    tau::Real=1e-2,
    debug::Bool=false)
    Q = copy(P)
    ls,bas,das = internalCoordinates(P)
    n = length(P)
    T = ftype(P)
    res = zeros(T,iter+1,3*n+3)
    res[1,:] = toArray(Q)
    for i in 1:iter
        println("iter = $i")
        newQ = PolygonalChain(res[i,:])
        res[i+1,:] = movement(newQ,ndens,ls,bas,das,alpha,beta,tau=tau,debug=debug)
    end
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    test = false
    if test
        a = Point()
        b = Point()
        c = Point()
        println(jacobian(bangle,hcat(toArray(a),toArray(b),toArray(c)),1e-3;df=hyperdualPartialDerivative))
        println(bangleDev(a,b,c))
        auxFunc(p::Point) = bangle(p,b,c)
        println(deriv(auxFunc,a))
        aux2(P::PolygonalChain) = bangle(P[1],P[2],P[3])
        println(deriv(aux2,PolygonalChain([a,b,c])))
    end
    #P = PolygonalChain(6)
    #P = PolygonalChain([Point(BigFloat,p) for p in P.vertices])
    #P = curveToChain(trefoilCurve,30,0,deg2rad(315))
    #ls,bas,das = internalCoordinates(P)
    #=
    mov = movement(P,true,ls,bas,das,debug=true)
    
    Q = PolygonalChain(mov)
    ls,bas,das = internalCoordinates(Q)
    =#
    P = knittingNeedle(2.23)
    P = PolygonalChain([Point(Float64,p) for p in P.vertices])
    println(P)
    ndens = [5,2,2,2,5]
    Q = helixify(P,ndens;r=5e-3)
    ls,bas,das = internalCoordinates(Q)
    mov = movement(Q,ndens,ls,bas,das,debug=false)
    println(mov)
    println(norm(mov-toArray(Q)))
    #res = sobolevGradientDescent(P,1000,renom=false)
    #saveTrajectory()
    #=
    println("\n")
    display(res)
    println("\n")
    println(PolygonalChain(res[end,:]))
    =#
end

