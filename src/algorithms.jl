include("intersections.jl")

using IterativeSolvers # julia's `\` operator is not good for solving this systems
using Statistics


function distToFlat(Q::AbstractChain,P::AbstractChain=flatten(Q))
    if length(P) == length(Q)
        return overlapedRmsd(Q,P)
    else
        error("Chains must be of same length")
    end
end

function squaredMaxSpan(Q::AbstractChain)
    v = Q[1] - Q[end]
    return -dot(v,v)
end

function fourKnotSpan(Q::AbstractChain)
    v = Q[2] - Q[end-1]
    return -dot(v,v)
end

function demaineEnergy1(Q::AbstractChain)
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in i+2:(n+1)
            sum += 1/pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
        end
    end
    return sum
end

function demaineEnergy2(Q::AbstractChain)
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in 1:(i-1)
            sum += 1/pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
        end
        for j in (i+2):(n+1)
            sum += pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
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

function tangentEnergy2(Q::AbstractChain;alpha::Real=3,beta::Real=6)
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


# we can use a single integer `k` to encode changes in both dihedral and internal ang_vals
# if the integer is of size 1 <= k <= n-3, then it representes a change in dihedral angle
function generateAngleIndex(n::Integer,internal::Bool)::Integer
    return internal ? rand(1:(2*n-3)) : rand(1:(n-2))
end

function mutateChain(P::AbstractChain,ang_idx::Integer,alpha::Real;debug::Bool=false)
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
    ang_vals = zeros(typeof(Q[1].x),max_iter)
    fun_vals = zeros(typeof(Q[1].x),max_iter)
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
    minFunc::Function,
    tolerance::Real=1e-2,
    alphamax::Real=pi/2,
    alphamin::Real=-pi/2,
    internal::Bool=false;
    temp_init::Float64=1.0,
    temp_f::Float64=1e-4,
    iter_per_temp::Integer=20,
    max_iter::Integer=1000,
    debug::Bool=false)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    ang_idxs  = zeros(Int16,max_iter)
    ang_vals  = zeros(typeof(Q[1].x),max_iter)
    fun_vals =  zeros(typeof(Q[1].x),max_iter)
    nq = length(Q)
    d = minFunc(Q)
    c = 1
    while d > tolerance && c <= max_iter
        alpha = rand()*(alphamax-alphamin) + alphamin
        (internal) && (alpha = alpha/2)
        ang_idx = generateAngleIndex(nq,internal)
        inter_flag , newQ = mutateChain(Q,ang_idx,alpha)
        if !inter_flag
            dnew = minFunc(newQ)
            if dnew < d
                Q = newQ
                d = dnew
                ang_idxs[c] = ang_idx
                ang_vals[c] = alpha
                fun_vals[c] = dnew
            end
        else
            fun_vals[c] = d
        end
        c += 1
    end
    c = c > max_iter ? max_iter : c
    return Q,ang_vals[1:c],ang_idxs[1:c],fun_vals[1:c]
end

function linearTemperature(temp::Real,k::Integer;b::Float64=1e-2)
    return temp - b
end

function exponentialTemperature(temp::Real,k::Integer;a::Float64=0.99)
    return a*temp
end

function simulatedAnnealing(Q::PolygonalChain,
    minFunc::Function,
    tolerance::Real=1e-2,
    alphamax::Real=pi/2,
    alphamin::Real=-pi/2,
    internal::Bool=false;
    temp_init::Float64=1.0,
    temp_f::Float64=1e-4,
    iter_per_temp::Integer=20,
    tempUpdate=exponentialTemperature,
    max_iter::Integer=1000,
    debug::Bool=false)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    ang_idxs = zeros(Int16,max_iter)
    ang_vals = zeros(typeof(Q[1].x),max_iter)
    fun_vals = zeros(typeof(Q[1].x),max_iter)
    nq = length(Q)
    d = minFunc(Q)
    temp = temp_init*abs(d)
    temp_f = temp_f*abs(d) 
    c = 1
    c2 = 1
    while d > tolerance && c <= max_iter && temp > temp_f
        for i in 1:iter_per_temp
            alpha = rand()*(alphamax-alphamin) + alphamin
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
            inter_flag,newQ = mutateChain(Q,ang_idx,alpha,debug=debug)
            debug && println("inter  = $(inter_flag)")
            debug && println("newQ = $(newQ)")
            if !inter_flag
                
                dnew = minFunc(newQ)
                r = log(rand())
                p = (-dnew + d)/temp
                if debug
                    println("# no inter")
                    println("d = $d")
                    println("dnew = $(dnew)")
                    println("p = $p")
                    println("r = $r")
                end
                if r < p
                    debug && println("# accepted")
                    Q = newQ
                    d = dnew
                    ang_idxs[c] = ang_idx
                    ang_vals[c] = alpha
                    fun_vals[c] = dnew

                end
            else
                fun_vals[c] = d
            end
            c += 1
            debug && println("\n\n")
        end
        c2 += 1
        temp = tempUpdate(temp,c2)
    end
    c = c > max_iter ? max_iter : c
    return Q,ang_vals[1:c],ang_idxs[1:c],fun_vals[1:c]
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
    T = typeof(P[1].x)
    B = zeros(T,n+1,n+1)
    for i in 1:n
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
    T = typeof(P[1].x)
    B0 = zeros(T,n+1,n+1)
    for i in 1:n
        for j in i+2:n
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
    T = typeof(P[1].x)
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
    T = typeof(P[1].x)
    
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
    T = typeof(P[1].x)
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

function distance(x)
    ax = x[1]
    ay = x[2]
    az = x[3]
    bx = x[4]
    by = x[5]
    bz = x[6]
    return sqrt((ax-bx)*(ax-bx) + (ay-by)*(ay-by) + (az-bz)*(az-bz))
end

function bangle(x)
    ax = x[1]
    ay = x[2]
    az = x[3]
    bx = x[4]
    by = x[5]
    bz = x[6]
    cx = x[7]
    cy = x[8]
    cz = x[9]
    u1x = bx-ax
    u1y = by-ay
    u1z = bz-az
    u2x = cx-bx
    u2y = cy-by
    u2z = cz-bz
    d = u1x*u2x +  u1y*u2y + u1z*u2z
    n1 = distance([u1x,u1y,u1z,0,0,0]) 
    n2 = distance([u2x,u2y,u2z,0,0,0]) 
    return pi-acos(d/(n1*n2))
end

function bangleDev(a,b,c)
    h = 1e-3
    ah = Point(Hyperdual(a.x,h),Hyperdual(a.y,0),Hyperdual(a.z,0))
    bh = Point(Hyperdual(b.x,0),Hyperdual(b.y,0),Hyperdual(b.z,0))
    ch = Point(Hyperdual(c.x,0),Hyperdual(c.y,0),Hyperdual(c.z,0))
    return bangle(ah,b,c).a1/h
end


function diff(f,p::Point,h::Real=1e-3)
    res = zeros(typeof(p.x),3)
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

function diff(f,P::PolygonalChain,h::Real=1e-3)
    n = length(P)
    T = typeof(P[1].x)
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
    T = typeof(P[1].x)
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
    T = typeof(P[1].x)
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


function constraintsJacobian(P::PolygonalChain,internal::Bool)
    n = length(P)
    T = typeof(P[1].x)
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
        arr = to2DArray(P)
        for i in 1:n-1
            lim1 = 3*i-2
            lim2 = 3*i+6
            x = transpose(arr)[lim1:lim2]
            dev = jacobian(bangle,x,1e-3;df=hyperdualPartialDerivative)
            C[n+i,lim1:lim2] = dev
        end
    end
    return C
end

function constraintsJacobian3(P::PolygonalChain,internal::Bool)
    n = length(P)
    T = typeof(P[1].x)
    if internal
        C = zeros(T,n,3*n+3)
        # length restrictions
        for i in 1:n
            grad1 = diff(p->distance(P[i+1],p),P[i])
            grad2 = diff(p->distance(p,P[i]),P[i+1])
            C[i,3*i-2:3*i]   = grad1
            C[i,3*i+1:3*i+3] = grad2
        end
    else
        C = zeros(T,2*n-1,3*n+3)
        # length restrictions
        for i in 1:n
            grad1 = diff(p->distance(P[i+1],p),P[i])
            grad2 = diff(p->distance(p,P[i]),P[i+1])
            C[i,3*i-2:3*i]   = grad1
            C[i,3*i+1:3*i+3] = grad2
        end
        # internal angle restrictions
        arr = to2DArray(P)
        for i in 1:n-1
            grad1 = diff(p->bangle(p,P[i+1],P[i+2]),P[i])
            grad2 = diff(p->bangle(P[i],p,P[i+2]),P[i+1])
            grad3 = diff(p->bangle(P[i],P[i+1],p),P[i+2])
            C[n+i,3*i-2:3*i]   = grad1
            C[n+i,3*i+1:3*i+3] = grad2
            C[n+i,3*i+4:3*i+6] = grad3
        end
    end
    return C
end

function constraintsJacobian4(P::PolygonalChain,internal::Bool)
    n = length(P)
    T = typeof(P[1].x)
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
        arr = to2DArray(P)
        for i in 1:n-1
            grad1 = diff(p->bangle(p,P[i+1],P[i+2]),P[i])
            grad2 = diff(p->bangle(P[i],p,P[i+2]),P[i+1])
            grad3 = diff(p->bangle(P[i],P[i+1],p),P[i+2])
            C[n+i,3*i-2:3*i]   = grad1
            C[n+i,3*i+1:3*i+3] = grad2
            C[n+i,3*i+4:3*i+6] = grad3
        end
    end
    return C
end


function constraintsJacobian2(P::PolygonalChain,internal::Bool)
    n = length(P)
    T = typeof(P[1].x)
    if internal
        C = zeros(T,n,3*n+3)
        # length restrictions
        for i in 1:n
            #println(i)
            grad = diff(Q->distance(Q[i+1],Q[i]),P)
            #println(grad)
            C[i,:] = grad
        end
    else
        C = zeros(T,2*n-1,3*n+3)
        # length restrictions
        for i in 1:n
            grad = diff(Q->distance(Q[i+1],Q[i]),P)
            #println(grad)
            C[i,:] = grad
        end
        # internal angle restrictions
        arr = to2DArray(P)
        for i in 1:n-1
            grad = diff(Q->bangle(Q[i],Q[i+1],Q[i+2]),P)
            C[n+i,:] = grad
        end
    end
    return C
end


function movement(P::PolygonalChain,internal::Bool,
    ls::Array{<:Real,1},bas::Array{<:Real,1},das::Array{<:Real,1},
    alpha::Real=2.0,beta::Real=4.5;
    tau::Real=1e-2,debug::Bool=false)
    sigma = ((beta-1)/alpha)-1
    n = length(P)
    T = typeof(P[1].x)
    k = internal ? n : 2*n-1
    ener = Q -> tangentEnergy(Q,alpha=alpha,beta=beta)
    ener_grad = diff(ener,P)
    #println(ener_grad)
    #ener_grad = reshape(ener_grad,(3,n+1))
    #ener_grad = transpose(ener_grad)
    #ener_grad = reshape(ener_grad,length(ener_grad))
    #println(ener_grad)
    Al = Aline(P,sigma)
    C = constraintsJacobian4(P,internal)
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

function sobolevGradientDescent(P::AbstractChain,iter::Integer,alpha::Real = 2.0, beta::Real = 4.5;tau::Real=1e-2,debug::Bool=false,renom::Bool=false)
    Q = copy(P)
    ls,bas,das = internalCoordinates(P)
    n = length(P)
    T = typeof(P[1].x)
    res = zeros(T,iter+1,3*n+3)
    res[1,:] = toArray(Q)
    for i in 1:iter
        if i % 100 == 0
            println("iter = $i")
        end
        newQ = PolygonalChain(res[i,:])
        res[i+1,:] = movement(newQ,true,ls,bas,das,alpha,beta,tau=tau,debug=debug)
    end
    return res
end


function assertEqualJac(n::Integer)
    eq_flag = false
    for i in 1:n
        C1 = constraintsJacobian(P,false)
        C2 = constraintsJacobian2(P,false)
        eq_flag = eq_flag || isapprox(C1-C2,1e-15)
        if eq_flag
            println()
            display(C1)
            println()
            display(C2)
            println()
            display(C1-C2)
            println()
        end
    end
end

function benchConsJac1()
    P = PolygonalChain(20)
    constraintsJacobian(P,false)
end

function benchConsJac2()
    P = PolygonalChain(20)
    constraintsJacobian2(P,false)
end

function benchConsJac3()
    P = PolygonalChain(20)
    constraintsJacobian3(P,false)
end

function benchConsJac4()
    P = PolygonalChain(20)
    constraintsJacobian4(P,false)
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
        println(diff(auxFunc,a))
        aux2(P::PolygonalChain) = bangle(P[1],P[2],P[3])
        println(diff(aux2,PolygonalChain([a,b,c])))
    end
    bench = false
    if bench
        using BenchmarkTools
        println()
        println()
        display(@benchmark benchConsJac1())
        println()
        println()
        display(@benchmark benchConsJac2())
        println()
        println()
        display(@benchmark benchConsJac3())
        println()
        println()
        display(@benchmark benchConsJac4())
        println()
        println()
    end
    #P = PolygonalChain(6)
    #P = PolygonalChain([Point(BigFloat,p) for p in P.vertices])
    P = parametricCurveChain(treefoil,30,0,deg2rad(315))
    
    ls,bas,das = internalCoordinates(P)
    #=
    mov = movement(P,true,ls,bas,das,debug=true)
    
    Q = PolygonalChain(mov)
    ls,bas,das = internalCoordinates(Q)
    
    P = knittingNeedle(2.23)
    P = PolygonalChain([Point(Float64,p) for p in P.vertices])
    println(P)
    =#
    res = sobolevGradientDescent(P,1000,renom=false)
    saveTrajectory()
    #=
    println("\n")
    display(res)
    println("\n")
    println(PolygonalChain(res[end,:]))
    =#
end