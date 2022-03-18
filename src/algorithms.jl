include("intersections.jl")

using Statistics


function stair(n::Integer)::PolygonalChain
    #vertices = Array{Point,1}(undef,n+1)
    vertices = [e0 for i in 1:n+1]
    vertices[1] = e0
    vertices[2] = ex
    for i in 3:(n+1)
        if mod(i,2) == 1
            vertices[i] += vertices[i-1]+ ey
        else
            vertices[i] += vertices[i-1]+ ex
        end
    end
    return PolygonalChain(vertices)
end

function knittingNeedle(l::Real=2.0;ep::Real=1/6)::PolygonalChain
    p0 = ex - ep*ey
    p1 = rotate(ez,-pi/6,ex)
    p2 = e0
    p3 = ex
    p4 = ex + ez
    p5 = ep*ey
    p0 = p1 + l*unitVector(p0-p1)
    p5 = p4 + l*unitVector(p5-p4)
    points = [p0,p1,p2,p3,p4,p5]
    return PolygonalChain(points)
end

function fourKnot(l::Real=sqrt(2);ep::Real=0.1)::PolygonalChain
    v0 = l*ey + ep*ex
    v1 = ep*ez
    v2 = ex
    v3 = ex + ey
    v4 = ey + ep*ez
    v5 = e0
    v6 = l*ex + ep*ey + ep*ez
    vertices = [v0,v1,v2,v3,v4,v5,v6]
    #vertices = [v + 1e-8*Point() for v in vertices]
    return PolygonalChain(vertices)
end

function exactFourKnot(l::Real)
    ang1 = 0.012566
    ang3 = ang1
    ang2 = 0.131947
    ang4 = ang2
    v2 = ey
    v3 = e0
    v4 = ex
    # first move
    v1 = v2 + rotate(ex,ang1,-1*ey)
    v0 = v1 + l*rotate(-1*ey,ang2,unitVector(v1-v2))
    # second move
    v5 = v4 + rotate(ey,ang3,-1*ex)
    v6 = v5 + l*rotate(-1*ex,ang4,unitVector(v5-v4))
    return PolygonalChain([v0,v1,v2,v3,v4,v5,v6])
end

function flatten(P::AbstractChain)::PolygonalChain
    lengths, angles, dihedrals = lengthsAndAngles(P)
    newDiheds = [pi for i in dihedrals]
    return PolygonalChain(lengths,angles,newDiheds)
end

function distToFlat(Q::AbstractChain,P::AbstractChain=flatten(Q))::T
    if length(P) == length(Q)
        return overlapedRmsd(Q,P)
    else
        error("Chains must be of same length")
    end
end

function squaredMaxSpan(Q::AbstractChain)::T
    v = Q[1] - Q[end]
    return -dot(v,v)
end

function fourKnotSpan(Q::AbstractChain)::T
    v = Q[2] - Q[end-1]
    return -dot(v,v)
end

function demaineEnergy1(Q::AbstractChain)::T
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in i+2:(n+1)
            sum += 1/pow2(norm(Q[j]-Q[i]) + norm(Q[j]-Q[i+1]) - norm(Q[i]-Q[i+1]))
        end
    end
    return sum
end

function demaineEnergy2(Q::AbstractChain)::T
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

function tangentPointKernel(p::Point,q::Point,tang::Point,alpha::Real,beta::Real)::T
    dir = p-q
    return norm(cross(tang,dir))^alpha/norm(dir)^beta
end

function tangentPointDiscrete(Q::AbstractChain,i::Integer,j::Integer,
                            alpha::Real,beta::Real)::T
    sum = 0.0
    tang = unitVector(Q[i+1] - Q[i])
    sum += tangentPointKernel(Q[i],Q[j],tang,alpha,beta)
    sum += tangentPointKernel(Q[i],Q[j+1],tang,alpha,beta)
    sum += tangentPointKernel(Q[i+1],Q[j],tang,alpha,beta)
    sum += tangentPointKernel(Q[i+1],Q[j+1],tang,alpha,beta)
    return sum/4
end

function tangentEnergy(Q::AbstractChain;alpha::Real=3,beta::Real=6)::T
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in 1:(i-2)
            sum += tangentPointDiscrete(Q,i,j,alpha,beta)
        end
        for j in (i+2):n
            sum += tangentPointDiscrete(Q,i,j,alpha,beta)
        end
    end
    return sum
end

function tangentEnergyFrac(Q::AbstractChain)
    return tangentEnergy(Q,alpha=2,beta=4.5)
end


# we can use a single integer `k` to encode changes in both dihedral and internal angles
# if the integer is of size 1 <= k <= n-3, then it representes a change in dihedral angle
function generateAngleIndex(n::Integer,internal::Bool)::Integer
    return internal ? rand(1:(2*n-3)) : rand(1:(n-2))
end

function mutateChain(P::AbstractChain,ang_idx::Integer,alpha::Real;debug::Bool=false)
    n = length(P)
    possible = false
    if 1 <= ang_idx <= n-2
        newP = moveBeforeDihedral(P,ang_idx)
        possible = checkRotationIntersection(newP,ang_idx,alpha,false,debug=debug)
        if possible; dihedralRotateFast!(newP,ang_idx,alpha); end;
    elseif n-1 <= ang_idx <= 2*n-3
        ang_idx = ang_idx - n + 2 # adjusting to be smaller 
        newP = moveBeforeInternal(P,ang_idx)
        possible = checkRotationIntersection(newP,ang_idx,alpha,true,debug=debug)
        if possible; internalRotateFast!(newP,ang_idx,alpha); end;
    else
        error("index $(ang_idx) is not supported")
    end
    return possible,newP
end

function localRandomSearchStep(Q,newQ,inter_flag,c,minf_val,minf_newval,
                               minvals,dihed,diheds,phi,angles)
    if !inter_flag && minf_newval < minf_val
        diheds[c] = dihed
        angles[c] = phi
        minvals[c] = minf_newval
        return minf_newval, newQ
    else
        minvals[c] = minf_val
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

    diheds = zeros(Int16,max_iter)
    angles = zeros(T,max_iter)
    minvals = zeros(T,max_iter)
    nq = length(Q)
    minf_val = minFunc(Q)
    c = 1
    while minf_val > tolerance && c <= max_iter
        phi = rand()*(phimax-phimin) + phimin
        dihed = rand(1:(nq-2))
        newQ = moveBeforeDihedral(Q,dihed)
        inter_flag = checkRotationIntersection(newQ,dihed,phi)
        newQ = dihedralRotate(newQ,dihed,phi)
        minf_newval = minFunc(newQ)
        minf_val,Q = advanceFunc!(Q,newQ,inter_flag,c,minf_val,minf_newval,
                                  minvals,dihed,diheds,phi,angles)
        c += 1
    end
    return Q,angles[1:(c-1)],diheds[1:(c-1)],minvals[1:(c-1)]
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
    phimax::Real=pi/2,
    phimin::Real=-pi/2,
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
    diheds  = zeros(Int16,max_iter)
    angles  = zeros(T,max_iter)
    minvals = zeros(T,max_iter)
    nq = length(Q)
    d = minFunc(Q)
    c = 1
    while d > tolerance && c <= max_iter
        phi = rand()*(phimax-phimin) + phimin
        (internal) && (phi = phi/2)
        dihed = generateAngleIndex(nq,internal)
        inter_flag , newQ = mutateChain(Q,dihed,phi)
        if !inter_flag
            dnew = minFunc(newQ)
            if dnew < d
                Q = newQ
                d = dnew
                diheds[c] = dihed
                angles[c] = phi
                minvals[c] = dnew
            end
        else
            minvals[c] = d
        end
        c += 1
    end
    c = c > max_iter ? max_iter : c
    return Q,angles[1:c],diheds[1:c],minvals[1:c]
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
    phimax::Real=pi/2,
    phimin::Real=-pi/2,
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
    diheds = zeros(Int16,max_iter)
    angles = zeros(T,max_iter)
    minvals = zeros(T,max_iter)
    nq = length(Q)
    d = minFunc(Q)
    temp = temp_init*d
    temp_f = temp_f*d 
    c = 1
    c2 = 1
    while d > tolerance && c <= max_iter && temp > temp_f
        for i in 1:iter_per_temp
            phi = rand()*(phimax-phimin) + phimin
            (internal) && (phi = phi/2)
            dihed = generateAngleIndex(nq,internal)
            inter_flag , newQ = mutateChain(Q,dihed,phi)
            debug && println(c)
            debug && println()
            debug && println(i)
            debug && println(Q)
            debug && println(linkLengths(Q))
            debug && println(newQ)
            debug && println(linkLengths(newQ))
            debug && println(c)
            debug && println(inter_flag)
            if !inter_flag
                dnew = minFunc(newQ)
                r = log(rand())
                p = (-dnew + d)/temp
                if r < p
                    Q = newQ
                    d = dnew
                    diheds[c] = dihed
                    angles[c] = phi
                    minvals[c] = dnew

                end
            else
                minvals[c] = d
            end
            c += 1
        end
        c2 += 1
        temp = tempUpdate(temp,c2)
    end
    c = c > max_iter ? max_iter : c
    return Q,angles[1:c],diheds[1:c],minvals[1:c]
end

#=
function simulatedAnnealing(Q::PolygonalChain,minFunc::Function,tolerance::Real=1e-2,phimax::Real=pi/2,phimin::Real=-pi/2; temp_init = 1,max_iter::Integer=1000,debug::Bool=false)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    diheds = zeros(Int16,max_iter)
    angles = zeros(T,max_iter)
    minvals = zeros(T,max_iter)
    nq = length(Q)
    debug && println("bien aca")
    temp = temp_init*minFunc(Q)
    delta_temp = (temp - 1e-6)/max_iter
    debug && println("intermedio")
    d = minFunc(Q)
    debug && println("mal aca")
    c = 1
    while d > tolerance && c <= max_iter
        phi = rand()*(phimax-phimin) + phimin
        dihed = rand(1:(nq-2))
        debug && println(c)
        debug && println()
        debug && println(i)
        debug && println(Q)
        debug && println(linkLengths(Q))
        newQ = moveBeforeDihedral(Q,dihed)
        debug && println(newQ)
        debug && println(linkLengths(newQ))
        inter_flag = checkRotationIntersection(newQ,dihed,phi)
        newQ = dihedralRotate(newQ,dihed,phi)
        dnew = minFunc(newQ)
        debug && println(c)
        debug && println(inter_flag)
        if !inter_flag && dnew < d
            Q = newQ
            d = dnew
            diheds[c] = dihed
            angles[c] = phi
            minvals[c] = dnew
        elseif !inter_flag && dnew >= d
            r = log(rand())
            p = (-dnew + d)/temp
            if r < p
                Q = newQ
                d = dnew
                diheds[c] = dihed
                angles[c] = phi
                minvals[c] = dnew
            else
                minvals[c] = d
            end
        else
            minvals[c] = d
        end
        c += 1
        temp -= delta_temp
    end
    return Q,angles[1:(c-1)],diheds[1:(c-1)],minvals[1:(c-1)]
end

function specialSimulatedAnnealing(Q::PolygonalChain,minFunc::Function,tolerance::Real=1e-2,phimax::Real=pi/2,phimin::Real=-pi/2; temp_init = 1,max_iter::Integer=1000,debug::Bool=false)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    diheds = zeros(Int16,max_iter)
    angles = zeros(T,max_iter)
    minvals = zeros(T,max_iter)
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
            minvals[c] = dnew
        elseif !inter_flag && dnew >= d
            r = log(rand())
            p = (-dnew + d)/temp
            if r < p
                Q = newQ
                d = dnew
                minvals[c] = dnew
            else
                minvals[c] = d
            end
        else
            minvals[c] = d
        end
        c += 1
        temp -= delta_temp
    end
    return Q,angles[1:(c-1)],diheds[1:(c-1)],minvals[1:(c-1)]
end


function linearSimulatedAnnealing(Q::PolygonalChain,minFunc::Function,tolerance::Real=1e-2,phimax::Real=pi/2,phimin::Real=-pi/2; temp_init = 1,max_iter::Integer=1000,debug::Bool=false)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    diheds = zeros(Int16,max_iter)
    angles = zeros(T,max_iter)
    minvals = zeros(T,max_iter)
    nq = length(Q)
    debug && println("bien aca")
    temp_range = range(temp_init*minFunc(Q),stop=1e-6,length=max_iter)
    debug && println("intermedio")
    d = minFunc(Q)
    debug && println("mal aca")
    c = 1
    while d > tolerance && c <= max_iter
        phi = rand()*(phimax-phimin) + phimin
        dihed = rand(1:(nq-2))
        debug && println(c)
        debug && println()
        debug && println(i)
        debug && println(Q)
        debug && println(linkLengths(Q))
        newQ = moveBeforeDihedral(Q,dihed)
        debug && println(newQ)
        debug && println(linkLengths(newQ))
        inter_flag = checkRotationIntersection(newQ,dihed,phi)
        newQ = dihedralRotate(newQ,dihed,phi)
        dnew = minFunc(newQ)
        debug && println(c)
        debug && println(inter_flag)
        if !inter_flag && dnew < d
            Q = newQ
            d = dnew
            diheds[c] = dihed
            angles[c] = phi
            minvals[c] = dnew
        elseif !inter_flag && dnew >= d
            r = log(rand())
            p = (-dnew + d)/temp_range[c]
            if r < p
                Q = newQ
                d = dnew
                diheds[c] = dihed
                angles[c] = phi
                minvals[c] = dnew
            else
                minvals[c] = d
            end
        else
            minvals[c] = d
        end
        c += 1
    end
    return Q,angles[1:(c-1)],diheds[1:(c-1)],minvals
end
=#

if abspath(PROGRAM_FILE) == @__FILE__
end