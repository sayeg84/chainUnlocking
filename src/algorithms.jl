include("intersections.jl")

using Statistics

function moveBeforeDihedral(P::PolygonalNew,i::Integer)::PolygonalNew
    newP = copy(P)
    n = length(P)
    if 1 <= i <= n-2
        u = unitVector(newP[i+2]-newP[i+1])
        mat = planeRotationXY(u)
        pivot = newP[i+2]
        for j in 1:n+1
            newP[j] = mat*(newP[j]-pivot)
        end
        return newP
    else
        error("i is not in range of dihedrals")
    end
end


function tryDihedralRotation!(P::PolygonalNew,i::Integer,theta::Real)
    n = length(P)
    Q = moveBeforeDihedral(P,i)
    if 1 <= i < n-1
        # TODO change to cells and make brute force
        b = checkRotationIntersection(Q[1:i],Q[i+1:end],theta)
        if b
            dihedralRotate!(Q,i,theta)
            return Q
        else
            println("Rotation by angle `theta' makes intersection")
            return Q
        end
    else
        error("``i is not between 1 and ``n")
    end
end

function stair(n::Integer)::PolygonalNew
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
    return PolygonalNew(vertices)
end

function knittingNeedle(l::Real=2.0;ep::Real=1/6)::PolygonalNew
    p0 = ex - ep*ey
    p1 = rotate(ez,-pi/6,ex)
    p2 = e0
    p3 = ex
    p4 = ex + ez
    p5 = ep*ey
    p0 = p1 + l*unitVector(p0-p1)
    p5 = p4 + l*unitVector(p5-p4)
    points = [p0,p1,p2,p3,p4,p5]
    return PolygonalNew(points)
end

function fourKnot(l::Real=sqrt(2);ep::Real=0.1)::PolygonalNew
    v0 = l*ey + ep*ex
    v1 = ep*ez
    v2 = ex
    v3 = ex + ey
    v4 = ey + ep*ez
    v5 = e0
    v6 = l*ex + ep*ey + ep*ez
    vertices = [v0,v1,v2,v3,v4,v5,v6]
    vertices = [v + 1e-8*Point() for v in vertices]
    return PolygonalNew(vertices)
end

function flatten(P::AbstractChain)::PolygonalNew
    lengths, angles, dihedrals = lengthsAndAngles(P)
    newDiheds = [pi for i in dihedrals]
    return PolygonalNew(lengths,angles,newDiheds)
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

function demaineEnergy1(Q::AbstractChain)::T
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in i+2:(n+1)
            sum += 1/pow2(norm(Q[j]-Q[i])+ norm(Q[j]-Q[i+1]) -norm(Q[i]-Q[i+1]))
        end
    end
    return sum
end

function demaineEnergy2(Q::AbstractChain)::T
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in 1:(i-1)
            sum += 1/pow2(norm(Q[j]-Q[i])+ norm(Q[j]-Q[i+1]) -norm(Q[i]-Q[i+1]))
        end
        for j in (i+2):(n+1)
            sum += pow2(norm(Q[j]-Q[i])+ norm(Q[j]-Q[i+1]) -norm(Q[i]-Q[i+1]))
        end
    end
    return sum
end

function tangentPointKernel(tang::Point,p::Point,q::Point;alpha::Real=2,beta=1)::T
    dir = p-q
    return norm(cross(tang,dir))^alpha/norm(dir)^beta
end

function tangentPointDiscrete(Q::AbstractChain,i::Integer,j::Integer)::T
    sum = 0.0
    tang = unitVector(Q[i+1] - Q[i])
    sum += tangentPointKernel(Q[i],Q[j],tang)
    sum += tangentPointKernel(Q[i],Q[j+1],tang)
    sum += tangentPointKernel(Q[i+1],Q[j],tang)
    sum += tangentPointKernel(Q[i+1],Q[j+1],tang)
    return sum/4
end

function tangentEnergy(Q::AbstractChain)::T
    n = length(Q)
    sum = 0
    for i in 1:n
        for j in 1:(i-2)
            sum += tangentPointDiscrete(Q,i,j)
        end
        for j in (i+2):n
            sum += tangentPointDiscrete(Q,i,j)
        end
    end
    return sum
end


function localSearchRandom(Q::PolygonalNew,minFunc::Function,tolerance::Real=1e-2,thetamax::Real=pi/2,thetamin::Real=-pi/2;max_iter::Integer=1000)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    diheds = zeros(Int8,max_iter)
    angles = zeros(T,max_iter)
    nq = length(Q)
    #println("bien aca")
    #println("intermedio")
    d = minFunc(Q)
    #println("mal aca")
    c = 1
    while d > tolerance && c <= max_iter
        theta = rand()*(thetamax-thetamin) + thetamin
        dihed = rand(1:(nq-2))
        #println(c)
        #println()
        #println(i)
        #println(Q)
        #println(linkLengths(Q))
        newQ = moveBeforeDihedral(Q,dihed)
        #println(newQ)
        #println(linkLengths(newQ))
        inter_flag = checkRotationIntersection(newQ,dihed,theta)
        newQ = dihedralRotate(newQ,dihed,theta)
        dnew = minFunc(newQ)
        #println(c)
        #println(inter_flag)
        #println()
        if !inter_flag #&& dnew < d
            Q = newQ
            d = dnew
            diheds[c] = dihed
            angles[c] = theta
        end
        c += 1
    end
    return Q,angles[1:(c-1)],diheds[1:(c-1)]
end

function simulatedAnnealing(Q::PolygonalNew,minFunc::Function,tolerance::Real=1e-2,thetamax::Real=pi/2,thetamin::Real=-pi/2; temp_init = 2,max_iter::Integer=1000)
    # preprocessing
    if minFunc == distToFlat
        auxQ = flatten(Q)
        minFunc = Q -> distToFlat(Q,auxQ)
    end
    diheds = zeros(Int8,max_iter)
    angles = zeros(T,max_iter)
    nq = length(Q)
    #println("bien aca")
    temp = temp_init
    delta_temp = (temp_init - 1e-6)/max_iter
    #println("intermedio")
    d = minFunc(Q)
    #println("mal aca")
    c = 1
    while d > tolerance && c <= max_iter
        theta = rand()*(thetamax-thetamin) + thetamin
        dihed = rand(1:(nq-2))
        #println(c)
        #println()
        #println(i)
        #println(Q)
        #println(linkLengths(Q))
        newQ = moveBeforeDihedral(Q,dihed)
        #println(newQ)
        #println(linkLengths(newQ))
        inter_flag = checkRotationIntersection(newQ,dihed,theta)
        newQ = dihedralRotate(newQ,dihed,theta)
        dnew = minFunc(newQ)
        #println(c)
        #println(inter_flag)
        if !inter_flag && dnew < d
            Q = newQ
            d = dnew
            diheds[c] = dihed
            angles[c] = theta
        elseif !inter_flag && dnew >= d
            r = log(rand())
            p = (-dnew + d)/temp
            if r < p
                Q = newQ
                d = dnew
                diheds[c] = dihed
                angles[c] = theta
            end
        end
        c += 1
        temp -= delta_temp
    end
    return Q,angles[1:(c-1)],diheds[1:(c-1)]
end


