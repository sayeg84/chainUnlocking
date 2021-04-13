include("intersections.jl")

function moveBeforeDihedral(P::PolygonalChain2,i::Int)
    newP = copy(P)
    n = length(P)
    if 1 <= i <= n-2
        u = unitVector(newP.endpoints[i+2]-newP.endpoints[i+1])
        mat = planeRotationXY(u)
        pivot = newP.endpoints[i+2]
        for j in 1:n+1
            newP.endpoints[j] = mat*(newP.endpoints[j]-pivot)
        end
        return newP
    else
        error("i is not in range of dihedrals")
    end
end


function tryDihedralRotation!(P::PolygonalChain2,i::Int,theta::Real)
    n = length(P)
    Q = moveBeforeDihedral(P,i)
    if 1 <= i < n-1
        # TODO change to cells and make brute force
        b = checkIntersection(Q.endpoints[1:i],Q.endpoints[i+1:end],theta)
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

function stair(n::Int)::PolygonalChain2
    #endpoints = Array{Point,1}(undef,n+1)
    endpoints = [e0 for i in 1:n+1]
    endpoints[1] = e0
    endpoints[2] = ex
    for i in 3:(n+1)
        if mod(i,2) == 1
            endpoints[i] += endpoints[i-1]+ ey
        else
            endpoints[i] += endpoints[i-1]+ ex
        end
    end
    return PolygonalChain2(endpoints)
end


function randomSearch(P::PolygonalChain2,Q::PolygonalChain2,tolerance::Real=1e-2,thetamax::Real=pi/2,thetamin::Real=-pi/2;max_iter::Integer=1000)
    np = length(P)
    nq = length(Q)
    angle_indexes = zeros(Int8,max_iter)
    angle_values = zeros(T,max_iter)
    if np == nq
        d = overlapedRmsd(P,Q)
        c = 1
        while d > tolerance && c <= max_iter
            theta = rand()*(thetamax-thetamin) + thetamin
            i = rand(1:(np-2))
            #println(c)
            #println()
            #println(i)
            #println(Q)
            #println(linkLengths(Q))
            newQ = moveBeforeDihedral(Q,i)
            #println(newQ)
            #println(linkLengths(newQ))
            b = checkIntersection(newQ.endpoints[1:i],newQ.endpoints[i+1:end],theta)
            newQ = dihedralRotate(newQ,i,theta)
            dnew = overlapedRmsd(P,newQ)
            if !b && dnew < d
                Q = newQ
                d = dnew
                angle_indexes[c] = i
                angle_values[c] = theta
            end
            c += 1
        end
        return Q,angle_values[1:(c-1)],angle_indexes[1:(c-1)]
    else
        error("Chains must be of same length")
    end
end