include("intersections.jl")
include("io.jl")

using Statistics

function moveBeforeDihedral(P::PolygonalChain2,i::Integer)
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


function tryDihedralRotation!(P::PolygonalChain2,i::Integer,theta::Real)
    n = length(P)
    Q = moveBeforeDihedral(P,i)
    if 1 <= i < n-1
        # TODO change to cells and make brute force
        b = checkIntersection(Q[1:i],Q[i+1:end],theta)
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

function stair(n::Integer)::PolygonalChain2
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
    return PolygonalChain2(vertices)
end

function knittingneedle(l::Real=2.0;ep::Real=1/6)
    p0 = ex - ep*ey
    p1 = rotate(ez,-pi/6,ex)
    p2 = e0
    p3 = ex
    p4 = ex + ez
    p5 = ep*ey
    p0 = p1 + l*unitVector(p0-p1)
    p5 = p4 + l*unitVector(p5-p4)
    points = [p0,p1,p2,p3,p4,p5]
    return PolygonalChain2(points)
end


function randomSearch(P::PolygonalChain2,Q::PolygonalChain2,tolerance::Real=1e-2,thetamax::Real=pi/2,thetamin::Real=-pi/2;max_iter::Integer=1000,log::Bool=false)
    np = length(P)
    nq = length(Q)
    diheds = zeros(Int8,max_iter)
    angles = zeros(T,max_iter)
    #println("bien aca")
    if np == nq
        #println("intermedio")
        d = overlapedRmsd(P,Q)
        #println("mal aca")
        c = 1
        while d > tolerance && c <= max_iter
            theta = rand()*(thetamax-thetamin) + thetamin
            dihed = rand(1:(np-2))
            #println(c)
            #println()
            #println(i)
            #println(Q)
            #println(linkLengths(Q))
            newQ = moveBeforeDihedral(Q,dihed)
            #println(newQ)
            #println(linkLengths(newQ))
            b = checkIntersection(newQ,dihed,theta)
            newQ = dihedralRotate(newQ,dihed,theta)
            dnew = overlapedRmsd(P,newQ)
            #println(c)
            #println(b)
            #println()
            if !b && dnew < d
                Q = newQ
                d = dnew
                diheds[c] = dihed
                angles[c] = theta
            end
            c += 1
        end
        return Q,angles[1:(c-1)],diheds[1:(c-1)]
    else
        error("Chains must be of same length")
    end
end

function lsimulation(ls,iter::Integer,angmax::Real=pi/20,angmin::Real=-pi/20;savename="")
    n = length(ls)
    rmsds_mean = zeros(n)
    rmsds_error = zeros(n)
    for i in 1:n
        temp_rmsds = zeros(iter)
        println(i)
        for j in 1:iter
            P = stair(5)
            Q = knittingneedle(ls[i])
            #println("creacion ok")
            lastQ, angles, diheds = randomSearch(P,Q,1e-5,angmax,angmin,max_iter=10_000)
            if !isempty(savename)
                n1zeros = Int(ceil(log10(length(ls)+1)))
                n1 = lpad(i,n1zeros,'0')
                n2zeros = Int(ceil(log10(iter+1)))
                n2 = lpad(j,n2zeros,"0")
                saveSimulation(string(savename,n1,"_",n2,),P,Q,lastQ,angles,diheds,saveTrajec=false)
            end
            temp_rmsds[j] = overlapedRmsd(P,lastQ)
        end
        rmsds_mean[i] = Statistics.mean(temp_rmsds)
        rmsds_error[i] = Statistics.std(temp_rmsds)
    end
    return ls,rmsds_mean,rmsds_error
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    const savename = "../extra/prueba/"
    const ls = LinRange(1.35,1.45,21)
    if !isempty(savename)
        saveLtable(string(savename,"l"),ls)
    end
    const ls,rmsds_mean,rmsds_error = lsimulation(ls,5,savename=savename)
    scatter(ls,rmsds_mean,yerror=rmsds_error,ylabel="RMSD sobrepuesto",xlabel="l",label=false)
    savefig("../extra/prueba/lsimul.png")
end