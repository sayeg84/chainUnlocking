include("algorithms.jl")

function testRotation(n,epsilon=1e-6)
    flag = true
    for i in 1:n
        p = Point()
        theta = 2*pi*rand()
        u = unitVector(Point())
        v = rotate(p,theta,u)
        w = rotateRodrigues(p,theta,u)
        b = distance(v,w) < epsilon
        if !b
            println(v)
            println(w)
        end
        flag = flag && (b)
    end
    return flag
end

function testDihedral(n,epsilon=1e-6)
    flag = true
    for i in 1:n
        a = Point()
        b = Point()
        c = Point()
        d = Point()
        phi1 = dihedralSO(a,b,c,d)
        phi2 = dihedral(a,b,c,d)
        b = abs(phi1-phi2) < epsilon
        if !b
            println(phi1)
            println(phi2)
        end
        flag = flag && (b)
    end
    return flag
end


function singleIntersection()
    p1 = Point()
    p2 = Point()
    t = rand()
    interPoint3D = p1 + t*(p2-p1)
    vq = Point()
    q1 = interPoint3D - 0.5*vq
    q2 = interPoint3D + 0.5*vq
    flag,p = intersectionJulia(p1,p2,q1,q2)
    return flag && distance(interPoint3D,p) < 1e-5
end

function singleIntersection2()
    p1 = Point()
    p2 = Point()
    t = rand()
    interPoint3D = p1 + t*(p2-p1)
    vq = Point()
    q1 = interPoint3D - 0.5*vq
    q2 = interPoint3D + 0.5*vq
    flag,p = intersectionStaticArrays(p1,p2,q1,q2)
    return flag && distance(interPoint3D,p) < 1e-5
end

function singleIntersection3()
    p1 = Point()
    p2 = Point()
    t = rand()
    interPoint3D = p1 + t*(p2-p1)
    vq = Point()
    q1 = interPoint3D - 0.5*vq
    q2 = interPoint3D + 0.5*vq
    flag, p = intersection(p1,p2,q1,q2)
    return flag && distance(interPoint3D,p) < 1e-5
end


function intersectionTest(n)
    res = true
    for i in 1:n
        p1 = Point()
        p2 = Point()
        t = rand()
        interPoint3D = p1 + t*(p2-p1)
        vq = Point()
        q1 = interPoint3D - 0.5*vq
        q2 = interPoint3D + 0.5*vq
        flag,p = intersection(p1,p2,q1,q2)
        if !flag
            println("debug")
            println(p1)
            println(p2)
            println(q1)
            println(q2)
            println(t)
            println(interPoint3D)
            println(p)
        end
        res = res && (flag && distance(interPoint3D,p) < 1e-5)
    end
    return res
end

function testSegmentPlaneIntersection()
    a = Point()
    b = Point(a.x,a.y,-a.z)
    u = ez
    inter = Point(a.x,a.y,0.0)
    
    randomRot = rotation(2*pi*rand(),unitVector(Point()))
    a = randomRot*a
    b = randomRot*b
    u = randomRot*u
    inter = randomRot*inter
    i,p = segmentPlaneIntersection(a,b,u)
    return distance(inter,p)
end

function testConversion()
    lengths = [rand(),rand()+1,rand()+2]
    angles = [pi*rand(),pi*rand()]
    dihed = (2*rand()-1)*pi
    a = Point(0.0,0.0,0.0)
    b = Point(lengths[1],0.0,0.0)
    c = b + lengths[2]*Point(cos(pi-angles[1]),sin(pi-angles[1]),0.0)
    aux1 = Point(cos(2*pi-(angles[1]+angles[2])),sin(2*pi-(angles[1]+angles[2])),0.0)
    aux2 = Point(cos(pi-angles[1]),sin(pi-angles[1]),0.0)
    d = rotate(aux1,dihed,aux2)
    d = c + lengths[3]*d
    P = PolygonalOld((a,b,c,d))
    arr = lengthsAndAngles(P)
    Q = PolygonalChainRosetta(arr[1],arr[2],arr[3])
    optimalRotation(P,Q)
    R = (norm(angles - arr[2]) < 1e-6) && (norm(angles - arr[2]) < 1e-6) && (abs(dihed - arr[3][1]) < 1e-6)
    S = rmsd(P,Q) < 1e-6
    if R && S
        return true
    else
        println((lengths,angles,dihed))
        println(arr)
        println()
        println(P.vertices)    
        println(Q)
        return false
    end
end

function testPolygonal(n)
    lengths = rand(T,n)
    angles = pi*rand(T,n-1)
    dihed = pi*([2*rand(T)-1 for i in 1:n-2])
    P = PolygonalOld(lengths,angles,dihed)
    Q = PolygonalChainRosetta(lengths,angles,dihed)
    e = rmsd(P,Q)
    println(e)
    if e < 1e-3
        return true
    else
        println(e)
        println(lengths)
        println(angles)
        println(dihed)
        arr1 = lengthsAndAngles(P)
        arr2 = lengthsAndAngles(Q)
        display(toArray(P))
        println()
        display(arr1)
        println()
        display(toArray(Q))
        println()
        display(arr2)
        println()
        println()
        display(arr2[2]+arr1[2])
        display(arr2[2]+arr1[2])
        println()
        return false
    end
end

function benchmarkCoefficientsAldo()
    q1 = Point()
    q2 = Point()
    return coefficientsAldo(q1,q2)
end

function benchmarkCoefficientsSO()
    q1 = Point()
    q2 = Point()
    return coefficientsSO(q1,q2)
end

function testCoefficients()
    q1 = Point()
    q2 = Point()
    coefAldo = surfaceCoefficientsAldo(q1,q2)
    coefSO = surfaceCoefficientsSO(q1,q2)
    if distance(coefAldo,coefSO) > 1e-6
        println(q1)
        println(q2)
        println()
        println(coefAldo)
        println(coefSO)
        return false
    end
    return true
end

function testSurfaceQuadric()
    p1 = Point()
    p2 = Point()
    z2,z1,z0 = rand(3)
    roots1 = surfaceQuadraticEquation(z2,z1,z0,p1,p2)
    roots2 = surfaceQuadraticEquation(z2,z1,z0,p1,p2)
    if abs(roots1[2]- roots2[2]) + abs(roots1[3]- roots2[3]) > 1e-6
        println(p1)
        println(p2)
        println()
        println(roots1)
        println(roots2)
        return false
    end
    return true
end

function testSurfaceIntersection()
    q1 = Point()
    q2 = Point()
    theta = pi/2*rand()
    zmax = abs(q1.z) > abs(q2.z) ? abs(q1.z) : abs(q2.z)
    #theta2 = pi/2*rand()
    #println(rad2deg(theta2))
    #p1 = rotate(1.5*ex,theta2,ez)
    #p2 = p1 + 3*ez
    p1 = Point((q1.x+q2.x)/2,(q1.y+q2.y)/2,0.0)
    p1 = rotate(p1,theta/2,ez)
    p1 = p1 - zmax*ez
    p2 = p1 + 2*zmax*ez
    int1,int2 = surfaceLineIntersection(p1,p2,q1,q2,theta)
    println(int1)
    println(int2)
end


function testClosest()
    c1 = Point(1.0,1.0,1.0)
    q1 = 2*Point() - c1
    q2 = 2*Point() - c1
    r1,dr1,r2,dr2 = closestAndFurthestPointToOrigin(q1,q2,q2-q1)
    s1,ds1,s2,ds2 = closestAndFurthestPointToOriginAlt(q1,q2,q2-q1)
    if distance(r1,s1) < 1e-6 && distance(r2,s2) < 1e-6
        return true
    else
        println(q1)
        println(q2)
        println(r1)
        println(r2)
        println(s1)
        println(s2)
        return false
    end
end

function benchmarkClosest1()
    c1 = Point(1.0,1.0,1.0)
    q1 = 2*Point() - c1
    q2 = 2*Point() - c1
    return closestAndFurthestPointToOrigin(q1,q2,q2-q1)
end

function benchmarkClosest2()
    c1 = Point(1.0,1.0,1.0)
    q1 = 2*Point() - c1
    q2 = 2*Point() - c1
    return closestAndFurthestPointToOriginAlt(q1,q2,q2-q1)
end

function testThirdCaseIntersection()
    #q1 = Point(rand(),rand(),0)
    #q2 = Point(rand(),rand(),0)
    q1 = e0
    q2 = ex
    theta = pi/2*rand()
    aux = (q1+q2)/2
    aux = rotate(aux,theta/2,ez)
    p1 = aux 
    p2 = aux + rotate(aux,pi/100*(2*rand()-1),ez)
    vp = p2-p1
    vq = q2-q1
    th = pi/4
    #println(bsms(q1,q2,vq,th))
    #println(p1)
    #println(vp)
    #println(eqs0(q1,q2,vq,p1,p2,vp,th))
    println()
    #println("rotation angle")
    #println(theta)
    return thirdCaseIntersection(p1,p2,vp,q1,q2,vq,theta)
end

function plotCircle!(c::Point,r::Real;lw::Real=1)
    xs = [r*cos(t) + c.x for t in LinRange(0,2*pi,101)]
    ys = [r*sin(t) + c.y for t in LinRange(0,2*pi,101)]
    plot!(xs,ys,lw=lw,aspect_ratio=1)
end

function plotLine!(p1::Point,p2::Point;lw::Real=1,label::String="")
    vp = p2-p1
    line = [p1+t*vp for t in LinRange(0,1,101)]
    line = [toArray(p) for p in line]
    line = hcat(line...)
    if label == ""
        plot!(line[1,:],line[2,:],lw=lw) 
    else
        plot!(line[1,:],line[2,:],lw=lw,label=label)
    end
end

function plotRotation!(q1::Point,q2::Point,theta::Real)
    vq = q2-q1
    angles = LinRange(0,theta,200)
    ts = LinRange(0,1,3)
    surf = [zrotation(a)*(q1 + t*vq) for a in angles for t in ts ]
    surf = [toArray(p) for p in surf]
    surf = hcat(surf...)
    plot!(surf[1,:],surf[2,:],lw=1,color="black",alpha=0.5)
    r1,d1,r2,d2 = closestAndFurthestPointToOrigin(q1,q2,vq)
    scatter!([r1.x,r2.x],[r1.y,r2.y])
    plotCircle!(e0,d1)
    plotCircle!(e0,d2)
end

function plotChain!(P::AbstractChain,c,lw=2)
    arr = toArray(P)
    scatter!(arr[:,1],arr[:,2],arr[:,3],color=c,xlabel="x",ylabel="y",zlabel="z",marker_z=1:(length(P)+1),seriescolor=:darkrainbow,markersize=1)
    plot!(arr[:,1],arr[:,2],arr[:,3],color="black",lw=lw)
end

function testParallelCheck(n)
    flag = true
    for i in 1:n
        p = plot()
        q1 = Point(2*rand()-1,2*rand()-1,0.0) 
        q2 = Point(2*rand()-1,2*rand()-1,0.0) 
        theta = pi/2*(2*rand()-1)
        p1 = zrotation(theta/3)*q1
        p2 = zrotation(theta/3)*q2
        vp = p2 - p1
        p1 = p1 + 0.6*vp #+ 0.1*ex
        p2 = p2 + 0.6*vp
        a = thirdCaseIntersection(p1,p2,p2-p1,q1,q2,q2-q1,theta)
        flag = flag && a
        if !a
            println(p1)
            println(p2)
            println(q1)
            println(q2)
            println(theta)
            plotRotation!(q1,q2,theta)
            plotLine!(p1,p2,lw=3,label="segmento")
            display(p)
        end
    end
    return flag
end

function benchmarkC1Quadratic1()
    p1 = Point()
    p2 = Point()
    q1 = Point()
    q2 = Point()
    theta = pi/4*(2*rand()-1)
    c,r,d = case1SurfaceCoefficients(q1,q2,q2-q1)
    return case1SurfaceQuadricRoots(c,r,d,p1,p2,q2-q1)
end

function benchmarkC1Quadratic2()
    p1 = Point()
    p2 = Point()
    q1 = Point()
    q2 = Point()
    theta = pi/4*(2*rand()-1)
    c,r,d = case1SurfaceCoefficients(q1,q2,q2-q1)
    return case1SurfaceQuadricRoots2(c,r,d,p1,p2,q2-q1)
end

function assertRootsEqual(n)
    flag = true
    for i in 1:n
        p1 = Point()
        p2 = Point()
        q1 = Point()
        q2 = Point()
        theta = pi/4*(2*rand()-1)
        c,r,d = case1SurfaceCoefficients(q1,q2,q2-q1)
        A =  case1SurfaceQuadricRoots(c,r,d,p1,p2,q2-q1)
        B =  case1SurfaceQuadricRoots(c,r,d,p1,p2,q2-q1)
        flag = flag && (A[1] == B[1] && isapprox(A[2],B[2],atol=1e-15) && isapprox(A[3],B[3],atol=1e-15))
    end
    return flag
end




if abspath(PROGRAM_FILE) == @__FILE__
    using BenchmarkTools, IntervalArithmetic, IntervalRootFinding
    println()
    println("entramos")
    println()  
    #println(testThirdCaseIntersection())
    display(@benchmark benchmarkC1Quadratic1())
    println()
    println()
    display(@benchmark benchmarkC1Quadratic2())
    println()
    println(assertRootsEqual(1000))
end