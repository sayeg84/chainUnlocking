include("scalarRootFinding.jl")
include("geomTypes.jl")



"""
intervalOverlap(a::RealOrDual,b::RealOrDual,c::RealOrDual,d::RealOrDual)::Bool

Function to check wheter the intervals of real numbers between a,b and c,d have an intersection.
"""
function intervalOverlap(a::RealOrDual,b::RealOrDual,c::RealOrDual,d::RealOrDual)::Bool
    a,b = minmax(a,b)
    c,d = minmax(c,d)
    return (b >= c) && (d >= a)
end


"""
segmentOverlap(p1::Point,p2::Point,q1::Point,q2::Point)::Bool

Function to check wheter segments on the same line `p1p2` and `q1q2` intersect
"""

function segmentOverlap(p1::Point,p2::Point,q1::Point,q2::Point)::Bool
    xtest = intervalOverlap(p1.x,p2.x,q1.x,q1.x)
    ytest = intervalOverlap(p1.y,p2.y,q1.y,q1.y)
    ztest = intervalOverlap(p1.z,p2.z,q1.z,q1.z)
    # segment overlaps if all their x, y and z intervals overlap
    return xtest && ytest && ztest
end





"""
quadraticEquationRoots(a::RealOrDual,b::RealOrDual,c::RealOrDual)

Returns a triple `n,r1,r2` containing the solutions of equation ``a x^2 + bx + c = 0``. `n` is the number of different roots that exist of the equation.
"""
function quadraticEquationRoots(a::RealOrDual,b::RealOrDual,c::RealOrDual)
    R = pow2(b) - 4*a*c
    A = isapprox(a,0,atol=1e-15)
    B = isapprox(b,0,atol=1e-15)
    C = isapprox(c,0,atol=1e-15)
    # possible infinite solutions    
    if (A && B && C) 
        return -1,0,0
    # no solutions
    elseif (A && B) || R < -1e-8
        return 0,0,0
    # equation is linear
    elseif isapprox(a,0,atol=1e-15)
        return 1, -c/b, -c/b
    # two roots are equal
    elseif -1e-8 <= R < 1e-15
        x = -b/(2*a)
        return 1,x,x
    # general case
    else
        v = sqrt(R)
        return 2, (-b - v)/(2*a), (-b + v)/(2*a)
    end
end


"""
linearIndependence(p::Point,q::Point)::Bool

Uses Julia's LinearAlgebra.rank to check if points p and q are linearly independent.
"""
function linearIndependence(p::Point,q::Point)::Bool
    A = hcat(to2DArray(p),to2DArray(q))
    return LinearAlgebra.rank(A) == 2
end


"""
xyIangle(u::Point,v::Point)

Finds the internal angle between the xy proyection of points `u` and `v`.
"""
function xyIangle(u::Point,v::Point)
    sint = u.x*v.y - u.y*v.x
    cost = u.x*v.x + u.y*v.y
    return atan(sint,cost)
end

#### intersections

"""
segmentIntersectionJulia(p1::Point,p2::Point,q1::Point,q2::Point)::Tuple{Bool,Point}

Finds the intersection point between two segments p1,p2 and q1,q2 by solving the linear equation system arising using Julia's `\\` operator. 

It returns a tuple regarding with a bool telling if there is intersection and the intersection point.
"""
function segmentIntersectionJulia(p1::Point,p2::Point,q1::Point,q2::Point)::Tuple{Bool,Point}
    vp = p2 - p1
    vq = q2 - q1
    B = [q1.x-p1.x , q1.y-p1.y , q1.z-p1.z]
    A = [vp.x -vq.x ; vp.y -vq.y ; vp.z -vq.z ]
    ts = A \ B
    if 0 < ts[1] < 1 && 0 < ts[2] < 1
        return true , p1 + ts[1]*vp
    else
        return false, Point(0,0,0)
    end
end

"""
segmentIntersectionStaticArrays(p1::Point,p2::Point,q1::Point,q2::Point)::Tuple{Bool,Point}

Finds the intersection point between two segments p1,p2 and q1,q2 by solving the linear equation system arising using Julia's `\\` operator but over staticArrays.

It returns a tuple regarding with a bool telling if there is intersection and the intersection point.
"""
function segmentIntersectionStaticArrays(p1::Point,p2::Point,q1::Point,q2::Point)::Tuple{Bool,Point}
    vp = p2 - p1
    vq = q2 - q1
    B = SVector{3,typeof(p1.x)}([q1.x-p1.x , q1.y-p1.y , q1.z-p1.z])
    A =  SMatrix{3,2,typeof(p1.x)}([vp.x -vq.x ; vp.y -vq.y ; vp.z -vq.z ])
    ts = A \ B
    if 0 < ts[1] < 1 && 0 < ts[2] < 1
        return true , p1 + ts[1]*vp
    else
        return false, Point(0,0,0)
    end
end


"""
twoxtwoLinearSystem(a11::RealOrDual,a12::RealOrDual,b1::RealOrDual,a21::RealOrDual,a22::RealOrDual,b2::RealOrDual)

Solves the 2x2 linear system defined by a11x + a12y = b1 and a21x + a22y = b2.

If the system is singular, it returns `0,0,0`. 

If there is an unique solution, it returns `1,x,y`
"""
function twoxtwoLinearSystem(a11::RealOrDual,a12::RealOrDual,b1::RealOrDual,a21::RealOrDual,a22::RealOrDual,b2::RealOrDual)
    det = a11*a22 - a12*a21
    if !isapprox(det,0,atol=1e-15)
        detx = b1*a22 - b2*a12
        dety = a11*b2 - a21*b1
        return 1, detx/det, dety/det
    else
        return 0, 0, 0
    end
end

"""
segmentIntersection(p1::Point,p2::Point,q1::Point,q2::Point)::Tuple{Bool,Point}

Checks if the segments defined by lines p1-p2 and q1-q2 intersect and returns a flag and the intersection point.

The intersection defines a the 3x2 linear system, which is solved by decomposing it into three different 2x2 linear systems which can be fastly solved. 
Note that some of the 2x2 systemscan be singular. After finding the solution, the function will compare the solutions between all the nonsingular systems to see it is the same.

Returns `false,e0` if there is not intersection and `true,p` if there is intersection at a point p
"""
function segmentIntersection(p1::Point,p2::Point,q1::Point,q2::Point)::Tuple{Bool,Point}
    vp = p2 - p1
    vq = q2 - q1
    # linear systems to determine parametric point of intersection.
    # returns triple of boolean, and tp and tq if such an intersection exists
    # x-y system
    ts1 = twoxtwoLinearSystem(vp.x,-vq.x,q1.x-p1.x,vp.y,-vq.y,q1.y-p1.y)
    # x-z system
    ts2 = twoxtwoLinearSystem(vp.x,-vq.x,q1.x-p1.x,vp.z,-vq.z,q1.z-p1.z)
    # y-z system
    ts3 = twoxtwoLinearSystem(vp.y,-vq.y,q1.y-p1.y,vp.z,-vq.z,q1.z-p1.z)
    ts = (ts1,ts2,ts3)
    # all systems are singular
    if isapprox(ts[1][1] + ts[2][1] + ts[3][1],0,atol=1e-15)
        return false, e0
    else
        # finding sum of distances between solutions of non-singular systems
        maxDist = 0
        for i in 1:2
            for j in i+1:3
                maxDist += ts[i][1]*ts[j][1]*distance(ts[i][2:3],ts[j][2:3])
            end
        end
        if maxDist < 1e-6
            if 0 < ts1[2] < 1 && 0 < ts1[3] < 1
                return true , p1 + ts1[2]*vp
            else
                return false, e0
            end
        else
            return false, e0
        end
    end
end

"""
segmentPlaneIntersection(p1::Point,p2::Point,u::Point)::Tuple{Int16,Point}

Finds the intersection of a segment defined by points `p1`,`p2` and a plane that goes thorugh the origin with normal vector `u`
"""
function segmentPlaneIntersection(p1::Point,p2::Point,u::Point)::Tuple{Int16,Point}
    v = p2 - p1
    discrim = dot(v,u)
    num = dot(-1*p1,u)
    t = num/discrim
    if isapprox(discrim, 0,atol=1e-15)
        if isapprox(num, 0,atol=1e-15)
            # line is contained in plane
            return 1, e0
        else
            # line and plane are parallel
            return 0, e0
        end
    elseif  0 <= t <= 1
        return 2, p1 + t*v
    else
        return 0, e0
    end
end



function angleTest(ang::RealOrDual,alpha::RealOrDual)::Bool
    return alpha > 0 ? (0 <= ang <= alpha) : (0 >= ang >= alpha)
end

#=
We want to find the intersection between a rotating 3D segment, denoted q1,q2 and another segment, denoted p1,p2. For simplicity, we will make appropiate transfomations so that the rotation is around the z axis.

In that case, we must distinguish three important cases that can arise between the segments


* Case 1: Neither of the segments is contained in a plane parallel to the XY plane. 

In this case, to find the intersection, we must use cylindrical coordinates to describe the surface of revolution arising from the rotation of segment q and solve a quadratic equation arising. The only troublesome case would be

** Subcase 1.1: the intersection between the rotating segment and the non rotating happens in an infinite set of points due to the segments becoming the same at one point


* Case 2: Segment q is in a plane parallel to the XY plane but segment p isn't. 

In this case, we must find the intersection of the segment p in the plane paralell to XY containing segment q and then test if its point is inside the region using the same parametrization than in case 1, which can be easily done. 


* Case 3: Segment q and segment p are in a plane parallel to the XY plane. 

In this case, the intersection only ocurrs if the plane is the same for both segments. If it is the same, we must test for the intersections by checking if the circles decribed by the rotations of vertices q1,q2 intersect the p1,p2 segment in apropriate regions
=#

### Case 1 functions

"""
case1SurfaceCoefficients(q1::Point,q2::Point,vq::Point)

Returns coeficients `d,r,c` such that the surface of revolution made by rotating segment q1q2 around the zaxis is described by ``r^2(z) = c z^2 + 2 r z + d``

For more information, check [this math SO post](https://math.stackexchange.com/questions/4082142/intersection-between-rotating-3d-line-and-3d-line)
"""
function case1SurfaceCoefficients(q1::Point,q2::Point,vq::Point)
    aux_1 = pow2(1/vq.z)
    aux_2 = 1/vq.z
    aux_3 = (pow2(vq.x) + pow2(vq.y))
    aux_4 = (q1.x*q2.z - q2.x*q1.z)
    aux_5 = (q1.y*q2.z - q2.y*q1.z)
    c = aux_1*aux_3
    r = aux_1*(vq.x*aux_4 + vq.y*aux_5)
    d = aux_1*(pow2(aux_4) + pow2(aux_5))
    return c,r,d
end

#=
function case1surfaceCoefficientsAldo(q1::Point,q2::Point)::Tuple{T,T,T}
    vq = q2-q1
    aux_1 = pow2(1/vq.z)
    aux_2 = 1/vq.z
    aux_3 = (pow2(vq.x) + pow2(vq.y))
    aux_4 = (q1.x*vq.x + q1.y*vq.y)
    c = aux_3*aux_1
    r = aux_2*aux_4 - q1.z*c
    d = q1.z*(q1.z*aux_1*aux_3 - 2*aux_4*aux_2) + pow2(q1.x) + pow2(q1.y)
    return c,r,d
end


function surfaceQuadraticEquation(c::RealOrDual,r::RealOrDual,d::RealOrDual,p1::Point,p2::Point)
    vp = p2-p1
    S = p1.x*vp.x + p1.y*vp.y - (r+c*p1.z)*vp.z
    D = c*pow2(vp.z) - pow2(vp.x) - pow2(vp.y)
    aux_1 = p1.x*p2.z - p2.x*p1.z
    aux_2 = p1.y*p2.z - p2.y*p1.z
    aux_3 = p1.x*p2.y - p2.x*p1.y
    R = c*(pow2(aux_1) + pow2(aux_2) - d*pow2(vp.z)) 
    R += pow2(r*vp.z) - pow2(aux_3) 
    R += 2*r*(-aux_1*vp.z-aux_2*vp.y)
    R += d*(pow2(vp.x) + pow2(vp.y))
    if R < 0 || D == 0
        return false,0,0
    elseif R == 0
        return true,S/D,S/D
    else
        return true, (S - sqrt(R))/D, (S + sqrt(R))/D
    end
end
=#




"""
function case1SurfaceQuadricRoots(c::RealOrDual,r::RealOrDual,d::RealOrDual,p1::Point,p2::Point,vp::Point)::Tuple{Int16,T,T}

Given a revolution surface described by ``r^2(z) = c z^2 + 2 r z + d``, it returns the roots ``t`` of the equation arising for substituing ``x^2(t) + y^2(t)= r^2(z(t))``.

It returns `a,r1,r2` where `a`is an Int16 representing the number of roots, `r` the first root and `r2`. It returns `-1,0,0` when the equation is poorly defined due to an infinite number of roots.
"""
function case1SurfaceQuadricRoots(c::RealOrDual,r::RealOrDual,d::RealOrDual,p1::Point,p2::Point,vp::Point)
    cp = cross(p1,p2)
    a = pow2(vp.x) +  pow2(vp.y) - c*pow2(vp.z)
    b = 2*(p1.x*vp.x + p1.y*vp.y - vp.z*(c*p1.z + r))
    c = pow2(p1.x) + pow2(p1.y) - c*pow2(p1.z) - 2*r*p1.z - d
    return quadraticEquationRoots(a,b,c)
end


"""
function case1SurfaceQuadricRootsAlt(c::RealOrDual,r::RealOrDual,d::RealOrDual,p1::Point,p2::Point,vp::Point)::Tuple{Int16,T,T}

Given a revolution surface described by ``r^2(z) = c z^2 + 2 r z + d``, it returns the roots ``t`` of the equation arising for substituing ``x^2(t) + y^2(t)= r^2(z(t))``.

It returns `a,r1,r2` where `a`is an Int16 representing the number of roots, `r` the first root and `r2`. It returns `-1,0,0` when the equation is poorly defined due to an infinite number of roots.

It computes the coefficientes of the equation using a simplified form obtained from `https://math.stackexchange.com/questions/4082142/intersection-between-rotating-3d-line-and-3d-line/`
"""
function case1SurfaceQuadricRootsAlt(c::RealOrDual,r::RealOrDual,d::RealOrDual,p1::Point,p2::Point,vp::Point)
    cp = cross(p1,p2)
    S = p1.x*vp.x + p1.y*vp.y - (r+c*p1.z)*vp.z
    D = c*pow2(vp.z) - pow2(vp.x) - pow2(vp.y)
    R = c*(pow2(cp.y) + pow2(cp.x) - d*pow2(vp.z)) 
    R += pow2(r*vp.z) - pow2(cp.z) 
    R += 2*r*(cp.y*vp.x - cp.x*vp.y)
    R += d*(pow2(vp.x) + pow2(vp.y))
    A = isapprox(D,0,atol=1e-15)
    B = isapprox(R,0,atol=1e-15)
    C = isapprox(S,0,atol=1e-15) 
    # infinite solutions
    if A && B && C 
        return -1,0,0
    #  no real roots
    elseif (A && B) || R < 1e-8
        return 0,0,0
    # equation is linear and only with one root
    elseif A
        n = pow2(p1.x) + pow2(p1.y) - c*pow2(p1.z) - 2*r*p1.z - d
        return 1, -n/sqrt(R),0
    # two roots are equal
    elseif -1e-8 <= R < 1e-15
        return 1, S/D, S/D
    # two different roots
    else
        return 2, (S - sqrt(R))/D, (S + sqrt(R))/D
    end
end


"""
alphaParameter(xp::RealOrDual,yp::RealOrDual,xq::RealOrDual,yq::RealOrDual)

returns the value of `alpha` solving the equation system:

xq cos(alpha) + yq sin(alpha) = xp 
-xq sin(alpha) + yq cos(alpha) = yp
"""
function alphaParameter(xp::RealOrDual,yp::RealOrDual,xq::RealOrDual,yq::RealOrDual)
    sinalpha = xq*yp - xp*yq
    cosalpha = xp*xq + yp*yq
    return atan(sinalpha,cosalpha)
end


"""
getParameters(t::RealOrDual,p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)

Given a value t such that the following equations hold 

q1.z + s*vq.z = p1.z + t*vp.z
(q1.x + s*vq.x) cos(alpha) + (q1.y + s*vq.y) sin(alpha) = p1.x + t*vp.x
-(q1.x + s*vq.x) sin(alpha) + (q1.y + s*vq.y) cos(alpha) = p1.y + t*vp.y

It returns the tuple of `s,alpha` such that the equations hold.
"""
function getParameters(t::RealOrDual,p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)
    inp = p1 + t*vp
    s = (inp.z- q1.z)/vq.z
    inq = q1 + s*vq
    alpha = alphaParameter(inp.x,inp.y,inq.x,inq.y)
    return s,alpha
end


"""
case1CheckRoot(t::RealOrDual,p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual)

Given a parameter `t` which is root of the equation ``x^2(t) + y^2(t)= r^2(z(t))`` arising from rotation segment q1q2, check if  the point (x(t),y(t),z(t)) is inside the surface of revolution of rotating segment q1q2 an angle alpha.

It returns a tuple `a,p` where `a` is a Bool saying if there is intersection and `p` is the point of such intersection.
"""
function case1CheckRoot(t::RealOrDual,p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual;debug=false)
    # checking the intersection is within the segment
    if 0 <= t <= 1
        debug && println("# roots are fine")
        # a,b is the interval of z values for the q1-q2 segment
        a,b = minmax(q1.z,q2.z)
        inp = p1 + t*vp
        if debug
            println(inp)
            println("minmax(q1.z,q2.z) = $a,$b")
            ((s,ang) = getParameters(t,p1,p2,vp,q1,q2,vq) )
            println("s = $s")
            println("ang = $(ang)")
        end
        # check intersection only if z value is in interval
        if a <= inp.z <= b 
            s,ang = getParameters(t,p1,p2,vp,q1,q2,vq)
            if debug
                println("# in z is fine")
                println("s = $s")
                println("ang = $(ang)")
                println(0 <= s <= 1)
                println(angleTest(ang,alpha))
            end
            if 0 <= s <= 1 && angleTest(ang,alpha)
                return true, inp
            end
        end 
        return false, e0
    else
        return false, e0
    end
end

function case1Intersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual;debug=false)::Tuple{Bool,Point,Point}
    c,r,d = case1SurfaceCoefficients(q1,q2,vq)
    roots = case1SurfaceQuadricRoots(c,r,d,p1,p2,vp)
    debug && println(roots)
    # special case for when the intersection happens in an infinite amount of values
    if roots[1] == -1
        A = isapprox(vq.x,0,atol = 1e-5)
        B = isapprox(vq.y,0,atol = 1e-5)
        # cylinder and hyperboloid have angles defined in different ways
        # cylinder angle of rotation 
        ang =  A && B ? xyIangle(p1,q1) : xyIangle(vq,vp)
        if angleTest(ang,alpha)
            mat = zrotation(ang)
            q1prime = mat*q1
            q2prime = mat*q2
            val = segmentOverlap(p1,p2,q1prime,q2prime)
            return val,e0,e0
        end
    elseif roots[1] == 1
        int = case1CheckRoot(roots[2],p1,p2,vp,q1,q2,vq,alpha,debug=debug)
        return int[1],int[2],e0
    elseif roots[1] == 2
        #println("aca")
        int1 = case1CheckRoot(roots[2],p1,p2,vp,q1,q2,vq,alpha,debug=debug)
        int2 = case1CheckRoot(roots[3],p1,p2,vp,q1,q2,vq,alpha,debug=debug)
        return int1[1] || int2[1], int1[2],int2[2]
    else
        return false,e0,e0
    end
end


"""
case2coefficients(inp::Point,q1::Point,q2::Point,vq::Point)::Tuple{T,T,T}

Returns the coefficients `a,b,c` of the quadric equation arising from the substitution of point ``x^2(s) + y^2(s) = inp.x^2 + inp.y^2``
"""
function case2coefficients(inp::Point,q1::Point,q2::Point,vq::Point)
    a = pow2(vq.x) + pow2(vq.y)
    b = q1.x*vq.x + q1.y*vq.y
    c = pow2(q1.x) + pow2(q1.y) - pow2(inp.x) - pow2(inp.y)
    return a,b,c
end


"""
case2checkRoot(s::RealOrDual,q1::Point,q2::Point,vq::Point,inp::Point,alpha::RealOrDual)::Tuple{Bool,Point}

For a point `s` such that it is root of the equations arising from ``x^2(s) + y^2(s) = inp.x^2 + inp.y^2``, it finds the angle `ang` and checks if the root represents a point of intersection between the 2D surface arising from rotating segment `q1q2` and the segment `p1p2`.
"""
function case2checkRoot(s::RealOrDual,q1::Point,q2::Point,vq::Point,inp::Point,alpha::RealOrDual)::Tuple{Bool,Point}
    if 0 <= s <= 1
        inq = q1 + s*vq
        ang = alphaParameter(inp.x,inp.y,inq.x,inq.y)
        if angleTest(ang,alpha)
            return true, inq
        else
            return false, e0
        end
    else
        return false, e0
    end
end


function case2Intersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual;debug=true)::Tuple{Bool,Point,Point}
    a,b = minmax(p1.z,p2.z)
    # plane parallel to XY must be in the z interval of segment p1p2
    if a < q1.z < b 
        t = (q1.z-p1.z)/vp.z
        inp = p1 + t*vp
        a,b,c = case2coefficients(inp,q1,q2,vq)
        roots = quadraticEquationRoots(a,b,c)
        int1 = case2checkRoot(roots[2],q1,q2,vq,inp,alpha)
        int2 = case2checkRoot(roots[3],q1,q2,vq,inp,alpha)
        if int1[1] && int2[1]
            return true, int1[2], int2[2]
        elseif int1[1]
            return true, int1[2], e0
        elseif int2[1]
            return true, int2[2], e0
        else
            return false, e0, e0   
        end
    else
        return false, e0, e0
    end
end



#### third case

"""
closestAndFurthestPointToOrigin(q1::Point,q2::Point,vq::Point)::Tuple{Point,T,Point,T}

Given a segment `q1q2`, it returns the tuple `pmin,dmin,pmax,dmax` such that `pmin,dmin` is the closest point in the segment to the origin and its distance, and such that `pmax,dmax` is the furthest point in the segment to the origin and its distance
"""
function closestAndFurthestPointToOrigin(q1::Point,q2::Point,vq::Point)
    d1 = norm(q1)
    d2 = norm(q2)
    # finding closest point to origin in the infinite line containing the segment
    t = -dot(vq,q1)/dot(vq,vq)
    point = q1 + t*vq
    if 0 <= t <= 1
        dp = norm(point)
        return d1 > d2 ? (point,dp,q1,d1) : (point,dp,q2,d2)
    else 
        return d1 > d2 ? (q2,d2,q1,d1) : (q1,d1,q2,d2)
    end
end

#=
function closestAndFurthestPointToOriginAlt(q1::Point,q2::Point,vq::Point)::Tuple{Point,T,Point,T}
    t = -dot(vq,q1)/dot(vq,vq)
    point = q1 + t*vq
    d1 = norm(q1)
    d2 = norm(q2)
    dp = norm(point)
    if t < 0
        return q1, d1, q2, d2
    elseif 0 <= t < 0.5
        return point, dp, q2,d2
    elseif 0.5 <= t < 1
        return point, dp, q1, d1
    else
        return q2, d2, q1, d1
    end
end
=#


"""
xyIntersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool

Function to compute the points the XY-contained lines containing  segments `p1p2` and `q1q2` intersect
"""
function xyIntersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)
    ts1 = twoxtwoLinearSystem(vp.x,-vq.x,q1.x-p1.x,vp.y,-vq.y,q1.y-p1.y)
    return ts1
end


"""
xySegmentIntersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool

Function to check if XY-contained segments `p1p2` and `q1q2` intersect
"""
function xySegmentIntersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool
    ts1 = twoxtwoLinearSystem(vp.x,-vq.x,q1.x-p1.x,vp.y,-vq.y,q1.y-p1.y)
    return Bool(ts1[1]) && (0 < ts1[2] < 1) && (0 < ts1[3] < 1)
end



"""
segmentCircleIntersection(p1::Point,p2::Point,vp::Point,r::Point)

Finds the intersection between an XY-contained segment `p1-p2` and a circle of radius `r` centered at the origin.
"""
function segmentCircleIntersection(p1::Point,p2::Point,vp::Point,r::RealOrDual)
    a = pow2(vp.x) + pow2(vp.y)
    b = 2*(p1.x*vp.x + p1.y*vp.y)
    c = pow2(p1.x) + pow2(p1.y) - pow2(r)
    roots = quadraticEquationRoots(a,b,c)
    if 0 <= roots[2] <= 1 && 0 <= roots[3] <= 1
        return 2,roots[2],roots[3]
    elseif 0 <= roots[2]  <= 1
        return 1, roots[2], 0
    elseif 0 <= roots[3]  <= 1
        return 1, roots[3], 0
    else
        return 0,0,0
    end
end

"""
checkSegmentCircleIntersection(roots,p1::Point,p2::Point,vp::Point,q::Point,alpha::RealOrDual)::Bool

Given the `roots` tuple representing the values of parameter `t` at which intersections between a segment `p1p2` and a circle centered at the origin and passing by `q`, checks if the intersection happens when rotating point `q` an angle smaller than `alpha`
"""
function checkSegmentCircleIntersection(roots,p1::Point,p2::Point,vp::Point,q::Point,alpha::RealOrDual)::Bool
    if roots[1] == 0
        return false
    elseif roots[1] == 1 
        inp = p1 + roots[2]*vp
        #scatter!([inp.x],[inp.y])
        # finding angle between intersection and circle-defining vector
        ang = xyIangle(q,inp)
        return angleTest(ang,alpha)
    else
        inp1 = p1 + roots[2]*vp
        ang1 = xyIangle(q,inp1)
        inp2 = p1 + roots[3]*vp
        ang2 = xyIangle(q,inp2)
        #scatter!([inp1.x,inp2.x],[inp1.y,inp2.y])
        return angleTest(ang1,alpha) || angleTest(ang2,alpha)
    end
end

function case3Intersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual;debug=false)
    # segments already intersect
    # this should never happen in practice but we add it for debugging reasons
    # println("inside case 3")
    if xySegmentIntersection(p1,p2,vp,q1,q2,vq)
        debug && println("# WARNING: segments intersect before rotating")
        return true
    # only test for intersection if points are in the same interval        
    elseif isapprox(p1.z,q1.z,atol=1e-15)
        debug && println("# special case")
        if (distance(p1,q1)<1e-40 || distance(p1,q2)<1e-40 || distance(p2,q1)<1e-40 || distance(p2,q2)<1e-40)
            debug && println("# special special case: endopoints coincide")

        end
        debug && println("# first pass")
        for q in (q1,q2)
            r = norm(q)
            roots = segmentCircleIntersection(p1,p2,vp,r)
            debug && println(roots)
            #if roots[1]==2 && ((roots[2]==1 && roots[3]==1) || (roots[2]==0 && roots[3]==0))
            #    return false
            #end
            test = checkSegmentCircleIntersection(roots,p1,p2,vp,q,alpha)
            if test
                return true
            end
        end
        debug && println("# second pass")
        for p in (p1,p2)
            r = norm(p)
            roots = segmentCircleIntersection(q1,q2,vq,r)
            test = checkSegmentCircleIntersection(roots,q1,q2,vq,p,-alpha)
            if test
                return true
            end
        end
        return false
    else
        return false
    end
end

"""
firstFilter(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool

Uses closestAndFurthestPointToOrigin to find closes and furthests points of segments `p1p2` and `q1q2` and tests if they satisfy the needed inequalities for posible intersection.
"""
function firstFilter(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool
    r1,dr1,r2,dr2 = closestAndFurthestPointToOrigin(q1,q2,vq)
    s1,ds1,s2,ds2 = closestAndFurthestPointToOrigin(p1,p2,vp)
    # if largest radius in segment p1p2 is smaller than the smallest radius
    # of q1q2 or if smallest radius of p1p2 is larger than large
    if  ds2 < dr1 || ds1 > dr2
        return false
    else
        return true
    end
end

"""
bsms(q1::Point,q2::Point,vq::Point,alpha)

Returns the values `b3x,b3y,m3x,m3y` obtained of representing the rotation of segment q1q2 as the parametric curve ``(b3x + s*m3x,b3x + s*m3x)`` with `s` in ``[0,1]``
"""
function bsms(q1::Point,q2::Point,vq::Point,alpha)
    a = cos(alpha)
    b = sin(alpha)
    b3x = q1.x*a - q1.y*b
    b3y = q1.x*b + q1.y*a
    m3x = vq.x*a - vq.y*b
    m3y = vq.x*b + vq.y*a
    return b3x,b3y,m3x,m3y
end

"""
st(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)

Returns the parameters `s,t` at which the intersection of the line containing segment `p1p2` and line containing segment `q1q2` rotated an angle `alpha` happens.
"""
function st(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)
    b3x,b3y,m3x,m3y = bsms(q1,q2,vq,alpha)
    s = ((b3y-p1.y)*vp.x - (b3x-p1.x)*vp.y)/(m3x*vp.y - m3y*vp.x)
    if !isapprox(vp.x,0,atol=1e-15)
        t = (b3x + s*m3x - p1.x)/vp.x
    else
        t = (b3y + s*m3y - p1.y)/vp.y
    end
    return s,t
end


"""
eqs0(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual

Explicit equation for when parameter `s` of segment `q1q2` satisfies `s=0` 
"""
function eqs0(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual
    b3x,b3y,m3x,m3y = bsms(q1,q2,vq,alpha)
    return ((b3y-p1.y)*vp.x - (b3x-p1.x)*vp.y)/(m3x*vp.y - m3y*vp.x)
end

"""
eqs1(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual

Explicit equation for when parameter `s` of segment `q1q2` satisfies `s=1` 
"""
function eqs1(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual
    b3x,b3y,m3x,m3y = bsms(q1,q2,vq,alpha)
    return (b3y - p1.y + m3y)*vp.x - (b3x - p1.x + m3x)*vp.y 
end

"""
eqt0(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual

Explicit equation for when parameter `t` of segment `p1p2` satisfies `t=0` 
"""
function eqt0(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual
    s,t = st(q1,q2,vq,p1,p2,vp,alpha)
    return t
end

"""
eqt1(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual

Explicit equation for when parameter `t` of segment `p1p2` satisfies `t=1` 
"""
function eqt1(q1::Point,q2::Point,vq::Point,p1::Point,p2::Point,vp::Point,alpha)::RealOrDual
    s,t = st(q1,q2,vq,p1,p2,vp,alpha)
    return t-1  
end


"""
rootsAldo(f,epsilon)

Function to obtain roots of function `f` up to a tolerance `epsilon` using Newton's method with derivative calculated using Hyperdual numbers
"""
function rootsAldo(f,epsilon)
    df(x) = hyperdualTotalDerivative(f,x,1e-8)
    root = newton(f,epsilon;x0=rand(),df=df)
    return root
end

# the following are auxiliary functions that are used for verifying that an angle is a point where an intersection is generated

function boolBasic(val1::RealOrDual,val2::RealOrDual,compar::RealOrDual)::Bool
    return xor(val1 > compar && val2 < compar, val1 < compar && val2 > compar)
end

function boolFuncs0(s1::RealOrDual,t1::RealOrDual,s2::RealOrDual,t2::RealOrDual)::Bool
    return boolBasic(s1,s2,0)
end

function boolFuncs1(s1::RealOrDual,t1::RealOrDual,s2::RealOrDual,t2::RealOrDual)::Bool
    return boolBasic(s1,s2,1)
end

function boolFunct0(s1::RealOrDual,t1::RealOrDual,s2::RealOrDual,t2::RealOrDual)::Bool
    return boolBasic(t1,t2,0)
end

function boolFunct1(s1::RealOrDual,t1::RealOrDual,s2::RealOrDual,t2::RealOrDual)::Bool
    return boolBasic(t1,t2,1)
end



"""
function checkPossibleEquality(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,

Checks if it is you can rotate the segment `q1q2` so that it becomes part of the same line as segment `p1p2`. It returns a `true,ang` with `ang` the possible angle if its possible and `false,0` if it isn't.
"""
function checkPossibleEquality(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual)
    ang_q = atan(vq.y/vq.x)
    ang_p = atan(vp.y/vp.x)
    ang = ang_p - ang_q
    # test if angle of paralellism is in interval of rotation
    if angleTest(ang,alpha)
        # using the `c` from the  `ax + by + c = 0` representation with (a,b,c) normalized  
        # to test if lines become the same at the angle they are parallel.
        # note that comparing the b from the y = mx + b fails for other cases
        mat = zrotation(ang)
        q1 = mat*q1
        q2 = mat*q2
        vq = q2 - q1
        a1,b1,c1 = vq.x, -vq.y, vq.y*q1.x - vq.x*q1.y
        nor1 = sqrt(pow2(a1) + pow2(b1) + pow2(c1))
        a2,b2,c2 = vp.x, -vp.y, vp.y*p1.x - vp.x*p1.y
        nor2 = sqrt(pow2(a2) + pow2(b2) + pow2(c2))
        return isapprox(c1/nor1,c2/nor2,atol=1e-15), ang
    end
    return false, 0
end


"""
function equalSegmentIntersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool

Checks if it is you can rotate the segment `q1q2` so that it becomes part of the same line as segment `p1p2`. It returns a `true,ang` with `ang` the possible angle if its possible and `false,0` if it isn't.
"""
function equalSegmentIntersection(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point)::Bool
    # to se if they intersect, it is enough to check wheter theyr x or w intervals overlap
    if isapprox(vp.x,0,atol=1e-15)
        a,b = minmax(p1.y,p2.y)
        c,d = minmax(q1.y,q2.y)
        return (b >= c) && (d >= a)
    else
        a,b = minmax(p1.x,p2.x)
        c,d = minmax(q1.x,q2.x)
        return (b >= c) && (d >= a)
    end
end

# using IntervalRootFinding

"""
case3IntersectionLong(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual)::Bool

Checks if XY-contained segment `p1p2` and `q1q2` intersect when rotating the segment `q1q2` an angle `alpha`. 

It checks for the intersection by first checking if lines already intersect. If not, then the intersection must happen when rotating. If they become a part of the same line when rotating, then that is subcase 3.1 and must be treated specially.
    
If they dont become part of the same line when rotating,then if they intersect they will  to intersect when endpoint must intersect with the other point. That is, for angle `ang`, the one of the equations `s=0,1` and `t=0,1` has a root in the interval `0,alpha` or `alpha,0`.

If the equations for `s` have an angle root, we must check that for that angle and for that intersection `t` is in ``0,1``. The same happens for when obtaining angles for the `t` equation: we must check that for that angle the `s` value is in `0,1`

For the rootfinding of the system, the IntervalRootFinding.jl package in order to ensure that the numerical method is validated and that all of the roots are found. There is no need to check for all roots of each equation, it is enough to check for the maximum and minimum. 
"""
function case3IntersectionLong(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual)::Bool
    # we first check if lines already are at an intersection
    # this case should never happen but we add it for caution
    if xySegmentIntersection(p1,p2,vp,q1,q2,vq)
        #println("# already intersected")
        return true
    # if intersection is posible
    elseif firstFilter(p1,p2,vp,q1,q2,vq)
        # check if segments become part of the same line at one point
        #println("# are not already intersected")
        test,ang = checkPossibleEquality(p1,p2,vp,q1,q2,vq,alpha)
        if test
            #println("# they become the same line")
            # they become part of the same line, 
            # intersection must be checked in a different way 
            mat = zrotation(ang)
            q1prime = mat*q1
            q2prime = mat*q2
            return  equalSegmentIntersection(p1,p2,vp,q1prime,q2prime,q2prime-q1prime)
        else
            # they are not part of the same line
            delta = pi/100
            gs = [eqs0,eqs1,eqt0,eqt1]
            boolFuncs = [boolFuncs0,boolFuncs1,boolFunct0,boolFunct1]
            for i in 1:4
                func(ang) = gs[i](q1,q2,vq,p1,p2,vp,ang)
                #alphas[i] = rootsAldo(func,1e-12)
                angleRoots = IntervalRootFinding.roots(func,-pi..pi)
                # converting from interval to float by taking midpoint
                angleRoots = [midpoint_radius(a.interval)[1] for a in angleRoots]
                # fixing choosing extremal angles
                a,b = minimum(angleRoots), maximum(angleRoots)
                for ang in (a,b)
                #for ang in angleRoots
                    if angleTest(ang,alpha)
                        s,t = st(q1,q2,vq,p1,p2,vp,ang)
                        # the first two equations guarantee s=0,1 so we need to check for t
                        # the last two have t=0,1 so we need to check for s
                        secondTest = (i <= 2) ? (0 <= t <= 1) : 0 <= s <= 1
                        if secondTest
                            #push!(angs,ang)
                            #push!(eqs,string(gs[i]))
                            s1,t1 = st(q1,q2,vq,p1,p2,vp,ang+delta)
                            s2,t2 = st(q1,q2,vq,p1,p2,vp,ang-delta)
                            if boolFuncs[i](s1,t1,s2,t2)
                                return true
                            end
                        end
                    end
                end
            end
            #return angs,eqs
            return false
        end
    else
        return false
    end
end

"""
case3IntersectionFast(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual;n::Integer=1000)::Bool

Testing if the the surface obtained by rotating the segment `q1q2` intersects the `p1p2` segment by decomposing the surface as `n+1` line segments and checking if one of them intersects the `p1p2` segment.

Before checking all the intersections, it checks for subcase 3.1: the segments become part of the same line at one point.
"""
function case3IntersectionFast(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual;n::Integer=1000)::Bool
    if firstFilter(p1,p2,vp,q1,q2,vq)
        test,ang = checkPossibleEquality(p1,p2,vp,q1,q2,vq,alpha)
        if test
            # they become part of the same line, 
            # intersection must be checked in a different way 
            mat = zrotation(ang)
            q1prime = mat*q1
            q2prime = mat*q2
            return  equalSegmentIntersection(p1,p2,vp,q1prime,q2prime,q2prime-q1prime)
        else
            flag = true
            alphas = LinRange(0,alpha,n)
            c = 1
            while c <= n && flag
                s,t = st(q1,q2,vq,p1,p2,vp,alphas[c])
                flag = flag && ( t<0 || t>1 || s<0 || s>1) 
                c += 1
            end
            return !flag
        end
    end    
end


"""
bm(p1::Point,p2::Point,vp::Point)

Returns the `b,m` values of the representation of de XY line passing though points `p1p2`
"""
function bm(p1::Point,p2::Point,vp::Point)
    m = vp.y/vp.x
    b = p1.y - m*p1.x
    return b,m
end

#= 

function experiment2(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual)
    if isapprox(vp.x,0,atol=1e-15)
        alphap = pi/2
        alphaq = atan(vq.y/vq.x)
        ang = alphap - alphaq
        # test if angle of change to paralellism is in interval of rotation
        if angleTest(ang,alpha)
            # lines become parallel, time to test for 
            b3x,b3y,m3x,m3y = bsms(q1,q2,vq,ang)
            if isapprox(b3x,p1.x,atol=1e-15)
                a,b = minmax(p1.y,p2.y)
                c,d = minmax(b3y,b3y + m3y)
                return (b >= c) && (d >= a)
            end
        end
        return false
    else
        ang_q = atan(vq.y/vq.x)
        ang_p = atan(vp.y/vp.x)
        ang = ang_p - ang_q
        # test if angle of change of paralellism is in interval of rotatoion
        # lines become parallel at some point we need to test if they are the same line
        if angleTest(ang,alpha)
            # using the `b` from the  `y = mx + b` representation with (a,b,c) normalized  
            # to test if lines become the same at the angle they are parallel.
            b3x,b3y,m3x,m3y = bsms(q1,q2,vq,ang)
            bp,mp = bm(p1,p2,vp)
            mq = m3y/m3x
            bq = b3y - b3x*mq
            # check wheter lines become the same
            if isapprox(bq,bp,atol=1e-15)
                # to se if they intersect, it is enough to check wheter their x intervals overlap
                a,b = minmax(p1.x,p2.x)
                c,d = minmax(b3x,b3x + m3x)
                return (b >= c) && (d >= a)
            end
        end
        return false
    end
end

function experiment3(p1::Point,p2::Point,vp::Point,q1::Point,q2::Point,vq::Point,alpha::RealOrDual)::Bool
    test,ang = checkPossibleEquality(vp,vq,alpha)
    if test
        mat = zrotation(ang)
        q1prime = mat*q1
        q2prime = mat*q2
        # using the `c` from the  `ax + by + c = 0` representation with (a,b,c) normalized  
        # to test if lines become the same at the angle they are parallel.
        # note that comparing the b from the y = mx + b fails for other cases
        
        b3x,b3y,m3x,m3y = bsms(q1,q2,vq,ang)
        a1,b1,c1 = m3x, -m3y, m3y*b3x - m3x*b3y
        nor1 = sqrt(pow2(a1) + pow2(b1) + pow2(c1))
        a2,b2,c2 = vp.x, -vp.y, vp.y*p1.x - vp.x*p1.y
        nor2 = sqrt(pow2(a2) + pow2(b2) + pow2(c2))
        # check wheter lines become the same
        if isapprox(c1/nor1,c2/nor2,atol=1e-15)
            # to se if they intersect, it is enough to check wheter theyr x or w intervals overlap
            if isapprox(vp.x,0,atol=1e-15)
                a,b = minmax(p1.y,p2.y)
                c,d = minmax(b3y,b3y + m3y)
                return (b >= c) && (d >= a)
            else
                a,b = minmax(p1.x,p2.x)
                c,d = minmax(b3x,b3x + m3x)
                return (b >= c) && (d >= a)
            end
        end
        
    end
    return false
end


=#


function surfaceSegmentIntersection(p1::Point,p2::Point,q1::Point,q2::Point,alpha::RealOrDual;debug=false)
    vp = p2-p1
    vq = q2-q1
    A = isapprox(vp.z,0,atol=1e-15)
    B = isapprox(vq.z,0,atol=1e-15)
    # case 1: the segments are not contained in an `xy` plane
    if !B
        debug && println("# case 1")
        return case1Intersection(p1,p2,vp,q1,q2,vq,alpha,debug=debug)
    # case 2: the segment
    elseif !A && B
        debug && println("# case 2")
        return case2Intersection(p1,p2,vp,q1,q2,vq,alpha,debug=debug)
    elseif A && B
        debug && println("# case 3")
        return case3Intersection(p1,p2,vp,q1,q2,vq,alpha,debug=debug)
    end
end

"""
specialIntersection(p1::Point,p2::Point,p3::Point,alpha::RealOrDual;debug=false)

Case for moving internal angles 
"""
function specialIntersection(p1::Point,p2::Point,p3::Point,alpha::RealOrDual;debug=false)
    ang = xyIangle(p1-p2,p3-p2)
    if alpha > 0
        return ang + alpha >= 2*pi
    else
        return -alpha >= ang
    end
end
"""
planeRotationXY(u::Point)::Matrix

returns the matrix representing a rotation of the system so that vector `u` becomes colineal with the z axis
"""
function planeRotationXY(u::Point)::Matrix
    if isapprox(distance(ez,u),0,atol=1e-15)
        return Matrix(1,0,0,0,1,0,0,0,1)
    elseif isapprox(distance(-1*ez,u),0,atol=1e-15)
        return Matrix(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
    else
        ang = iangle(u,ez)
        u = unitVector(cross(u,ez))
        return rotation(ang,u)
    end
end

function planeRotationYZ(u::Point)::Matrix
    if isapprox(distance(ex,u),0,atol=1e-15)
        return Matrix(1,0,0,0,1,0,0,0,1)
    elseif isapprox(distance(-1*ex,u),0,atol=1e-15)
        return Matrix(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0)
    else
        theta = iangle(u,ex)
        u = unitVector(cross(u,ex))
        return rotation(theta,u)
    end
end

function moveBeforeDihedral(P::PolygonalChain,i::Integer)::PolygonalChain
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

function moveBeforeInternal(P::PolygonalChain,i::Integer)
    newP = copy(P)
    n = length(P)
    if 1 <= i <= n-1
        u = unitVector(newP[i]-newP[i+1])
        mat1 = planeRotationYZ(u)
        pivot = newP[i+1]
        for j in 1:n+1
            newP[j] = mat1*(newP[j]-pivot)
        end
        # second rotation
        v = unitVector(newP[i+2])  # newP[i+1] = e0
        n1 = cross(v,ex)
        theta = atan(newP[i+2].y,newP[i+2].z)-pi/2
        mat2 = xrotation(theta)
        for j in 1:n+1
            newP[j] = mat2*newP[j]
        end
        return newP
    else
        error("i is not in range of internals")
    end
end

function checkRotationIntersection(P::AbstractChain,k::Integer,alpha::RealOrDual,internal::Bool;debug::Bool=false)::Bool
    n = length(P)
    flag = false
    # testing first the intersection between the adjacent edges
    # for internal angle modifications only
    if internal
        flag = specialIntersection(P[k],P[k+1],P[k+2],alpha,debug=debug)
        if debug && flag
            println("==========================")
            println("# INTERSECTION (internal special case): ")
            println("k = $k")
            println("theta = $(alpha)")
            println()
            println("=========================")
        end
    end
    debug && println("k = $k")
    # index of segments that will rotate
    i = internal ? k+1 : k+2
    debug && println("i = $i")
    while i <= n && !flag
        # indexes of static semgments
        if internal
            suplim = i==k+1 ? k-1 : k
        else
            suplim = i==k+2 ? k : k+1
        end
        j = 1
        while j <= suplim && !flag
            debug && println("i=$i, j=$j")
            vpj = P[j+1]-P[j]
            vpi = P[i+1]-P[i]
            inter = surfaceSegmentIntersection(P[j],P[j+1],P[i],P[i+1],alpha,debug=debug)[1]
            if debug && inter
                println("===================================")
                println(" #### INTERSECTION #### ")
                println("k = $(k)")
                println("i = $(i)")
                println("j = $(j)")
                println("alpha = $(alpha)")
                println()
                println("===================================")
            end
            flag = flag || inter
            j += 1
        end
        i += 1
    end
    return flag
end


if abspath(PROGRAM_FILE) == @__FILE__

end