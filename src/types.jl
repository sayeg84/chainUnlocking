using LinearAlgebra, BenchmarkTools, StaticArrays

const T = BigFloat

abstract type AbstractPoint end

struct Point <: AbstractPoint
    x::T
    y::T
    z::T
end

const ex = Point(1.0,0.0,0.0)
const ey = Point(0.0,1.0,0.0)
const ez = Point(0.0,0.0,1.0)
const e0 = Point(0.0,0.0,0.0)

function pow2(x::Number)::Number
    return x*x
end

function minmax(x::Real,y::Real)
    return x < y ? (x,y) : (y,x)
end

function Point()
    return Point(rand(),rand(),rand())
end

#import Base.:+, Base.:-, Base.:*, Base.transpose,  Base.copy, Base.rand

function Base.:+(p::Point,q::Point)::Point
    return Point(p.x + q.x, p.y + q.y, p.z + q.z)
end

function Base.:*(a::Real,p::Point)::Point
    return Point(a*p.x, a*p.y, a*p.z)
end

function Base.:*(p::Point,a::Real)::Point
    return a*p
end

function Base.:/(p::Point,a::Real)::Point
    return Point(p.x/a, p.y/a, p.z/a)
end

function Base.:-(p::Point,q::Point)::Point
    return Point(p.x - q.x, p.y - q.y, p.z - q.z)
end

function Base.copy(p::Point)::Point
    return Point(p.x, p.y, p.z)
end
 
"""
dot(p::Point,q::Point)::T

Dot product for points
"""
function dot(p::Point,q::Point)::T
    return p.x*q.x + p.y*q.y + p.z*q.z
end

"""
cross(p::Point,q::Point)::Point

Cross product for two Points
"""
function cross(p::Point,q::Point)::Point
    return Point(p.y*q.z - p.z*q.y,  p.z*q.x - p.x*q.z,  p.x*q.y - p.y*q.x)
end

import LinearAlgebra.norm

function norm(p::Point)::T
    return sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z))
end

function snorm(q::Point;s::T=2.0)::T
    return (q.x^s + q.y^s + q.z^s)^(1/s)
end

function distance(p::Point,q::Point)::T
    return norm(p-q)
end

function sdistance(p::Point,q::Point;s::T=2.0)::T
    return snorm(p-q,s=s)
end

function unitVector(p::Point)::Point
    return p/norm(p)
end

"""
iangle(u::Point,v::Point)::T 

returns the internal angle between vectors `u` and `v`.
"""
function iangle(u::Point,v::Point)::T
    sint = norm(cross(u,v))
    cost = dot(u,v)
    return atan(sint,cost)    
end

function toArray(p::Point)::Array{T,1}
    return [p.x,p.y,p.z]
end

struct Matrix 
    xx::T
    xy::T
    xz::T
    yx::T
    yy::T
    yz::T
    zx::T
    zy::T
    zz::T
end


function Matrix(A::Array{<:Real,2}) 
    return Matrix(A[1,1],A[1,2],A[1,3],A[2,1],A[2,2],A[2,3],A[3,1],A[3,2],A[3,3])
end

function Matrix()
    return Matrix(rand(),rand(),rand(),rand(),rand(),rand(),rand(),rand(),rand())
end

function matrixFromRows(row1::Array{<:Real,1},row2::Array{<:Real,1},row3::Array{<:Real,1})
    if length(row1)!=3
        error("row 1 is not of length 3")
    elseif length(row2)!=3
        error("row 2 is not of length 3")
    elseif length(row3)!=3
        error("row 3 is not of length 3")
    end
    return Matrix(row1[1],row1[2],row1[3],row2[1],row2[2],row2[3],row3[1],row3[2],row3[3])
end

function matrixFromCols(col1::Array{<:Real,1},col2::Array{<:Real,1},col3::Array{<:Real,1})
    if length(col1)!=3
        error("col 1 is not of length 3")
    elseif length(col2)!=3
        error("col 2 is not of length 3")
    elseif length(col3)!=3
        error("col 3 is not of length 3")
    end
    return Matrix(col1[1],col2[1],col3[1],col1[2],col2[2],col3[2],col1[3],col2[3],col3[3])
end

function matrixFromRows(row1::Point,row2::Point,row3::Point)
    return Matrix(row1.x,row1.y,row1.z,row2.x,row2.y,row2.z,row3.x,row3.y,row3.z)
end

function matrixFromCols(col1::Point,col2::Point,col3::Point)
    return Matrix(col1.x,col2.x,col3.x,col1.y,col2.y,col3.y,col1.z,col2.z,col3.z)
end


function Base.:+(A::Matrix,B::Matrix)::Matrix
    return Matrix(A.xx + B.xx , A.xy + B.xy , A.xz + B.xz, 
    A.yx + B.yx , A.yy + B.yy , A.yz + B.yz, 
    A.zx + B.zx , A.zy + B.zy , A.zz + B.zz)
end

function Base.:*(a::Real,A::Matrix)::Matrix
    return Matrix(a*A.xx , a*A.xy , a*A.xz, 
    a*A.yx , a*A.yy , a*A.yz, 
    a*A.zx , a*A.zy , a*A.zz)
end

function Base.:-(A::Matrix,B::Matrix)::Matrix
    return Matrix(A.xx - B.xx , A.xy - B.xy , A.xz - B.xz, 
    A.yx - B.yx , A.yy - B.yy , A.yz - B.yz, 
    A.zx - B.zx , A.zy - B.zy , A.zz - B.zz)
end

function Base.:*(A::Matrix,B::Matrix)::Matrix
    return Matrix(A.xx*B.xx + A.xy*B.yx + A.xz*B.zx , A.xx*B.xy + A.xy*B.yy + A.xz*B.zy , A.xx*B.xz + A.xy*B.yz + A.xz*B.zz , 
    A.yx*B.xx + A.yy*B.yx + A.yz*B.zx , A.yx*B.xy + A.yy*B.yy + A.yz*B.zy , A.yx*B.xz + A.yy*B.yz + A.yz*B.zz ,
    A.zx*B.xx + A.zy*B.yx + A.zz*B.zx , A.zx*B.xy + A.zy*B.yy + A.zz*B.zy , A.zx*B.xz + A.zy*B.yz + A.zz*B.zz )
end

function Base.:*(A::Matrix,p::Point)::Point
    return Point(A.xx*p.x + A.xy*p.y + A.xz*p.z , A.yx*p.x + A.yy*p.y + A.yz*p.z , A.zx*p.x + A.zy*p.y + A.zz*p.z)
end

function Base.:*(p::Point,A::Matrix)::Point
    return Point(p.x*A.xx + p.y*A.yx + p.z*A.zx , p.x*A.xy + p.y*A.yy + p.z*A.zy , p.x*A.xz + p.y*A.yz + p.z*A.zz)
end

function Base.copy(A::Matrix)::Matrix
    return  Matrix(A.xx, A.xy, A.xz, 
    A.yx, A.yy, A.yz,
    A.zx, A.zy, A.zz)
end

function Base.transpose(A::Matrix)::Matrix
    return  Matrix(A.xx, A.yx, A.zx,
    A.xy, A.yy, A.zy,
    A.xz, A.yz, A.zz)
end

"""
dot(A::Matrix,B::Matrix)::T

Dot product for matrices
"""
function dot(A::Matrix,B::Matrix)::T
    return A.xx*B.xx + A.xy*B.xy + A.xz*B.xz + A.yx*B.yx + A.yy*B.yy + A.yz*B.yz + A.zx*B.zx + A.zy*B.zy + A.zz*B.zz
end

function xrotation(ang::Real)::Matrix
    a = cos(ang)
    b = sin(ang)
    return Matrix(1, 0, 0,
    0, a, -b,
    0, b, a)
end

function yrotation(ang::Real)::Matrix
    a = cos(ang)
    b = sin(ang)
    return Matrix(a, 0, b,
    0, 1, 0,
    -b, 0, a)
end

function zrotation(ang::Real)::Matrix
    a = cos(ang)
    b = sin(ang)
    return Matrix(a, -b, 0,
    b, a, 0,
    0, 0, 1)
end


"""
rotation(ang::T,u::Point)::Matrix

Computes the rotation matrix around unit vector `u` for an angle `ang` using the explicit formula in https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
"""
function rotation(ang::Real,u::Point)::Matrix
    a = cos(ang)
    b = sin(ang)
    c = 1-a
    return Matrix(a + pow2(u.x)*c, u.x*u.y*c - u.z*b, u.x*u.z*c + u.y*b,
    u.y*u.x*c + u.z*b, a + pow2(u.y)*c, u.y*u.z*c - u.x*b,
    u.z*u.x*c - u.y*b, u.z*u.y*c + u.x*b, a + pow2(u.z)*c)
end


"""
rotate(p::Point,ang::T,u::Point)::Point

Rotate a point `p` an angle `ang` using the plane defined by unit vector `u` by multiplying point `p` by a rotation matrix.
"""
function rotate(p::Point,ang::Real,u::Point)::Point
    return rotation(ang,u)*p
end

"""
rotateRodrigues(p::Point,ang::T,u::Point)::Point

Rotate a point `p` an angle `ang` using the plane defined by unit vector `u` by using Rodrigues rotation forumla
"""
function rotateRodrigues(p::Point,ang::Real,u::Point)::Point
    a = cos(ang)
    b = sin(ang)
    return a*p + b*cross(u,p) + (1-a)*dot(u,p)*u
end

"""
bangle(a::Point,b::Point,c::Point)::T

Calculates the bond angle between a link defined by points b-a, c-b. 
The angle is the angle between the links (Colinear points will return pi)
"""
function bangle(a::Point,b::Point,c::Point)::T
    u1 = b-a
    u2 = c-b
    ang = dot(u1,u2)/(norm(u1)*norm(u2))
    return pi - acos(ang)
end


"""
dihedral(a::Point,b::Point,c::Point,d::Point)::T

Calculates dihedral angle between the planes defined by points `abc` and `bcd` using the formula in the [Wikipedia article](https://en.wikipedia.org/wiki/Dihedral_angle)
"""
function dihedral(a::Point,b::Point,c::Point,d::Point)::T
    u1 = b-a
    u2 = c-b
    u3 = d-c
    v = cross(u2,u3)
    y = norm(u2)*dot(u1,v)
    x = dot(cross(u1,u2),v)
    return atan(y,x)
end

"""
dihedralSO(a::Point,b::Point,c::Point,d::Point)::T

Calculates dihedral angle between the planes defined by points `abc` and `bcd` using the formula in this [Math SE question](https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates)
"""
function dihedralSO(a::Point,b::Point,c::Point,d::Point)::T
    b1 = b-a
    b2 = c-b
    b3 = d-c
    n1 = unitVector(cross(b1,b2))
    n2 = unitVector(cross(b2,b3))
    m1 = cross(n1,unitVector(b2))
    x = dot(n1,n2)
    y = dot(m1,n2)
    return -atan(y,x)
end

function distance(p::Tuple{<:Real,<:Real},q::Tuple{<:Real,<:Real})
    return sqrt(pow2(p[1]-q[1]) + pow2(p[2]-q[2]))
end

function distance(p::Tuple{<:Real,<:Real,<:Real},q::Tuple{<:Real,<:Real,<:Real})
    return sqrt(pow2(p[1]-q[1]) + pow2(p[2]-q[2]) + pow2(p[3]-q[3]))
end

abstract type AbstractChain end

struct PolygonalOld{N} <: AbstractChain
    vertices::NTuple{N,Point}
end

function PolygonalOld(A::Array{Point,1})
    return PolygonalOld(Tuple(A))
end

function PolygonalOld(n::Integer)
    A = [Point() for i in 1:n]
    return PolygonalOld(Tuple(A))
end

struct PolygonalChain <: AbstractChain
    vertices::Array{Point,1}
end

function PolygonalChain(n::Integer)
    A = [Point() for i in 1:n]
    return PolygonalChain(A)
end

function PolygonalChain(arr::Array{<:Real,2})
    if size(arr)[2] == 3
        A = [Point(arr[i,1],arr[i,2],arr[i,3]) for i in 1:(size(arr)[1])]
        return PolygonalChain(A)
    else
        error("Array dimensions must be (n,3)")
    end
end

import Base.length, Base.copy, Base.getindex, Base.setindex!, Base.lastindex

function Base.length(P::AbstractChain)::Int64
    return length(P.vertices)-1
end

function Base.copy(P::AbstractChain)
    return typeof(P)(copy(P.vertices))
end

function Base.getindex(P::AbstractChain,idx::Integer)
    return P.vertices[idx]
end

function Base.setindex!(P::AbstractChain,p::AbstractPoint,idx::Integer)
    P.vertices[idx] = p
end

function Base.lastindex(P::AbstractChain)
    return length(P.vertices)
end

function toArray(P::AbstractChain)::Array{T,2}
    return vcat([[p.x p.y p.z] for p in P.vertices]...)
end

function centroid(P::AbstractChain)::Point
    n = length(P)
    c = e0
    for i in 1:n+1
        c += P[i]
    end
    return 1/(n+1)*c
end

function centerChain(P::PolygonalOld)::PolygonalOld
    c = centroid(P)
    arr = [p-c for p in P.vertices]
    return PolygonalOld(arr)
end

function centerChain(P::PolygonalChain)::PolygonalChain
    c = centroid(P)
    arr = [p-c for p in P.vertices]
    return PolygonalChain(arr)
end


function centerChain!(P::PolygonalChain)
    c = centroid(P)
    n = length(P)
    for i in 1:(n+1)
        P[i] -= c
    end
end

using GenericLinearAlgebra

"""
optimalRotation(P::AbstractChain,Q::AbstractChain)::Matrix

Taking as an input two centered (i.e. with centroid 0) AbstractChains P,Q, it uses the [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm) to find the optimal rotation to overlap P into q
"""

function optimalRotation(P::AbstractChain,Q::AbstractChain)::Matrix
    A = toArray(P)
    B = toArray(Q)
    H = transpose(A)*B
    #println("Aca")
    #=
    if T == BigFloat
        display(A)
        println()
        display(B)
        println()
        println("Big")
        display(H)
        println()
        S = GenericLinearAlgebra.svd(H,full=false)    
        println("No salgo")
    else
        S = LinearAlgebra.svd(H)    
    end
    =#
    S = LinearAlgebra.svd(Float64.(H))
    d = sign(LinearAlgebra.det(transpose(S.U*S.Vt)))
    mat = transpose(S.Vt)*[1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 d ]*transpose(S.U)
    return Matrix(mat)
end

"""
simpleRmsd(P::AbstractChain,Q::AbstractChain)::T

Calculates RMSD for two AbstractChains. 
"""
function simpleRmsd(P::AbstractChain,Q::AbstractChain)::T
    n = length(P)
    if n != length(Q)
        error("Chains must be of same length")
    end
    dist = 0
    for i in 1:n+1
        dist += pow2(distance(P[i],Q[i]))
    end
    return sqrt(dist/(n+1))
end

"""
overlapedRmsd(P::AbstractChain,Q::AbstractChain)::T

Calculates RMSD for two AbstractChains by centering then and rotating them for maximum overlap
"""
function overlapedRmsd(P::AbstractChain,Q::AbstractChain)::T
    P = centerChain(P)
    Q = centerChain(Q)
    mat = optimalRotation(P,Q)
    arr = [mat*p for p in P.vertices]
    S = typeof(P)(arr)
    n = length(S)
    if n != length(Q)
        error("Chains must be of same length")
    end
    dist = 0
    for i in 1:n+1
        dist += pow2(distance(S[i],Q[i]))
    end
    return sqrt(dist/(n+1))
end


## Geometrical functions

function linkLengths(P::AbstractChain)::Array{T,1}
    n = length(P)
    lengths = zeros(T,n)
    for i in 1:n
        lengths[i] = distance(P[i],P[i+1])
    end
    return lengths
end

function linkAngles(P::AbstractChain)::Array{T,1}
    n = length(P)
    angles = zeros(T,n-1)
    for i in 1:n-1
        angles[i] = bangle(P[i],P[i+1],P[i+2])
    end
    return angles
end

function dihedralAngles(P::AbstractChain)::Array{T,1}
    n = length(P)
    dihedrals = zeros(T,n-2)
    for i in 1:n-2
        dihedrals[i] = dihedral(P[i],P[i+1],P[i+2],P[i+3])
    end
    return dihedrals
end

function lengthsAndAngles(P::AbstractChain)::Tuple{Array{T,1},Array{T,1},Array{T,1}}
    n = length(P)
    lengths = zeros(T,n)
    angles = zeros(T,n-1)
    dihedrals = zeros(T,n-2)
    lengths[1] = distance(P[1],P[2])
    lengths[2] = distance(P[2],P[3])
    angles[1] = bangle(P[1],P[2],P[3])
    for i in 3:n
        lengths[i] = distance(P[i],P[i+1])
        angles[i-1] = bangle(P[i-1],P[i],P[i+1])
        dihedrals[i-2] = dihedral(P[i-2],P[i-1],P[i],P[i+1])
    end
    return (lengths,angles,dihedrals)
end

function PolygonalChain(linkLengths::Array{<:Real,1},linkAngles::Array{<:Real,1},dihedralAngles::Array{<:Real,1})::PolygonalChain
    n = length(linkLengths)
    if n == (length(linkAngles)+1) && n == (length(dihedralAngles)+2)
        vertices = Array{Point,1}(undef,n+1)
        # origin as arbitrary starting point
        vertices[1] = e0
        vertices[2] = Point(linkLengths[1],0.0,0.0)
        vertices[3] = vertices[2] + Point(-linkLengths[2]*cos(linkAngles[1]),linkLengths[2]*sin(linkAngles[1]),0.0)
        for i in 3:n
            #ab = vertices[i-1]-vertices[i-2]
            bc = vertices[i]-vertices[i-1]
            #bcnorm = unitVector(bc)
            bcnorm = bc/linkLengths[i-1]
            d0 = linkLengths[i]*bcnorm
            n = unitVector(cross(vertices[i-1]-vertices[i-2],bcnorm))
            d0 = rotate(d0,pi-linkAngles[i-1],n)
            d0 = rotate(d0,dihedralAngles[i-2],bcnorm)
            vertices[i+1] = vertices[i] + d0
        end
        return PolygonalChain(vertices)
    else
        error("Arrays size $n , $(length(linkAngles)+1) and $(length(linkAngles)+2) do not match")
    end
end

function PolygonalChainRosetta(linkLengths::Array{<:Real,1},linkAngles::Array{<:Real,1},dihedralAngles::Array{<:Real,1})::PolygonalChain
    n = length(linkLengths)
    if n == (length(linkAngles)+1) && n == (length(dihedralAngles)+2)
        vertices = Array{Point,1}(undef,n+1)
        # origin as arbitrary starting point
        vertices[1] = e0
        vertices[2] = Point(linkLengths[1],0.0,0.0)
        vertices[3] = vertices[2] + Point(-linkLengths[2]*cos(linkAngles[1]),linkLengths[2]*sin(linkAngles[1]),0.0)
        for i in 3:n
            a = sin(pi-linkAngles[i-1])
            d0 = linkLengths[i]*Point(cos(pi-linkAngles[i-1]),cos(dihedralAngles[i-2])*a,sin(dihedralAngles[i-2])*a)
            bcnorm  = (vertices[i]-vertices[i-1])/linkLengths[i-1]
            n = unitVector(cross(vertices[i-1]-vertices[i-2],bcnorm))
            M = matrixFromCols(bcnorm,cross(n,bcnorm),n)
            #println(M*d0)
            vertices[i+1] = vertices[i] + M*d0
            #display(vertices)
            #println()
        end
        return PolygonalChain(vertices)
    else
        error("Arrays size $n , $(length(linkAngles)+1) and $(length(linkAngles)+2) do not match")
    end
end



function dihedralRotate(P::PolygonalChain,i::Integer,phi::Real)::PolygonalChain
    n = length(P)
    if 1 <= i <= n-2
        rotN = unitVector(P[i+2]-P[i+1])
        rotMat = rotation(phi,rotN)
        points = [p for p in P.vertices]
        for j in i+3:n+1
            points[j] = rotMat*(points[j]-P[i+2]) + P[i+2]
        end
        return PolygonalChain(points)
    else
        error("$i must be lower than $n")
    end
end

"""
Same as `dihedralRotate` but for case where chain `P` has already passed
by `moveBeforeDihedral`
"""
function dihedralRotateFast(P::PolygonalChain,i::Integer,phi::Real)::PolygonalChain
    n = length(P)
    if 1 <= i <= n-2
        rotMat = xrotation(phi)
        points = [p for p in P.vertices]
        for j in i+3:n+1
            points[j] = rotMat*(points[j]-P[i+2]) + P[i+2]
        end
        return PolygonalChain(points)
    else
        error("$i must be lower than $n")
    end
end

"""
Rotating chain inpalce
"""
function dihedralRotate!(P::PolygonalChain,i::Integer,phi::Real)
    n = length(P)
    if 1 <= i <= n-2
        rotN = unitVector(P[i+2]-P[i+1])
        rotMat = rotation(phi,rotN)
        for j in i+3:n+1
            P[j] = rotMat*(P[j]-P[i+1]) + P[i+1]
        end
    else
        error("$i must be lower than $n")
    end
end

function dihedralRotateFast!(P::PolygonalChain,i::Integer,phi::Real)
    n = length(P)
    if 1 <= i <= n-2
        rotMat = xrotation(phi)
        for j in i+3:n+1
            P[j] = rotMat*(P[j]-P[i+1]) + P[i+1]
        end
    else
        error("$i must be lower than $n")
    end
end

function internalRotate(P::PolygonalChain,i::Integer,theta::Real)::PolygonalChain
    n = length(P)
    if 1 <= i <= n-1
        rotN = unitVector(cross(P[i]-P[i+1],P[i+2]-P[i+1]))
        rotMat = rotation(theta,rotN)
        points = [p for p in P.vertices]
        for j in i+2:n+1
            points[j] = rotMat*(points[j]-P[i+1]) + P[i+1]
        end
        return PolygonalChain(points)
    else
        error("$i must be lower than $n")
    end
end

"""
Same as `internalRotate` but for case where chain `P` has already passed
by `moveBeforeInternal`
"""
function internalRotateFast(P::PolygonalChain,i::Integer,theta::Real)::PolygonalChain
    n = length(P)
    if 1 <= i <= n-1
        rotMat = xrotation(theta)
        points = [p for p in P.vertices]
        for j in i+2:n+1
            points[j] = rotMat*(points[j]-P[i+1]) + P[i+1]
        end
        return PolygonalChain(points)
    else
        error("$i must be lower than $n")
    end
end



"""
Rotating chain inpalce
"""
function internalRotate!(P::PolygonalChain,i::Integer,theta::Real)
    n = length(P)
    if 1 <= i <= n-1
        rotN = unitVector(cross(P[i]-P[i+1],P[i+2]-P[i+1]))
        rotMat = rotation(theta,rotN)
        for j in i+2:n+1
            P[j] = rotMat*(P[j]-P[i+1]) + P[i+1]
        end
    else
        error("$i must be lower than $n")
    end
end

function internalRotateFast!(P::PolygonalChain,i::Integer,theta::Real)
    n = length(P)
    if 1 <= i <= n-1
        rotMat = xrotation(theta)
        for j in i+2:n+1
            P[j] = rotMat*(P[j]-P[i+1]) + P[i+1]
        end
    else
        error("$i must be lower than $n")
    end
end

function polygonal(n)
    lengths = rand(T,n)
    angles = pi*rand(T,n-1)
    dihed = pi*([2*rand(T)-1 for i in 1:n-2])
    P = PolygonalOld(lengths,angles,dihed)
    return P
end

function polygonal2(n)
    lengths = rand(T,n)
    angles = pi*rand(T,n-1)
    dihed = pi*([2*rand(T)-1 for i in 1:n-2])
    P = PolygonalChainRosetta(lengths,angles,dihed)
    return P
end

if abspath(PROGRAM_FILE) == @__FILE__
    println()
    println("entramos")
    println()  
    const m = 300
    #display(@benchmark polygonal(m))
    #println()
    #display(@benchmark polygonal2(m))
    #println()
    P = PolygonalOld(m)
    display(@benchmark dihedralRotate(P,10,pi/2))
    println()
    Q = PolygonalChain(m)    
    display(@benchmark dihedralRotate!(Q,10,pi/2))
    println()
    P = PolygonalOld(m)
    Q = PolygonalChain([p for p in P.vertices])
    println(simpleRmsd(P,Q))
    P = dihedralRotate(P,10,pi/2)
    dihedralRotate!(Q,10,pi/2)
    println(simpleRmsd(P,Q))
    a = e0
    b = ex
    c = Point(1.0,1.0,0.0)
    d = Point(1.0,1.0,1.0)
    e = Point(0.0,1.0,1.0)
    P = PolygonalChain([a,b,c,d,e])
    dihedralRotate!(P,2,-pi/4)
    println(P.vertices)
    a = rand(6,3)
    b = PolygonalChain(a)
end