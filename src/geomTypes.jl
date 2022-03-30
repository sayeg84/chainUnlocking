using LinearAlgebra, BenchmarkTools, StaticArrays

include("derivatives.jl")

abstract type AbstractPoint end

struct Point{T<:FloatOrDual} <: AbstractPoint
    x::T
    y::T
    z::T
end

Point(a,b,c)                = Point{typeof(a)}(a,b,c)
Point(U::DataType,p::Point) = Point{U}(p.x,p.y,p.z)
Point()                     = Point{Float64}(rand(),rand(),rand())
Point(arr)                  = Point{typeof(arr[1])}(arr[1],arr[2],arr[3])
Point{Hyperdual}(p::Point{<:AbstractFloat}) = Point{Hyperdual}(Hyperdual(p.x),Hyperdual(p.y),Hyperdual(p.z))

const ex_l = Point{Float64}(1.0,0.0,0.0)
const ey_l = Point{Float64}(0.0,1.0,0.0)
const ez_l = Point{Float64}(0.0,0.0,1.0)
const e0_l = Point{Float64}(0.0,0.0,0.0)

const ex = Point{BigFloat}(1.0,0.0,0.0)
const ey = Point{BigFloat}(0.0,1.0,0.0)
const ez = Point{BigFloat}(0.0,0.0,1.0)
const e0 = Point{BigFloat}(0.0,0.0,0.0)

function pow2(x::Number)
    return x*x
end

function pow2(x::Hyperdual)::Hyperdual
    return x*x
end

function minmax(x::RealOrDual,y::RealOrDual)
    return x < y ? (x,y) : (y,x)
end

#import Base.:+, Base.:-, Base.:*, Base.transpose,  Base.copy, Base.rand

function Base.:+(p::Point,q::Point)::Point
    return Point(p.x + q.x, p.y + q.y, p.z + q.z)
end

function Base.:*(a::RealOrDual,p::Point)::Point
    return Point(a*p.x, a*p.y, a*p.z)
end

function Base.:*(p::Point,a::RealOrDual)::Point
    return a*p
end

function Base.:/(p::Point,a::RealOrDual)::Point
    return Point(p.x/a, p.y/a, p.z/a)
end

function Base.:-(p::Point,q::Point)::Point
    return Point(p.x - q.x, p.y - q.y, p.z - q.z)
end

function Base.copy(p::Point)::Point
    return Point(p.x, p.y, p.z)
end
 
"""
dot(p::Point,q::Point)

Dot product for points
"""
function dot(p::Point,q::Point)
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

function norm(p::Point)
    return sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z))
end

function snorm(q::Point;s::RealOrDual=2.0)
    return (q.x^s + q.y^s + q.z^s)^(1/s)
end

function distance(p::Point,q::Point)
    return norm(p-q)
end

function sdistance(p::Point,q::Point;s::RealOrDual=2.0)
    return snorm(p-q,s=s)
end

function unitVector(p::Point)::Point
    return p/norm(p)
end

"""
iangle(u::Point,v::Point) 

returns the internal angle between vectors `u` and `v`.
"""
function iangle(u::Point,v::Point)
    sint = norm(cross(u,v))
    cost = dot(u,v)
    return atan(sint,cost)    
end

function toArray(p::Point)::Array{<:RealOrDual,1}
    return [p.x,p.y,p.z]
end

struct Matrix{T<:RealOrDual}
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

Matrix(xx,xy,xz,yx,yy,yz,zx,zy,zz) = Matrix{typeof(xx)}(xx,xy,xz,yx,yy,yz,zx,zy,zz)
Matrix(U::DataType,M::Matrix)      = Matrix{U}(M.xx,M.xy,M.xz,M.yx,M.yy,M.yz,M.zx,M.zy,M.zz)

function Matrix(A::Array{<:RealOrDual,2}) 
    return Matrix(A[1,1],A[1,2],A[1,3],A[2,1],A[2,2],A[2,3],A[3,1],A[3,2],A[3,3])
end

function Matrix()
    return Matrix{Float64}(rand(),rand(),rand(),rand(),rand(),rand(),rand(),rand(),rand())
end

function matrixFromRows(row1::Array{<:RealOrDual,1},row2::Array{<:RealOrDual,1},row3::Array{<:RealOrDual,1})
    if length(row1)!=3
        error("row 1 is not of length 3")
    elseif length(row2)!=3
        error("row 2 is not of length 3")
    elseif length(row3)!=3
        error("row 3 is not of length 3")
    end
    return Matrix(row1[1],row1[2],row1[3],row2[1],row2[2],row2[3],row3[1],row3[2],row3[3])
end

function matrixFromCols(col1::Array{<:RealOrDual,1},col2::Array{<:RealOrDual,1},col3::Array{<:RealOrDual,1})
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

function Base.:*(a::RealOrDual,A::Matrix)::Matrix
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
dot(A::Matrix,B::Matrix)

Dot product for matrices
"""
function dot(A::Matrix,B::Matrix)
    return A.xx*B.xx + A.xy*B.xy + A.xz*B.xz + A.yx*B.yx + A.yy*B.yy + A.yz*B.yz + A.zx*B.zx + A.zy*B.zy + A.zz*B.zz
end

function xrotation(ang::RealOrDual)::Matrix
    a = cos(ang)
    b = sin(ang)
    return Matrix{typeof(a)}(1, 0, 0,
    0, a, -b,
    0, b, a)
end

function yrotation(ang::RealOrDual)::Matrix
    a = cos(ang)
    b = sin(ang)
    return Matrix{typeof(a)}(a, 0, b,
    0, 1, 0,
    -b, 0, a)
end

function zrotation(ang::RealOrDual)::Matrix
    a = cos(ang)
    b = sin(ang)
    return Matrix{typeof(a)}(a, -b, 0,
    b, a, 0,
    0, 0, 1)
end


"""
rotation(ang::RealOrDual,u::Point)::Matrix

Computes the rotation matrix around unit vector `u` for an angle `ang` using the explicit formula in https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
"""
function rotation(ang::RealOrDual,u::Point)::Matrix
    a = cos(ang)
    b = sin(ang)
    c = 1-a
    return Matrix(a + pow2(u.x)*c, u.x*u.y*c - u.z*b, u.x*u.z*c + u.y*b,
    u.y*u.x*c + u.z*b, a + pow2(u.y)*c, u.y*u.z*c - u.x*b,
    u.z*u.x*c - u.y*b, u.z*u.y*c + u.x*b, a + pow2(u.z)*c)
end


"""
rotate(p::Point,ang::RealOrDual,u::Point)::Point

Rotate a point `p` an angle `ang` using the plane defined by unit vector `u` by multiplying point `p` by a rotation matrix.
"""
function rotate(p::Point,ang::RealOrDual,u::Point)::Point
    return rotation(ang,u)*p
end

"""
rotateRodrigues(p::Point,ang::RealOrDual,u::Point)::Point

Rotate a point `p` an angle `ang` using the plane defined by unit vector `u` by using Rodrigues rotation forumla
"""
function rotateRodrigues(p::Point,ang::RealOrDual,u::Point)::Point
    a = cos(ang)
    b = sin(ang)
    return a*p + b*cross(u,p) + (1-a)*dot(u,p)*u
end

"""
bangle(a::Point,b::Point,c::Point)

Calculates the bond angle between a link defined by points b-a, c-b. 
The angle is the angle between the links (Colinear points will return pi)
"""
function bangle(a::Point,b::Point,c::Point)
    u1 = b-a
    u2 = c-b
    ang = dot(u1,u2)/(norm(u1)*norm(u2))
    return pi - acos(ang)
end


"""
dihedral(a::Point,b::Point,c::Point,d::Point)

Calculates dihedral angle between the planes defined by points `abc` and `bcd` using the formula in the [Wikipedia article](https://en.wikipedia.org/wiki/Dihedral_angle)
"""
function dihedral(a::Point,b::Point,c::Point,d::Point)
    u1 = b-a
    u2 = c-b
    u3 = d-c
    v = cross(u2,u3)
    y = norm(u2)*dot(u1,v)
    x = dot(cross(u1,u2),v)
    return atan(y,x)
end

"""
dihedralSO(a::Point,b::Point,c::Point,d::Point)

Calculates dihedral angle between the planes defined by points `abc` and `bcd` using the formula in this [Math SE question](https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates)
"""
function dihedralSO(a::Point,b::Point,c::Point,d::Point)
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

function distance(p::Tuple{<:RealOrDual,<:RealOrDual},q::Tuple{<:RealOrDual,<:RealOrDual})
    return sqrt(pow2(p[1]-q[1]) + pow2(p[2]-q[2]))
end

function distance(p::Tuple{<:RealOrDual,<:RealOrDual,<:RealOrDual},q::Tuple{<:RealOrDual,<:RealOrDual,<:RealOrDual})
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

function PolygonalChain(arr::Array{<:RealOrDual,1})
    np = length(arr) / 3
    if isinteger(np)
        A = [Point(arr[3*i-2],arr[3*i-1],arr[3*i]) for i in 1:Int(np)]
        return PolygonalChain(A)
    else
        error("Array dimensions must be divisible by 3")
    end
end

function PolygonalChain(arr::Array{<:RealOrDual,2})
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

function to2DArray(P::AbstractChain)
    return vcat([[p.x p.y p.z] for p in P.vertices]...)
end

function toArray(P::AbstractChain)
    arr = to2DArray(P)
    return reshape(transpose(arr),length(arr))
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

#using GenericLinearAlgebra

"""
optimalRotation(P::AbstractChain,Q::AbstractChain)::Matrix

Taking as an input two centered (i.e. with centroid 0) AbstractChains P,Q, it uses the [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm) to find the optimal rotation to overlap P into q
"""

function optimalRotation(P::AbstractChain,Q::AbstractChain)::Matrix
    A = to2DArray(P)
    B = to2DArray(Q)
    H = transpose(A)*B
    S = LinearAlgebra.svd(Float64.(H))
    d = sign(LinearAlgebra.det(transpose(S.U*S.Vt)))
    mat = transpose(S.Vt)*[1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 d ]*transpose(S.U)
    return Matrix(mat)
end

"""
simpleRmsd(P::AbstractChain,Q::AbstractChain)

Calculates RMSD for two AbstractChains. 
"""
function simpleRmsd(P::AbstractChain,Q::AbstractChain)
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
overlapedRmsd(P::AbstractChain,Q::AbstractChain)

Calculates RMSD for two AbstractChains by centering then and rotating them for maximum overlap
"""
function overlapedRmsd(P::AbstractChain,Q::AbstractChain)
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

function linkLengths(P::AbstractChain)
    n = length(P)
    lengths = zeros(typeof(P[1].x),n)
    for i in 1:n
        lengths[i] = distance(P[i],P[i+1])
    end
    return lengths
end

function linkAngles(P::AbstractChain)
    n = length(P)
    angles = zeros(typeof(P[1].x),n-1)
    for i in 1:n-1
        angles[i] = bangle(P[i],P[i+1],P[i+2])
    end
    return angles
end

function dihedralAngles(P::AbstractChain)
    n = length(P)
    dihedrals = zeros(typeof(P[1].x),n-2)
    for i in 1:n-2
        dihedrals[i] = dihedral(P[i],P[i+1],P[i+2],P[i+3])
    end
    return dihedrals
end

function internalCoordinates(P::AbstractChain)
    n = length(P)
    T = typeof(P[1].x)
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

function PolygonalChain(linkLengths::Array{<:RealOrDual,1},linkAngles::Array{<:RealOrDual,1},dihedralAngles::Array{<:RealOrDual,1})::PolygonalChain
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

function PolygonalChainRosetta(linkLengths::Array{<:RealOrDual,1},linkAngles::Array{<:RealOrDual,1},dihedralAngles::Array{<:RealOrDual,1})::PolygonalChain
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

function flatten(P::AbstractChain)::PolygonalChain
    lengths, ang_vals, dihedrals = internalCoordinates(P)
    newDiheds = [pi for i in dihedrals]
    return PolygonalChain(lengths,ang_vals,newDiheds)
end

function dihedralRotate(P::PolygonalChain,i::Integer,phi::RealOrDual)::PolygonalChain
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
function dihedralRotateFast(P::PolygonalChain,i::Integer,phi::RealOrDual)::PolygonalChain
    n = length(P)
    if 1 <= i <= n-2
        rotMat = zrotation(phi)
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
function dihedralRotate!(P::PolygonalChain,i::Integer,phi::RealOrDual)
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

function dihedralRotateFast!(P::PolygonalChain,i::Integer,phi::RealOrDual)
    n = length(P)
    if 1 <= i <= n-2
        rotMat = zrotation(phi)
        for j in i+3:n+1
            P[j] = rotMat*(P[j]-P[i+1]) + P[i+1]
        end
    else
        error("$i must be lower than $n")
    end
end

function internalRotate(P::PolygonalChain,i::Integer,theta::RealOrDual)::PolygonalChain
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
function internalRotateFast(P::PolygonalChain,i::Integer,theta::RealOrDual)::PolygonalChain
    n = length(P)
    if 1 <= i <= n-1
        rotMat = zrotation(theta)
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
function internalRotate!(P::PolygonalChain,i::Integer,theta::RealOrDual)
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

function internalRotateFast!(P::PolygonalChain,i::Integer,theta::RealOrDual)
    n = length(P)
    if 1 <= i <= n-1
        rotMat = zrotation(theta)
        for j in i+2:n+1
            P[j] = rotMat*(P[j]-P[i+1]) + P[i+1]
        end
    else
        error("$i must be lower than $n")
    end
end

########################################################

# Important polygonal chains

########################################################

function circle(t::Real)
    return (cos(t),sin(t),0)
end

function treefoil(t::Real)
    return (cos(t)+2*cos(2*t),sin(t)-2*sin(2*t),2*sin(3t))
end

function eightKnot(t::Real)
    return (3*cos(t)+5*cos(3*t),3*sin(t)+5*sin(3*t),sin(5*t/2)*sin(3*t) + sin(4*t) - sin(6*t))
end

function eightKnotTight(t::Real)
    return (sin(t) + t/10, sin(t)*cos(t/2), sin(2*t)*sin(t/2)/4)
end

function polynomialTrefoil(t::Real)
    return (t^3 - 3*t, t^4 - 4*t^2, t^5/5 - 2*t)
end

function polynomialEight(t::Real)
    return (2/5*t*(t^2 - 7)*(t^2 - 10), t^4 - 13*t^2, 1/10*(t^7 - 31*t^5 + 164*t^3 + 560*t))
end

function parametricCurveChain(curv::Function,n::Integer,a::Real=0,b::Real=pi)
    angs = LinRange(a,b,n)
    points = [Point{Float64}(curv(ang)...) for ang in angs]
    return PolygonalChain(points)
end

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


function GenericKnittingNeedle(l::Real=2.0;ep::Real=1/6)::PolygonalChain
    p0 = ex - ep*ey
    p1 = rotate(ez,-pi/6,ex)
    p2 = e0
    p3 = ex
    p4 = ex + ez
    p5 = ep*ey
    p0 = p1 + l*unitVector(p0-p1)
    p5 = p4 + l*unitVector(p5-p4)
    points = [p0,p1,p2,p3,p4,p5]
    randP = Point()
    points = [p + 0.1*(2*randP-Point(1.0,1.0,1.0)) for p in points]
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

function makeGeneric(P::PolygonalChain)
    tras = 2*Point() - (ex_l + ey_l + ez_l)
    rot = 2*Point() - (ex_l + ey_l + ez_l)
    thetarand = 2*pi*rand()
    rotmat = rotation(thetarand,unitVector(rot))
    return PolygonalChain([rotmat*(p + tras) for p in P.vertices])
end

function makeGeneric!(P::PolygonalChain)
    tras = 2*Point() - (ex_l + ey_l + ez_l)
    rot = 2*Point() - (ex_l + ey_l + ez_l)
    thetarand = 2*pi*rand()
    rotmat = rotation(thetarand,unitVector(rot))
    n = length(P)
    for i in 1:n+1
        P[i] = rotmat*(P[i] + tras)
    end
end

function polygonal(n)
    lengths = rand(BigFloat,n)
    angles = pi*rand(BigFloat,n-1)
    diheds = pi*([2*rand(BigFloat)-1 for i in 1:n-2])
    P = PolygonalOld(lengths,angles,diheds)
    return P
end

function polygonal2(n)
    lengths = rand(BigFloat,n)
    angles = pi*rand(BigFloat,n-1)
    diheds = pi*([2*rand(BigFloat)-1 for i in 1:n-2])
    P = PolygonalChainRosetta(lengths,angles,diheds)
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
    Q = PolygonalChain(m)    
    display(@benchmark dihedralRotate!(Q,10,pi/2))
    println()
    P = PolygonalChain(m)
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
    println(typeof(e0))
end