
struct Hyperdual{T<:AbstractFloat}
    a0::T
    a1::T
    a2::T
    a3::T
    Hyperdual{T}() where {T<:AbstractFloat} = new()
    Hyperdual{T}(a0,a1,a2,a3) where {T<:AbstractFloat} = new(a0,a1,a2,a3)
    Hyperdual(a0::AbstractIrrational,a1,a2,a3) = Hyperdual{Float64}(a0,a1,a2,a3)
    Hyperdual(a0::AbstractFloat,a1,a2,a3) = Hyperdual{typeof(a0)}(a0,a1,a2,a3)
    Hyperdual(a0::AbstractFloat,h) = Hyperdual{typeof(a0)}(a0,h,h,0)
    Hyperdual(a0,a1,a2,a3) = Hyperdual{Float64}(a0,a1,a2,a3)
    Hyperdual(a0) = Hyperdual{Float64}(a0,0.0,0.0,0.0)
end

const FloatOrDual = Union{AbstractFloat,Hyperdual}
const RealOrDual  = Union{Real,Hyperdual}

function Base.abs(x::Hyperdual)
    return abs(x.a0)
end

function Base.:+(x::Hyperdual,y::Hyperdual)
    return Hyperdual(x.a0+y.a0,x.a1+y.a1,x.a2+y.a2,x.a3+y.a3)
end

function Base.:+(x::Real,y::Hyperdual)
    return Hyperdual(x+y.a0,y.a1,y.a2,y.a3)
end

function Base.:+(x::Hyperdual,y::Real)
    return y+x
end

function Base.:-(x::Hyperdual,y::Hyperdual)
    return Hyperdual(x.a0-y.a0,x.a1-y.a1,x.a2-y.a2,x.a3-y.a3)
end

function Base.:-(x::Real,y::Hyperdual)
    return Hyperdual(x-y.a0,-y.a1,-y.a2,-y.a3)
end

function Base.:-(y::Hyperdual)
    return 0-y
end

function Base.:-(x::Hyperdual,y::Real)
    return -y+x
end

function Base.:*(x::Hyperdual,y::Hyperdual)
    return Hyperdual(x.a0*y.a0, x.a1*y.a0 + x.a0*y.a1, x.a2*y.a0 + x.a0*y.a2, x.a3*y.a0 + x.a2*y.a1 + x.a1*y.a2 + x.a0*y.a3)
end

function Base.:*(x::Real,y::Hyperdual)
    return Hyperdual(x*y.a0, x*y.a1,  x*y.a2, x*y.a3)
end

function Base.:*(x::Hyperdual,y::Real)
    return y*x
end

function Base.inv(x::Hyperdual)
    if x.a0 == 0
        error("real part must be non zero")
    end
    return Hyperdual(1/x.a0,-x.a1/(x.a0^2), -x.a2/(x.a0^2), -(2*x.a1*x.a2)/(x.a0^3) + x.a3/(x.a0^2))
end

function Base.:/(x::Hyperdual,y::Hyperdual)
    return x*inv(y)
end

function Base.:/(x::Real,y::Hyperdual)
    return x*inv(y)
end

function Base.:/(x::Hyperdual,y::Real)
    return x*(1/y)
end

function Base.:^(x::Hyperdual,n::Integer)
    if n==0
        return Hyperdual(1,0,0,0)
    elseif n>0
        res = Hyperdual(1,0,0,0)
        for i in 1:n
            res *= x
        end
        return res
    else
        return 1/(x^(-n))
    end
end

function Base.sqrt(x::Hyperdual)
    return Hyperdual(sqrt(x.a0),x.a1/(2*sqrt(x.a0)), x.a2/(2*sqrt(x.a0)) ,x.a3/(2*sqrt(x.a0))-x.a1*x.a2/(4*sqrt(x.a0)^3))
end

function Base.exp(x::Hyperdual)
    return Hyperdual(exp(x.a0),x.a1*exp(x.a0), x.a2*exp(x.a0) ,x.a3*exp(x.a0)+x.a1*x.a2*exp(x.a0))
end

function Base.log(x::Hyperdual)
    return Hyperdual(log(x.a0),x.a1/x.a0,x.a2/x.a0,x.a3/x.a0-x.a1*x.a2/(x.a0^2))
end

function Base.:^(x::Hyperdual,y::Real)
    return exp(y*log(x))
end

function Base.:^(x::Real,y::Hyperdual)
    return exp(y*log(x))
end

function generalFunctionEval(x::Hyperdual,f,df,ddf)
    b0 = f(x.a0)
    b1 = x.a1*df(x.a0)
    b2 = x.a2*df(x.a0)
    b3 = x.a3*df(x.a0) + x.a1*x.a2*ddf(x.a0)
    return Hyperdual(b0,b1,b2,b3)
end

# for computing trigonometric functions, using only sine and cosine is generally faster
function Base.sin(x::Hyperdual)
    return generalFunctionEval(x,sin,cos,x->-sin(x))
end

function Base.cos(x::Hyperdual)
    return generalFunctionEval(x,cos,x->-sin(x),x->-cos(x))
end

function Base.tan(x::Hyperdual)
    df(x) = 1/(cos(x))^2
    ddf(x) = 2*sin(x)/(cos(x)^3)
    return generalFunctionEval(x,x->sin(x)/cos(x),df,ddf)
end

function Base.cot(x::Hyperdual)
    df(x) = -1/(sin(x))^2
    ddf(x) = 2*cos(x)/(sin(x)^3)
    return generalFunctionEval(x,x->cos(x)/sin(x),df,ddf)
end

function Base.sec(x::Hyperdual)
    df(x) = sin(x)/(cos(x))^2 
    ddf(x) = (1+sin(x)^2)/(cos(x)^3)
    return generalFunctionEval(x,x->1/cos(x),df,ddf)
end

function Base.csc(x::Hyperdual)
    df(x) = -cos(x)/(sin(x))^2
    ddf(x) = (1+cos(x)^2)/(sin(x)^3)
    return generalFunctionEval(x,x->1/sin(x),df,ddf)
end

#inverse trigonometric functions
function Base.asin(x::Hyperdual)
    df(x) = 1/sqrt(1-x^2)
    ddf(x) = x/(1-x^2)^(3/2)
    return generalFunctionEval(x,asin,df,ddf)
end

function Base.acos(x::Hyperdual)
    df(x)= -1/sqrt(1-x^2)
    ddf(x) = -x/(1-x^2)^(3/2)
    return generalFunctionEval(x,acos,df,ddf)
end

function Base.atan(x::Hyperdual)
    df(x) = 1/(1+x^2)
    ddf(x) = -2*x/(1+x^2)^2
    return generalFunctionEval(x,atan,df,ddf)
end

function Base.acot(x::Hyperdual)
    df(x) = -1/(1+x^2)
    ddf(x) = 2*x/(1+x^2)^2
    return generalFunctionEval(x,acot,df,ddf)
end

function Base.asec(x::Hyperdual)
    df(x) = 1/(sqrt(1-1/x^2)*x^2)
    ddf(x) = -1/((1-1/x^2)^(3/2)*x^5) - 2/((1-1/x^2)*x^3)
    return generalFunctionEval(x,asec,df,ddf)
end

function Base.acsc(x::Hyperdual)
    df(x) = -1/(sqrt(1-1/x^2)*x^2)
    ddf(x) = 1/((1-1/x^2)^(3/2)*x^5) + 2/((1-1/x^2)*x^3)
    return generalFunctionEval(x,acsc,df,ddf)
end

# hyperbolic trigonometric functions

function Base.sinh(x::Hyperdual)
    f(x) = (exp(x)-exp(-x))/2
    df(x) = (exp(x)+exp(-x))/2
    ddf(x) =(exp(x)-exp(-x))/2
    return generalFunctionEval(x,f,df,ddf)
end

function Base.cosh(x::Hyperdual)
    f(x) = (exp(x)+exp(-x))/2
    df(x) = (exp(x)-exp(-x))/2
    ddf(x) =(exp(x)+exp(-x))/2
    return generalFunctionEval(x,f,df,ddf)
end

function Base.tanh(x::Hyperdual)
    f(x) = (exp(x)-exp(-x))/((exp(x)+exp(-x)))
    df(x) = 4/(exp(x)+exp(-x))^2
    ddf(x) =-8*(exp(x)-exp(-x))/(exp(x)+exp(-x))^3
    return generalFunctionEval(x,f,df,ddf)
end

function Base.coth(x::Hyperdual)
    f(x) = (exp(x)+exp(-x))/((exp(x)-exp(-x)))
    df(x) = -4/(exp(x)-exp(-x))^2
    ddf(x) = 8*(exp(x)+exp(-x))/(exp(x)-exp(-x))^3
    return generalFunctionEval(x,f,df,ddf)
end

function Base.sech(x::Hyperdual)
    f(x) = 2/(exp(x)+exp(-x))
    df(x) = -2*(exp(x)-exp(-x))/((exp(x)+exp(-x))^2)
    ddf(x) = (-8+2*(exp(x)-exp(-x))^2/4)/((exp(x)+exp(-x))^3)
    return generalFunctionEval(x,f,df,ddf)
end

function Base.csch(x::Hyperdual)
    f(x) = 2/(exp(x)-exp(-x))
    df(x) = -2*(exp(x)+exp(-x))/((exp(x)-exp(-x))^2)
    ddf(x) = (8+2*(exp(x)+exp(-x))^2/4)/((exp(x)-exp(-x))^3)
    return generalFunctionEval(x,f,df,ddf)
end

#inverse hyperbolic trigonometric functions

function Base.asinh(x::Hyperdual)
    df(x) = 1/sqrt(1+x^2)
    ddf(x) = -x/(1+x^2)^(3/2)
    return generalFunctionEval(x,asinh,df,ddf)
end

function Base.acosh(x::Hyperdual)
    df(x)= 1/(sqrt(-1+x)*sqrt(1+x))
    ddf(x) = -1/(2*sqrt(-1+x)*(1+x)^(3/2)) -1/(2*sqrt(1+x)*(-1+x)^(3/2))
    return generalFunctionEval(x,acosh,df,ddf)
end

function Base.atanh(x::Hyperdual)
    df(x) = 1/(1-x^2)
    ddf(x) = 2*x/(1-x^2)^2
    return generalFunctionEval(x,atanh,df,ddf)
end

function Base.acoth(x::Hyperdual)
    df(x) = -1/(1+x^2)
    ddf(x) = 2*x/(1+x^2)^2
    return generalFunctionEval(x,acot,df,ddf)
end

function Base.asech(x::Hyperdual)
    df(x) = sqrt((1-x)/(1+x))/((x-1)*x)
    ddf(x) = (-1+2*x^2)/((x-1)*x^2*sqrt((1-x)/(1+x))*(1+x)^2)
    return generalFunctionEval(x,asech,df,ddf)
end

function Base.acsch(x::Hyperdual)
    df(x) = -1/(sqrt(1+1/x^2)*x^2)
    ddf(x) = -1/((1+1/x^2)^(3/2)*x^5) + 2/((1+1/x^2)*x^3)
    return generalFunctionEval(x,acsch,df,ddf)
end


# logical comparisons:
function Base.:>(x::Hyperdual,y::Hyperdual)
    return x.a0 > y.a0
end

function Base.:>=(x::Hyperdual,y::Hyperdual)
    return x.a0 >= y.a0
end

function Base.:<(x::Hyperdual,y::Hyperdual)
    return x.a0 < y.a0
end

function Base.:<=(x::Hyperdual,y::Hyperdual)
    return x.a0 <= y.a0
end

#=
function Base.:!=(x::Hyperdual,y::Hyperdual)
    return x.a0 != y.a0
end
=#

import Base.==
function ==(x::Hyperdual,y::Hyperdual)
    return x.a0 == y.a0
end



if abspath(PROGRAM_FILE) == @__FILE__
    x = Hyperdual(1,2,3,4)
    y = Hyperdual(5,6,7,8)
    z = Hyperdual(2,0,0,0)
    w = Hyperdual(pi,2*pi,0,0)
    w = Hyperdual(3.0,2*pi,0,0)
    w = Hyperdual(3.0,BigFloat(2*pi),0,0)
    u = Hyperdual(0.5,0.5,0.5,0.5)
    println(x+y)
    println(x-y)
    println(x*y)
    println(x/y)
    println(exp(x))
    println(sin(w))
    println(cos(w))
    println(tan(w))
    println(cot(w))
    println(sec(w))
    println(csc(w))
    println(asin(u))
    println(acos(u))
    println(atan(u))
    println(acot(u))
    println(asec(w))
    println(acsc(w))
    println(sinh(w))
    println(cosh(w))
    println(tanh(w))
    println(coth(w))
    println(sech(w))
    println(csch(w))
    println(asinh(u))
    println(acosh(w))
    println(atanh(u))
    println(acoth(u))
    println(asech(u))
    println(acsch(u))
end
