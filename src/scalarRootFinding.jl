
const m = 1000

function generalBracketingMethod(f,epsilon,newPointMethod;a,b,nmax=m)
    if f(a)*f(b) >= 0
        error("function does not change values in interval $a,$b")
    end   
    c = newPointMethod(f,a,b)
    n = 1
    while (abs(f(c)) > epsilon) && (n < nmax)
        if f(a)*f(c) < 0
            a,b = a,c
        else
            a,b = c,b
        end
        c = newPointMethod(f,a,b)
        n += 1
    end
    return c
end
function generalBracketingMethodChain(f,epsilon,newPointMethod;a,b,nmax=m)
    if f(a)*f(b) >= 0
        error("function does not change values in interval $a,$b")
    end   
    cs = zeros(nmax+1)
    c = newPointMethod(f,a,b)
    n = 1
    cs[1] = c
    while (abs(f(c)) > epsilon) && (n < nmax)
        if f(a)*f(c) < 0
            a,b = a,c
        else
            a,b = c,b
        end
        c = newPointMethod(f,a,b)
        n += 1
        cs[n] = c
    end
    return cs[1:(n)]
end

function middlePoint(f,a,b)
    return (a+b)/2
end
function linearInterpolantZero(f,a,b)
    if f(a)*f(b) >= 0
        error("function does not change values in interval $a,$b")
    end
    return -f(a)*((b-a)/(f(b)-f(a)))+a
end

function bisection(f,epsilon;a,b,nmax=m,kwargs...)
    return generalBracketingMethod(f,epsilon,middlePoint;a=a,b=b,nmax=nmax)
end

function bisectionChain(f,epsilon;a,b,nmax=m,kwargs...)
    return generalBracketingMethodChain(f,epsilon,middlePoint;a=a,b=b,nmax=nmax)
end

function falsePosition(f,epsilon;a,b,nmax=m,kwargs...)
    return generalBracketingMethod(f,epsilon,linearInterpolantZero;a=a,b=b,nmax=nmax)
end

function falsePositionChain(f,epsilon;a,b,nmax=m,kwargs...)
    return generalBracketingMethodChain(f,epsilon,linearInterpolantZero;a=a,b=b,nmax=nmax)
end

function generalIterativeMethod(f,epsilon,newPointMethod;x0,nmax=m,kwargs...)
    # no need for copying as everything will be scalar
    x = x0
    n = 1
    while (abs(f(x)) > epsilon) && (n < nmax)
        q = newPointMethod(f,x)
        x += -inv(q)*f(x)
        n += 1
    end
    return x
end
function generalIterativeMethodChain(f,epsilon,newPointMethod;x0,nmax=m,kwargs...)
    # no need for copying as everything will be scalar
    xs = zeros(nmax)
    x = x0
    xs[1] = x0
    n = 1
    while (abs(f(x)) > epsilon) && (n < nmax)
        q = newPointMethod(f,x)
        x += -inv(q)*f(x)
        n += 1
        xs[n]=x
    end
    return xs[1:n]
end

function chordPoint(f,x,a,b)
    return (f(b)-f(a))/(b-a)
end

function newtonPoint(f,x,df)
    return df(x)
end

function chord(f,epsilon;a,b,nmax=m,kwargs...)
    func = (f,x) -> chordPoint(f,x,a,b)
    return generalIterativeMethod(f,epsilon,func,x0=(a+b)/2,nmax=nmax)
end

function chordChain(f,epsilon;a,b,nmax=m,kwargs...)
    func = (f,x) -> chordPoint(f,x,a,b)
    return generalIterativeMethodChain(f,epsilon,func,x0=(a+b)/2,nmax=nmax)
end

function newton(f,epsilon;x0,df,nmax=m,kwargs...)
    func = (f,x) -> newtonPoint(f,x,df)
    return generalIterativeMethod(f,epsilon,func,x0=x0,nmax=nmax,kwargs...)
end

function newtonChain(f,epsilon;x0,df,nmax=m,kwargs...)
    func = (f,x) -> newtonPoint(f,x,df)
    return generalIterativeMethodChain(f,epsilon,func;x0=x0,nmax=nmax)
end

# secant method requires special treatment due to its dependence on two last values
function secant(f,epsilon;x0,x1,nmax=m,kwargs...)
    n = 1
    while (abs(f(x1)) > epsilon) && (n < nmax)
        x0,x1 = x1, x1 -(x1-x0)*inv(f(x1)-f(x0))*f(x1)
        n += 1
    end
    return x1
end
function secantChain(f,epsilon;x0,x1,nmax=m,kwargs...)
    n = 1
    xs = zeros(nmax+2)
    xs[1] = x0
    xs[2] = x1
    while (abs(f(x1)) > epsilon) && (n < nmax)
        x0,x1 = x1, x1 -(x1-x0)*inv(f(x1)-f(x0))*f(x1)
        n += 1
        xs[n+1] = x1
    end
    return xs[1:(n+1)]
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    f = x -> x^2 - x - 2
    epsilon = 1e-5
    values = ["bisection","falsePosition","chord","newton","secant"]
    kwargs = Dict(:a => 1 + 0.2*randn(), :b => 3 + 0.2*randn(),:x0 => 2*rand()+1,:x1 => 2*rand()+1,:df=> x -> 2*x-1)
    p = plot()
    for i in 1:length(values)
        #println(values[i])
        method = getfield(Main,Symbol(string(values[i],"Chain")))
        A =  method(f,epsilon;kwargs...)
        plot!([abs(a-2) for a in A],markershape = :circle,markersize=4,label=values[i],xscale=:log10,yscale=:log10)
    end
    title!("Roots methods")
    xlabel!("Step")
    ylabel!("Absolut error")
    display(p)
    f = x -> sin(x) - 0.1*x^2 + 1
    df = x -> cos(x) - 0.2*x
    epsilon=1e-8
    values = ["bisection","falsePosition","chord","newton","secant"]
    kwargs = Dict(:a => 2 + 0.2*randn(), :b => 4 + 0.2*randn(),:x0 => 2*rand()+1,:x1 => 2*rand()+1,:df=> df)
    p = plot()
    for i in 1:length(values)
        #println(values[i])
        method = getfield(Main,Symbol(string(values[i],"Chain")))
        A =  method(f,epsilon;kwargs...)
        plot!([abs(a-3.1495) for a in A],markershape = :circle,markersize=4,label=values[i],xscale=:log10,yscale=:log10)
    end
    title!("Roots methods")
    xlabel!("Step")
    ylabel!("Absolut error")
    display(p)
end



