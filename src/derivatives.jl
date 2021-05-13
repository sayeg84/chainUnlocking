include("hyperdual.jl")

function forwardDifference(f,x,h)
    return (f(x+h)-f(x))/h
end
function backwardDifference(f,x,h)
    return (f(x)-f(x-h))/h
end
function centralDifference(f,x,h)
    return (f(x+h)-f(x-h))/(2*h)
end
function centralPartialDerivative(f,x,h,i)
    e = zeros(length(x))
    e[i] = 1.0
    return (f(x+h*e)- f(x-h*e))/(2*h)
end

function centralSecondPartialDerivative(f,x,h,i,j)
    e1 = zeros(length(x))
    e2 = copy(e1)
    e1[i] = 1.0
    e2[j] = 1.0
    return (f(x+h*(e1+e2)) + f(x-h*(e1+e2))- f(x+h*(e1-e2))-f(x+h*(-e1+e2)))/(4*h^2)
end



function complexForwardPartialDerivative(f,x,h,i)
    e = zeros(Complex,length(x))
    e[i] = im
    return imag(f(x+h*e)/(h))
end

function complexForwardSecondPartialDerivative(f,x,h,i,j)
    e1 = zeros(Complex,length(x))
    e2 = copy(e1)
    e1[i] = im
    e2[j] = im
    return (real(f(x+h*e1+h*e2)) + f(x) - real(f(x+h*e1)) - real(f(x+h*e2)))/h^2
end

function complexCentralTotalDerivative(f,x,h)
    return imag((f(x+h*im)-f(x-h*im))/(2*h))
end

function complexCentralPartialDerivative(f,x,h,i)
    e = zeros(Complex,length(x))
    e[i] = im
    return imag((f(x+h*e)-f(x-h*e))/(2*h))
end

function complexCentralSecondPartialDerivative(f,x,h,i,j)
    e1 = zeros(Complex,length(x))
    e2 = copy(e1)
    e1[i] = im
    e2[j] = im
    return real(f(x+h*(e1-e2))+ f(x+h*(-e1+e2)) -f(x+h*(e1+e2))-f(x-h*(e1+e2)))/(4*h^2)
end

function hyperdualTotalDerivative(f,x,h)
    y = Hyperdual(x,h,h,0)
    return f(y).a1/h
end

function hyperdualPartialDerivative(f,x,h,i)
    d = length(x)
    y = [Hyperdual(x[j],0,0,0) for j in 1:d]
    y[i] = Hyperdual(x[i],h,h,0)
    return f(y).a1/h
end

function hyperdualSecondPartialDerivative(f,x,h,i,j)
    d = length(x)
    y = [Hyperdual(x[j],0,0,0) for j in 1:d]
    y[i] = Hyperdual(x[i],h,h,0)
    return (f(y).a3)/(h^2)
end

function jacobian(f,x,h;df=centralPartialDerivative)
    return [df(f,x,h,i) for i in 1:length(x)]
end

function hessian(f,x,h;ddf=centralSecondPartialDerivative)
    n = length(x)
    return [ddf(f,x,h/2,i,j) for i in 1:n, j in 1:n]
end

if abspath(PROGRAM_FILE) == @__FILE__
    f(x) = cos(x[1])*exp(-x[2])
    h = 1e-4
    jacobian(f,[0.0,0.0],h)
    hessian(f,[0.0,0.0],h)
    jacobian(f,[0.0,0.0],h,df=complexCentralPartialDerivative)
    hessian(f,[0.0,0.0],h,ddf=complexCentralSecondPartialDerivative)
    hessian(f,[0.0,0.0],1e-3,ddf=complexForwardSecondPartialDerivative)
    hyperdualPartialDerivative(x->-exp(-x[1]),[1],1e-3,1)
    hyperdualSecondPartialDerivative(x->-exp(-x[1]),[1],1e-3,1,1)
    """
    Computation must be done following certain order to avoid numerical problems.

    Check Jeffrey Alan Fike PhD thesis.

    https://stacks.stanford.edu/file/druid:jw107zn5044/JeffFike_thesis_online_twosided-augmented.pdf
    """
    function prueba(x)
        t1 = exp(x[1])
        t2 = sin(x[1])
        t3 = t2^3
        t4 = cos(x[1])
        t5 = t4^3
        t6 = t3+t5
        t7 = t6^(-0.5)
        return t1*t7
    end

    
    # this way also avoids problems
    function prueba(x)
        return exp(x[1])*(sin(x[1])^3 + cos(x[1])^3)^(-0.5)
    end

    function dprueba(x)
        return exp(x[1])*(3*cos(x[1])+5*cos(3*x[1])+9*sin(x[1])+sin(3*x[1]))/(8*(cos(x[1])^3+sin(x[1])^3)^(3/2))    
    end

    function ddprueba(x)
        a = -12*cos(2*x[1]) + 30*cos(4*x[1]) + 12*cos(6*x[1])
        b = -111*sin(2*x[1]) + 48*sin(4*x[1]) + 5*sin(6*x[1])
        c = 64*(cos(x[1])^3 + sin(x[1])^3)^(5/2)
        return exp(x[1])*(130+a+b)/c
    end

    using Plots

    x0 = [1.5]
    r = dprueba(x0)
    hs = range(-18,stop=-1,length=40)
    hs = [10^(h) for h in hs]
    p = plot()
    testFuncs = ["central","complexCentral","hyperdual"]
    for f in testFuncs
        method = getfield(Main,Symbol(string(f,"PartialDerivative")))
        A =  [jacobian(prueba,x0,h,df=method)[1] for h in hs]
        plot!(hs,[abs(a-r)+1e-15 for a in A],markershape = :circle,markersize=4,label=f,xscale=:log10,yscale=:log10)
    end
    display(p)

    x0 = [1.5]
    r = ddprueba(x0)
    hs = range(-15,stop=-1,length=30)
    hs = [10^a for a in hs]
    p = plot()
    testFuncs = ["central","complexCentral","hyperdual"]
    for f in testFuncs
        method = getfield(Main,Symbol(string(f,"SecondPartialDerivative")))
        A =  [hessian(prueba,x0,h,ddf=method)[1,1] for h in hs]
        plot!(hs,[abs(a-r)+1e-15 for a in A],markershape = :circle,markersize=4,label=f,xscale=:log10,yscale=:log10)
    end
    display(p)

end
