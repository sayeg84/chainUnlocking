include("argumentParsing.jl")
include("io.jl")

function runLSimulation(SimulType::MHAlgorithm,parsed_args)
    minFunc = getfield(Main,Symbol(parsed_args["minFunc"]))
    if parsed_args["log_l"]
        exps = LinRange(log10(parsed_args["lmin"]),log10(parsed_args["lmax"]),parsed_args["lvals"])
        ls = [10^x for x in exps]
    else
        ls = LinRange(parsed_args["lmin"],parsed_args["lmax"],parsed_args["lvals"])
    end
    n = length(ls)
    minfs_mean = zeros(n)
    minfs_error = zeros(n)
    ts_mean = zeros(n)
    ts_error = zeros(n)
    for i in 1:n
        n1zeros = Int(ceil(log10(n+1)))
        n1 = lpad(i,n1zeros,'0')
        P = makeChain(parsed_args["chain"],ls[i])
        initQs,finalQs,fun_vals,_ = runSingle(P,SimulType,joinpath(parsed_args["path"],string(n1,"_",)),parsed_args)
        per = round((i-1)*100/n; digits=2)
        prog = "Progress: $(per) % "
        println()
        println(prog)
        #println(current)
        println()
        temp_minfs = [minFunc(Q) for Q in finalQs]
        temp_ts = [size(fun_vals)[1] for _ in finalQs]
        minfs_mean[i] = Statistics.mean(temp_minfs)
        minfs_error[i] = Statistics.std(temp_minfs)
        ts_mean[i] = Statistics.mean(temp_ts)
        ts_error[i] = Statistics.std(temp_ts)
    end
    return ls,ts_mean,ts_error,minfs_mean,minfs_error
end

function runLSimulation(SimulType::GDAlgorithm,parsed_args)
    if parsed_args["log_l"]
        exps = LinRange(log10(parsed_args["lmin"]),log10(parsed_args["lmax"]),parsed_args["lvals"])
        ls = [10^x for x in exps]
    else
        ls = LinRange(parsed_args["lmin"],parsed_args["lmax"],parsed_args["lvals"])
    end
    n = length(ls)
    minfs_mean = zeros(n)
    minfs_error = zeros(n)
    ts_mean = zeros(n)
    ts_error = zeros(n)
    for i in 1:n
        n1zeros = Int(ceil(log10(n+1)))
        n1 = lpad(i,n1zeros,'0')
        P = makeChain(parsed_args["chain"],ls[i])
        initQ,finalQ = runSingle(P,SimulType,joinpath(parsed_args["path"],string(n1,"_",)),parsed_args)
        per = round((i-1)*100/n; digits=2)
        prog = "Progress: $(per) % "
        println()
        println(prog)
        #println(current)
        println()
        minfs_mean[i] = tangentEnergyFrac(finalQ) 
        ts_mean[i] = parsed_args["max_iter"]
    end
    return ls,ts_mean,ts_error,minfs_mean,minfs_error
end

function main()
    parsed_args = parse_commandline()
    simulation = getfield(Main,Symbol(parsed_args["algorithm"]))()
    if !isdir(parsed_args["path"])
        mkdir(parsed_args["path"])
    end
    saveMetaParams(parsed_args["path"],simulation,parsed_args)
    ls,ts_mean,ts_error,minfs_mean,minfs_error = runLSimulation(simulation,parsed_args)
    println()
    println("Progress: 100 % ")
    #println(current)
    println()
    open(joinpath(parsed_args["path"],"results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,minfs_mean,minfs_error)
        write(io,"l,t_mean,t_std,minf_mean,minf_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    saveLtable(parsed_args["path"],ls)
    scatter(ls,minfs_mean,yerror=minfs_error./2,ylabel="minFunc",xlabel="l",label=false)
    savefig(joinpath(parsed_args["path"],"lsimulation.pdf"))
end

using Plots
println("Simple simulator")
println()
main()   
