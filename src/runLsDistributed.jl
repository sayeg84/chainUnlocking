include("argumentParsing.jl")
# importing `io` and other settings have to be done after getting the 
# args due to variable number of processes.
const parsed_args = parse_commandline()
using Distributed, SharedArrays
Distributed.addprocs(parsed_args["processes"])

@everywhere include("io.jl")

let z = getfield(Main,Symbol(parsed_args["algorithm"]))
    @sync @everywhere const simulation = $z()
end

let w = getfield(Main,Symbol(parsed_args["minFunc"]))
    @sync @everywhere const minFunc = $w
end

let u = getfield(Main,Symbol(parsed_args["chain"]))
    @sync @everywhere const chainFunc = $u
end

function runLSimulationDistributed(SimulType::MHAlgorithm,parsed_args)
    if parsed_args["log_l"]
        exps = LinRange(log10(parsed_args["lmin"]),log10(parsed_args["lmax"]),parsed_args["lvals"])
        ls = [10^x for x in exps]
    else
        ls = LinRange(parsed_args["lmin"],parsed_args["lmax"],parsed_args["lvals"])
    end
    n = length(ls)
    minfs_mean  = SharedArray(zeros(n))
    minfs_error = SharedArray(zeros(n))
    ts_mean  = SharedArray(zeros(n))
    ts_error = SharedArray(zeros(n))
    @sync @distributed for i in 1:n
        n1zeros = Int(ceil(log10(n+1)))
        n1 = lpad(i,n1zeros,'0')
        P = chainFunc(ls[i])
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

function runLSimulationDistributed(SimulType::GDAlgorithm,parsed_args)
    if parsed_args["log_l"]
        exps = LinRange(log10(parsed_args["lmin"]),log10(parsed_args["lmax"]),parsed_args["lvals"])
        ls = [10^x for x in exps]
    else
        ls = LinRange(parsed_args["lmin"],parsed_args["lmax"],parsed_args["lvals"])
    end
    n = length(ls)
    minfs_mean  = SharedArray(zeros(n))
    minfs_error = SharedArray(zeros(n))
    ts_mean  = SharedArray(zeros(n))
    ts_error = SharedArray(zeros(n))
    @sync @distributed for i in 1:n
        n1zeros = Int(ceil(log10(n+1)))
        n1 = lpad(i,n1zeros,'0')
        P = chainFunc(ls[i])
        initQ,finalQ = runSingle(P,SimulType,joinpath(parsed_args["path"],string(n1,"_",)),parsed_args)
        per = round((i-1)*100/n; digits=2)
        prog = "Progress: $(per) % "
        println()
        println(prog)
        #println(current)
        println()
        minfs_mean[i] = minFunc(finalQ) 
        ts_mean[i] = parsed_args["max_iter"]
    end
    return ls,ts_mean,ts_error,minfs_mean,minfs_error
end

function main()
    if !isdir(parsed_args["path"])
        mkdir(parsed_args["path"])
    end
    ls,ts_mean,ts_error,minfs_mean,minfs_error = runLSimulationDistributed(simulation,parsed_args)
    saveMetaParams(parsed_args["path"],simulation,parsed_args)
    open(joinpath(parsed_args["path"],"results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,minfs_mean,minfs_error)
        write(io,"l,t_mean,t_std,minf_mean,minf_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    println("Progress: 100 % ")
    # saving info
    saveLtable(parsed_args["path"],ls)
    # making plots
    scatter(ls,minfs_mean,yerror=minfs_error./2,ylabel="minFunc",xlabel="l",label=false)
    savefig(joinpath(parsed_args["path"],"lsimulation.pdf"))
end

using Plots
println("Distributed simulator")
println()
main()    
