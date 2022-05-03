
using ArgParse
using Plots

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--path"
            help = "Folder to save the simulations"
            arg_type = String
            required = true
        "--burnout"
            help = "Fraction of burnout values. Value of 1 saves only the last one"
            arg_type = Float64
            default  = 1.0
        "--processes" 
            help = "number of processs to run"
            arg_type = Int
            default = 4
        "--minFunc"
            help = "minFunc to make analysis of"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

const parsed_args = parse_commandline()

using Distributed
using SharedArrays
Distributed.addprocs(parsed_args["processes"])

@everywhere include("io.jl")

function readLSimulationPar(name::AbstractString, burnout::Real,minFuncStr::AbstractString=""; verbose::Bool=true)
    ls = DelimitedFiles.readdlm(joinpath(name,"ls.csv"))
    ls = reshape(ls,length(ls))
    ln = length(ls)
    metaParams = readMetaParams(name)
    simul   = getfield(Main,Symbol(metaParams["algorithm"]))()
    if isempty(minFuncStr)
        minFunc = getfield(Main,Symbol(metaParams["minFunc"]))
    else
        minFunc = getfield(Main,Symbol(minFuncStr))
    end
    if verbose
        println("Minimizing function")
        println(string(minFunc))
    end
    indep_simuls = typeof(simul) <: MHAlgorithm ? metaParams["indep_simuls"] : 1
    # saving `_lastQ.csv` files values
    ts_table             = SharedArray(zeros(ln,indep_simuls))
    minfs_table          = SharedArray(zeros(ln,indep_simuls))
    accepted_moves_table = SharedArray(zeros(ln,indep_simuls))
    @sync @distributed for i in 1:ln
        println("Reading lval = $i")
        n1zeros = Int(ceil(log10(ln+1)))
        n1 = lpad(i,n1zeros,'0')
        Qs,funvals,accepted,ts = readSingleSimulation(joinpath(name,n1),simul,minFunc,burnout)
        minfs_table[i,:] = Statistics.mean(funvals,dims=1)
        accepted_moves_table[i,:] = accepted
        ts_table[i,:] = ts
    end
    return ls, ts_table, minfs_table, accepted_moves_table
end

function main()    
    ls,ts_table,minfs_table,accepted_moves_table = readLSimulationPar(parsed_args["path"],parsed_args["burnout"],parsed_args["minFunc"])
    println("saving results")
    ts_mean = Statistics.mean(ts_table,dims=2)
    ts_error = Statistics.std(ts_table,dims=2)
    minfs_mean = Statistics.mean(minfs_table,dims=2)
    minfs_error = Statistics.std(minfs_table,dims=2)
    open(joinpath(parsed_args["path"],"new_results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,minfs_mean,minfs_error)
        write(io,"l,t_mean,t_std,minf_mean,minf_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    open(joinpath(parsed_args["path"],"ts_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ts_table,',')
    end
    open(joinpath(parsed_args["path"],"minfs_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,minfs_table,',')
    end
    open(joinpath(parsed_args["path"],"accepted_moves_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,accepted_moves_table,',')
    end
    
    println("Making Plots")
    scatter(ls,ts_mean,yerror=ts_error./2,label=false,title="L analysis",ylabel="steps",xlabel="L")
    savefig(joinpath(parsed_args["path"],"l-t.pdf"))
    scatter(ls,minfs_mean,yerror=minfs_error./2,label=false,title="L analysis",ylabel="minf",xlabel="L")
    savefig(joinpath(parsed_args["path"],"l-minf.pdf"))
    println("Done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end