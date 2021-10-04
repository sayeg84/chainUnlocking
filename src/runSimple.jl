include("argumentParsing.jl")
include("io.jl")

function lsimulationPar(ls,iter::Integer,angmax::Real=pi/20,angmin::Real=-pi/20;parsed_args,savename="")
    n = length(ls)
    algoFunc = getfield(Main,Symbol(parsed_args["algorithm"]))
    minFunc = getfield(Main,Symbol(parsed_args["minFunc"]))
    chainFunc = getfield(Main,Symbol(parsed_args["chain"]))
    minfs_mean = zeros(n)
    minfs_error = zeros(n)
    ts_mean = zeros(n)
    ts_error = zeros(n)
    for i in 1:n
        temp_minfs = zeros(iter)
        temp_ts = zeros(iter)
        for j in 1:iter
            Q = chainFunc(ls[i])
            #println("creacion ok")
            lastQ, angles, diheds, minvals = algoFunc(Q,minFunc,parsed_args["tolerance"],angmax,angmin,temp_init=parsed_args["temp_init"],max_iter=parsed_args["max_iter"])
            per = round(((i-1)*iter+(j-1))*100/(parsed_args["indep_simuls"]*n); digits= 2)
            prog = "Progress: $(per) % "
            println()
            println(prog)
            #println(current)
            println()
            if !isempty(savename)
                n1zeros = Int(ceil(log10(n+1)))
                n1 = lpad(i,n1zeros,'0')
                n2zeros = Int(ceil(log10(iter+1)))
                n2 = lpad(j,n2zeros,"0")
                saveSimulation(joinpath(savename,string(n1,"_",n2,)),Q,lastQ,angles,diheds,saveTrajec=false)
            end
            temp_ts[j] = length(diheds)
	    #println(temp_ts[j])
            temp_minfs[j] = minFunc(lastQ)
	    #println(temp_minfs[j])
        end
	#println(temp_minfs)
	#println(temp_ts)
        minfs_mean[i] = Statistics.mean(temp_minfs)
        minfs_error[i] = Statistics.std(temp_minfs)
        ts_mean[i] = Statistics.mean(temp_ts)
        ts_error[i] = Statistics.std(temp_ts)
    end
    return ls,ts_mean,ts_error,minfs_mean,minfs_error
end

function main()
    parsed_args = parse_commandline()

    if !isdir(parsed_args["path"])
        mkdir(parsed_args["path"])
    end
    if parsed_args["log_l"]
        exps = LinRange(log10(parsed_args["lmin"]),log10(parsed_args["lmax"]),parsed_args["lvals"])
        ls = [10^x for x in exps]
    else
        ls = LinRange(parsed_args["lmin"],parsed_args["lmax"],parsed_args["lvals"])
    end
    ls,ts_mean,ts_error,minfs_mean,minfs_error = lsimulationPar(ls,parsed_args["indep_simuls"];parsed_args=parsed_args,savename=parsed_args["path"])
    open(joinpath(parsed_args["path"],"results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,minfs_mean,minfs_error)
        write(io,"l,t_mean,t_std,minf_mean,minf_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    println("Progress: 100 % ")
    # saving info
    saveMetaParams(parsed_args["path"],parsed_args)
    saveLtable(parsed_args["path"],ls)
    scatter(ls,minfs_mean,yerror=minfs_error./2,ylabel="minFunc",xlabel="l",label=false)
    savefig(joinpath(parsed_args["path"],"lsimulation.pdf"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    main()    
end
