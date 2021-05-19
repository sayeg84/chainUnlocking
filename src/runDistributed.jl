using  ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--path"
            help = "Folder to save the simulations"
            arg_type = String
            required = true
        "--indep_simuls"
            help = "Number of simulations per l value"
            arg_type = Int
            default = 4
        "--processes"
            help = "Parallel processes to repart the simulations"
            arg_type = Int
            default = 4
        "--lmin"
            help = "minimum value for l"
            arg_type = Float64
            default = 1.4
        "--lmax"
            help = "maximum value for l"
            arg_type = Float64
            default = 1.42
        "--l_vals"
            help = "Values for the L interval"
            arg_type = Int
            default = 12
        "--log_l"
            help = "flag to tell if logarithmic space must be filled for the l Values"
            action = :store_true
        "--max_iter"
            help = "maximum number of iterations"
            arg_type = Int
            default = 10_000
        
        
    end
    return parse_args(s)
end

const parsed_args = parse_commandline()
using Distributed, SharedArrays
Distributed.addprocs(parsed_args["processes"])
@everywhere include("algorithms.jl")

function lsimulationPar(ls,iter::Integer,angmax::Real=pi/20,angmin::Real=-pi/20;savename="")
    n = length(ls)
    rmsds_mean = SharedArray(zeros(n))
    rmsds_error = SharedArray(zeros(n))
    ts_mean = SharedArray(zeros(n))
    ts_error = SharedArray(zeros(n))
    @sync @distributed for i in 1:n
        temp_rmsds = zeros(iter)
        temp_ts = zeros(iter)
        for j in 1:iter
            P = stair(5)
            Q = knittingneedle(ls[i])
            #println("creacion ok")
            lastQ, angles, diheds = randomSearch(P,Q,0.5,angmax,angmin,max_iter=parsed_args["max_iter"])
            nwork = myid()-1
            if nwork==1
                per = round(((i-1)*iter+(j-1))*parsed_args["processes"]*100/(parsed_args["indep_simuls"]*n); digits= 2)
                prog = "Progress: $(per) % "
                println()
                println(prog)
                #println(current)
                println()
            end
            if !isempty(savename)
                n1zeros = Int(ceil(log10(n+1)))
                n1 = lpad(i,n1zeros,'0')
                n2zeros = Int(ceil(log10(iter+1)))
                n2 = lpad(j,n2zeros,"0")
                saveSimulation(string(savename,n1,"_",n2,),P,Q,lastQ,angles,diheds,saveTrajec=false)
            end
            temp_ts[j] = length(diheds)
	    #println(temp_ts[j])
            temp_rmsds[j] = overlapedRmsd(P,lastQ)
	    #println(temp_rmsds[j])
        end
	#println(temp_rmsds)
	#println(temp_ts)
        rmsds_mean[i] = Statistics.mean(temp_rmsds)
        rmsds_error[i] = Statistics.std(temp_rmsds)
        ts_mean[i] = Statistics.mean(temp_ts)
        ts_error[i] = Statistics.std(temp_ts)
    end
    return ls,ts_mean,ts_error,rmsds_mean,rmsds_error
end

function main()
    if !isdir(dirname(parsed_args["path"]))
        mkdir(dirname(parsed_args["path"]))
    end
    if parsed_args["log_l"]
        exps = LinRange(log10(parsed_args["lmin"]),log10(parsed_args["lmax"]),parsed_args["l_vals"])
        lvals = [10^x for x in exps]
    else
        lvals = LinRange(parsed_args["lmin"],parsed_args["lmax"],parsed_args["l_vals"])
    end
    lvals,ts_mean,ts_error,rmsds_mean,rmsds_error = lsimulationPar(lvals,parsed_args["indep_simuls"];savename=parsed_args["path"])
    open(string(parsed_args["path"],"results.csv"),"w+") do io
        table = hcat(lvals,ts_mean,ts_error,rmsds_mean,rmsds_error)
        write(io,"l,t_mean,t_std,rmsd_mean,rmsd_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    println("Progress: 100 % ")
    saveLtable(parsed_args["path"],lvals)
    scatter(lvals,rmsds_mean,yerror=rmsds_error,ylabel="RMSD sobrepuesto",xlabel="l",label=false)
    savefig(joinpath(parsed_args["path"],"lsimulation.pdf"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    main()    
end
