using Distributed, SharedArrays
const cores = parse(Int16,ARGS[2])
Distributed.addprocs(cores)
@everywhere include("simulations.jl")

if !isdir(dirname(ARGS[1]))
    mkdir(dirname(ARGS[1]))
end

function lsimulationPar(ls,iter::Integer,angmax::Real=pi/20,angmin::Real=-pi/20;savename="")
    n = length(ls)
    rmsds_mean = SharedArray(zeros(n))
    rmsds_error = SharedArray(zeros(n))
    ts_mean = SharedArray(zeros(n))
    ts_error = SharedArray(zeros(n))
    max_iter = parse(Int64,ARGS[5])
    @sync @distributed for i in 1:n
        temp_rmsds = zeros(iter)
        temp_ts = zeros(iter)
        for j in 1:iter
            P = stair(5)
            Q = knittingneedle(ls[i])
            #println("creacion ok")
            lastQ, angles, diheds = randomSearch(P,Q,0.5,angmax,angmin,max_iter=max_iter)
            nwork = myid()-1
            if nwork==1
                per = round(((i-1)*iter+(j-1))*100/(n); digits= 2)
                prog = "Progress: $(per) % "
                #current="Current: B, J, C, T = $(simulParam) "
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


if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    const lvals,ts_mean,ts_error,rmsds_mean,rmsds_error = lsimulationPar(LinRange(1.4,1.41,parse(Int16,ARGS[4])),parse(Int16,ARGS[3]);savename=ARGS[1])
    open(string(ARGS[1],"results.csv"),"w+") do io
        table = hcat(lvals,ts_mean,ts_error,rmsds_mean,rmsds_error)
        DelimitedFiles.writedlm(io,table,',')
    end
    println("Progress: 100 % ")
    saveLtable(joinpath(ARGS[1]),lvals)
    scatter(lvals,rmsds_mean,yerror=rmsds_error,ylabel="RMSD sobrepuesto",xlabel="l",label=false)
    savefig(joinpath(ARGS[1],"lsmiulation.pdf"))
end
