include("simulations.jl")

using DelimitedFiles
using Statistics

const base_vars = ["algorithm","minFunc","chain","indep_simuls",
"processes","max_iter","tolerance"]

function saveMetaParams(name::AbstractString,simul::MHAlgorithm,parsed_args)
    open(joinpath(name,"metaParams.csv"),"w+") do io
        sim_vars = ["max_angle","mut_k","selection",
        "temp_init","temp_f","iter_per_temp","temp_program"]
        for var in vcat(base_vars,sim_vars)
            str = string(var,",",parsed_args[var],"\n") 
            write(io,str)
        end
    end
end

function saveMetaParams(name::AbstractString,simul::GDAlgorithm,parsed_args)
    open(joinpath(name,"metaParams.csv"),"w+") do io
        sim_vars = ["time_step","ndens"]
        for var in vcat(base_vars,sim_vars)
            str = string(var,",",parsed_args[var],"\n") 
            write(io,str)
        end
    end
end

function saveChain(name::AbstractString,P::AbstractChain)
    ns = length(P)+1
    nzeros = Int(ceil(log10(ns+1)))
    pnames = [lpad(i,nzeros,'0') for i in 1:ns]
    pnames = [string("p",i) for i in pnames]
    pnames = [string(i,c) for i in pnames for c in ("x","y","z")]
    firstrow = join(pnames,',')
    firstrow = string(firstrow,"\n")
    arr = to2DArray(P)
    open(name,"w+") do io
        write(io,firstrow)
        coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
        coordinates = join(coordinates,",")
        coordinates = string(coordinates,"\n")
        write(io,coordinates)
    end
end

function saveChains(name::AbstractString,Ps::Array{<:AbstractChain,1})
    ns = length(Ps[1])+1
    nzeros = Int(ceil(log10(ns+1)))
    pnames = [lpad(i,nzeros,'0') for i in 1:ns]
    pnames = [string("p",i) for i in pnames]
    pnames = [string(i,c) for i in pnames for c in ("x","y","z")]
    firstrow = join(pnames,',')
    firstrow = string(firstrow,"\n")
    open(name,"w+") do io
        write(io,firstrow)
        for P in Ps
            arr = to2DArray(P)
            coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
            coordinates = join(coordinates,",")
            coordinates = string(coordinates,"\n")
            write(io,coordinates)
        end
    end
end

function funcValues(Q::AbstractChain,ang_idxs::Array{<:Real,1},ang_vals::Array{<:Real,1},minFunc)
    newQ = copy(Q)
    funcvals = zeros(ftype(Q),length(ang_vals))
    funcvals[1] = minFunc(newQ)
    for i in 1:(length(ang_vals)-1)
        if ang_idxs[i]!= 0
            rotate!(newQ,ang_idxs[i],ang_vals[i])
        end
        funcvals[i+1] = minFunc(newQ)
    end
    #display(toArray(newQ))
    #println()
    #display(toArray(lastQ))
    #println()
    return funcvals
end

function saveTrajectory(name::AbstractString,Q::AbstractChain,ang_vals::Array{<:Real,1},ang_idxs::Array{<:Real,1})
    ns = length(Q)+1
    nzeros = Int(ceil(log10(ns+1)))
    # padding the number with zeros to the right in order to have no problems if sorting
    pnames = [lpad(i,nzeros,'0') for i in 1:ns]
    pnames = [string("p",i) for i in pnames]
    pnames = [string(i,c) for i in pnames for c in ("x","y","z")]
    firstrow = join(pnames,',')
    firstrow = string(firstrow,"\n")
    newQ = copy(Q)
    open(name,"w+") do io
        write(io,firstrow)
        arr = to2DArray(newQ)
        coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
        coordinates = join(coordinates,",")
        coordinates = string(coordinates,"\n")
        write(io,coordinates)
        for i in 1:length(ang_vals)
            if ang_idxs[i]!= 0 
                rotate!(newQ,ang_idxs[i],ang_vals[i])
                arr = to2DArray(newQ)
                coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
                coordinates = join(coordinates,",")
                coordinates = string(coordinates,"\n")
                write(io,coordinates)
            end
        end
    end
end

function saveTrajectory(name::AbstractString,Q::AbstractChain,res::Array{<:Real,2})
    ns = length(Q)+1
    nzeros = Int(ceil(log10(ns+1)))
    # padding the number with zeros to the right in order to have no problems if sorting
    pnames = [lpad(i,nzeros,'0') for i in 1:ns]
    pnames = [string("p",i) for i in pnames]
    pnames = [string(i,c) for i in pnames for c in ("x","y","z")]
    firstrow = join(pnames,',')
    firstrow = string(firstrow,"\n")
    newQ = copy(Q)
    open(name,"w+") do io
        write(io,firstrow)
        DelimitedFiles.writedlm(io,res,',')
    end
end

function saveSimulation(name::AbstractString,Q::AbstractChain,lastQ::AbstractChain,
    ang_vals::Array{<:Real,1},ang_idxs::Array{<:Real,1};saveTrajec=true)
    saveChain(string(name,"init_Q.csv"),Q)
    saveChain(string(name,"_lastQ.csv"),lastQ)
    if saveTrajec
        saveTrajectory(string(name,"_trajectory.csv"),Q,ang_vals,ang_idxs)
    end
    open(string(name,"_ang_idxs.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ang_idxs,',')
    end
    open(string(name,"_ang_vals.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ang_vals,',')
    end
end

function saveLtable(name::AbstractString,ls)
    open(joinpath(name,"ls.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ls,',')
    end
end



function readChain(name::AbstractString)
    chain,_ = DelimitedFiles.readdlm(name,',',BigFloat,header=true)
    return PolygonalChain(chain[1,:])
end

function readSingleSimulation(name::AbstractString)
    Q = readChain(string(name,"init_Q.csv"))
    lastQ = readChain(string(name,"final_Q.csv"))
    ang_idxs = DelimitedFiles.readdlm(string(name,"_ang_idxs.csv"),',',Int16)
    ang_idxs = reshape(ang_idxs,length(ang_idxs))
    ang_vals = DelimitedFiles.readdlm(string(name,"_ang_vals.csv"),',',BigFloat)
    ang_vals = reshape(ang_vals,length(ang_vals))
    return Q, lastQ, ang_idxs, ang_vals
end


# `name` must be name of folder
function readMetaParams(name::AbstractString)
    metaParams = DelimitedFiles.readdlm(joinpath(name,"metaParams.csv"),',')
    metaParams = Dict(zip(metaParams[:,1],metaParams[:,2]))
    return metaParams
end


function readLSimulation(name::AbstractString, burnout::Real; verbose::Bool=true)
    ls = DelimitedFiles.readdlm(joinpath(name,"ls.csv"))
    ls = reshape(ls,length(ls))
    ln = length(ls)
    metaParams = readMetaParams(name)
    if verbose
        println("Minimizing function")
        println(metaParams["minFunc"])
    end
    minFunc = getfield(Main,Symbol(metaParams["minFunc"]))
    simuls  = [file for file in readdir(name)]
    # saving `_lastQ.csv` files values
    simuls = [file for file in simuls if length(file)>7]
    simuls = [file for file in simuls if file[end-8:end]=="lastQ.csv"]
    # number of total simulations
    lt = length(simuls)
    # filter independent interations
    simuls = [split(file,"_")[1] for file in simuls]
    simuls = unique(simuls)
    simuls = sort(simuls)
    ts_mean     = zeros(ln)
    ts_error    = zeros(ln)
    minfs_mean  = zeros(ln)
    minfs_error = zeros(ln)
    tentativeIters = div(lt,ln,RoundUp)
    if verbose
        println("Different simulations")
        println(lt)
        println("L vals")
        println(ln)
        println("Assumed iterations")
        println(tentativeIters)
    end
    ts_table             = zeros(ln,tentativeIters)
    minfs_table          = zeros(ln,tentativeIters)
    accepted_moves_table = zeros(ln,tentativeIters)

    for (i,sim) in enumerate(simuls)
        
        iterations = [file for file in readdir(name) if split(file,"_")[1]==sim]
        iterations = [split(file,"_")[2] for file in iterations]
        iterations = unique(iterations)
        ts_sim     = zeros(length(iterations))
        minfs_sim  = zeros(length(iterations))

        for (j,iter) in enumerate(iterations)
            verbose && println("reading simulation $(sim),$(iter)")
            lastQ  = readChain(joinpath(name,string(sim,"_",iter,"_lastQ.csv")))
            ang_idxs = DelimitedFiles.readdlm(joinpath(name,string(sim,"_",iter,"_ang_idxs.csv")),',',Int16)
            ang_idxs = reshape(ang_idxs,length(ang_idxs))
            ts_sim[j]        = length(ang_idxs)
            ts_table[i,j]    = length(ang_idxs)
            if burnout < 1
                fun_vals = DelimitedFiles.readdlm(joinpath(name,string(sim,"_",iter,"_fun_vals")),',',BigFloat)
                ncut    = Int(ceil(burnout*length(ang_idxs)))
                minfs_sim[j]     = Statistics.mean(fun_vals[ncut:end])
                minfs_table[i,j] = Statistics.mean(fun_vals[ncut:end])
            else
                minfs_sim[j]     = minFunc(lastQ)
                minfs_table[i,j] = minFunc(lastQ)
            end
            accepted_moves_table[i,j] = length([k for k in ang_idxs if k!=0])
        end
        ts_mean[i] = Statistics.mean(ts_sim)
        ts_error[i] = Statistics.std(ts_sim)
        minfs_mean[i] = Statistics.mean(minfs_sim)
        minfs_error[i] = Statistics.std(minfs_sim)
    end
    return ls, ts_mean, ts_error, minfs_mean, minfs_error, ts_table, minfs_table,accepted_moves_table
end



function makeCoordinates(name::AbstractString)
    Q = readChain(string(name,"init_Q.csv"))
    ang_vals = DelimitedFiles.readdlm(string(name,"_ang_vals.csv"))
    ang_vals = reshape(ang_vals,length(ang_vals))
    ang_idxs = DelimitedFiles.readdlm(string(name,"_ang_idxs.csv"),',',Int16)
    ang_idxs = reshape(ang_idxs,length(ang_idxs))
    saveTrajectory(string(name,"trajectory.csv"),Q,ang_vals,ang_idxs)
 end