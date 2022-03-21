include("algorithms.jl")

using DelimitedFiles, Statistics

function rotate!(P::AbstractChain,ang_idx::Integer,alpha::Real)
    n = length(P)
    if 1 <= ang_idx <= n-2
        dihedralRotate!(P,ang_idx,alpha)
    elseif n-1 <= ang_idx <= 2*n-3
        ang_idx = ang_idx - n + 2 # adjusting to be smaller 
        internalRotate!(P,ang_idx,alpha)
    else
        error("index $(ang_idx) is not supported")
    end
end

function saveChain(name::AbstractString,P::AbstractChain)
    arr = toArray(P)
    open(name,"w+") do io
        DelimitedFiles.writedlm(io,arr,',')
    end
end

function funcValues(Q,ang_idxs,ang_vals,minFunc,lastQ)
    newQ = copy(Q)
    funcvals = zeros(T,length(ang_vals))
    funcvals[1] = minFunc(newQ)
    for i in 1:(length(ang_vals)-1)
        if ang_idxs[i]!= 0
            newQ = moveBeforeDihedral(newQ,ang_idxs[i])
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


function saveTrajectory(name,Q::AbstractChain,ang_vals,ang_idxs)
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
        arr = toArray(newQ)
        coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
        coordinates = join(coordinates,",")
        coordinates = string(coordinates,"\n")
        write(io,coordinates)
        for i in 1:length(ang_vals)
            if ang_idxs[i]!= 0 
                rotate!(newQ,ang_idxs[i],ang_vals[i])
                arr = toArray(newQ)
                coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
                coordinates = join(coordinates,",")
                coordinates = string(coordinates,"\n")
                write(io,coordinates)
            end
        end
    end
end

function saveSimulation(name,Q,lastQ,ang_vals,ang_idxs;saveTrajec=true)
    saveChain(string(name,"_Q.csv"),Q)
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

function saveLtable(name::AbstractString,ls::AbstractArray)
    open(joinpath(name,"ls.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ls,',')
    end
end

function saveMetaParams(name::AbstractString,parsed_args)
    open(joinpath(name,"metaParams.csv"),"w+") do io
        variables = ["chain","algorithm","minFunc","tolerance",
                     "temp_init","max_iter","max_angle"]
        for var in variables
            str = string(var,",",parsed_args[var],"\n") 
            write(io,str)
        end
    end
end

function readChain(name::AbstractString)
    chain = DelimitedFiles.readdlm(name,',',T)
    return PolygonalChain(chain)
end

function readSingleSimulation(name::AbstractString)
    Q = readChain(string(name,"_Q.csv"))
    lastQ = readChain(string(name,"_lastQ.csv"))
    ang_idxs = DelimitedFiles.readdlm(string(name,"_ang_idxs.csv"),',',Int16)
    ang_idxs = reshape(ang_idxs,length(ang_idxs))
    ang_vals = DelimitedFiles.readdlm(string(name,"_ang_vals.csv"),',',T)
    ang_vals = reshape(ang_vals,length(ang_vals))
    return Q, lastQ, ang_idxs, ang_vals
end


# `name` must be name of folder
function readMetaParams(name::AbstractString)
    metaParams = DelimitedFiles.readdlm(joinpath(name,"metaParams.csv"),',')
    metaParams = Dict(zip(metaParams[:,1],metaParams[:,2]))
    return metaParams
end


function readLSimulation(name::AbstractString,burnout::Real; verbose::Bool=true)
    ls = DelimitedFiles.readdlm(joinpath(name,"ls.csv"))
    ls = reshape(ls,length(ls))
    ln = length(ls)
    metaParams = readMetaParams(name)
    verbose && println("minimizing function")
    verbose && println(metaParams["minFunc"])
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
    verbose && println("Diferent simulations")
    verbose && println(lt)
    verbose && println("L vals")
    verbose && println(ln)
    verbose && println("Assumed iterations")
    verbose && println(tentativeIters)

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
                fun_vals = DelimitedFiles.readdlm(joinpath(name,string(sim,"_",iter,"_fun_vals")),',',T)
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



function makeCoordinates(name)
    Q = readChain(string(name,"_Q.csv"))
    ang_vals = DelimitedFiles.readdlm(string(name,"_ang_vals.csv"))
    ang_vals = reshape(ang_vals,length(ang_vals))
    ang_idxs = DelimitedFiles.readdlm(string(name,"_ang_idxs.csv"),',',Int16)
    ang_idxs = reshape(ang_idxs,length(ang_idxs))
    saveTrajectory(string(name,"_trajectory.csv"),Q,ang_vals,ang_idxs)
 end