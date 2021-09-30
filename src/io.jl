include("algorithms.jl")

using DelimitedFiles, Statistics

function saveChain(name::AbstractString,P::AbstractChain)
    arr = toArray(P)
    open(name,"w+") do io
        DelimitedFiles.writedlm(io,arr,',')
    end
end

function saveTrajectory(name,Q::AbstractChain,angles,diheds)
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
        for i in 1:length(angles)
            if diheds[i]!= 0 
                dihedralRotate!(newQ,diheds[i],angles[i])
                arr = toArray(newQ)
                coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
                coordinates = join(coordinates,",")
                coordinates = string(coordinates,"\n")
                write(io,coordinates)
            end
        end
    end
end

function saveSimulation(name,Q,lastQ,angles,diheds;saveTrajec=true)
    saveChain(string(name,"_Q.csv"),Q)
    saveChain(string(name,"_lastQ.csv"),lastQ)
    if saveTrajec
        saveTrajectory(string(name,"_trajectory.csv"),Q,angles,diheds)
    end
    open(string(name,"_angles-indexes.csv"),"w+") do io
        DelimitedFiles.writedlm(io,diheds,',')
    end
    open(string(name,"_angles-values.csv"),"w+") do io
        DelimitedFiles.writedlm(io,angles,',')
    end
end

function saveLtable(name::AbstractString,ls::AbstractArray)
    open(joinpath(name,"ls.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ls,',')
    end
end

function saveMetaParams(name::AbstractString,parsed_args)
    open(joinpath(name,"metaParams.csv"),"w+") do io
        write(io,"algorithm,$(parsed_args["algorithm"])\n")
        write(io,"temp_init,$(parsed_args["temp_init"])\n")
        write(io,"minFunc,$(parsed_args["minFunc"])\n")
        write(io,"chain,$(parsed_args["chain"])\n")
    end
end

function readChain(name::AbstractString)
    chain = DelimitedFiles.readdlm(name,',',T)
    return PolygonalNew(chain)
end

function readSingleSimulation(name::AbstractString)
    Q = readChain(string(name,"_Q.csv"))
    lastQ = readChain(string(name,"_lastQ.csv"))
    diheds = DelimitedFiles.readdlm(string(name,"_angles-indexes.csv"),',',Int16)
    diheds = reshape(diheds,length(diheds))
    angles = DelimitedFiles.readdlm(string(name,"_angles-values.csv"),',',T)
    angles = reshape(angles,length(angles))
    return Q, lastQ, diheds, angles
end


# `name` must be name of folder
function readMetaParams(name::AbstractString)
    metaParams = DelimitedFiles.readdlm(joinpath(name,"metaParams.csv"),',')
    metaParams = Dict(zip(metaParams[:,1],metaParams[:,2]))
    return metaParams
end


function readLSimulation(name::AbstractString; verbose::Bool=true)
    ls = DelimitedFiles.readdlm(joinpath(name,"ls.csv"))
    ls = reshape(ls,length(ls))
    ln = length(ls)
    metaParams = readMetaParams(name)
    verbose && println("minimizing function")
    verbose && println(metaParams["minFunc"])
    minFunc = getfield(Main,Symbol(metaParams["minFunc"]))
    simuls = [file for file in readdir(name)]
    # saving `_lastQ.csv` files values
    simuls = [file for file in simuls if length(file)>7]
    simuls = [file for file in simuls if file[end-8:end]=="lastQ.csv"]
    # number of total simulations
    lt = length(simuls)
    # filter independent interations
    simuls = [split(file,"_")[1] for file in simuls]
    simuls = unique(simuls)
    simuls = sort(simuls)
    ts_mean = zeros(ln)
    ts_error = zeros(ln)
    minfs_mean = zeros(ln)
    minfs_error = zeros(ln)
    tentativeIters = div(lt,ln,RoundUp)
    verbose && println("Diferent simulations")
    verbose && println(lt)
    verbose && println("L vals")
    verbose && println(ln)
    verbose && println("Assumed iterations")
    verbose && println(tentativeIters)
    ts_table = zeros(ln,tentativeIters)
    minfs_table = zeros(ln,tentativeIters)
    accepted_moves_table = zeros(ln,tentativeIters)
    for (i,sim) in enumerate(simuls)
        iterations = [file for file in readdir(name) if split(file,"_")[1]==sim]
        iterations = [split(file,"_")[2] for file in iterations]
        iterations = unique(iterations)
        ts_sim = zeros(length(iterations))
        minfs_sim = zeros(length(iterations))
        for (j,iter) in enumerate(iterations)
            verbose && println("reading simulation $(sim),$(iter)")
            lastQ = readChain(joinpath(name,string(sim,"_",iter,"_lastQ.csv")))
            diheds = DelimitedFiles.readdlm(joinpath(name,string(sim,"_",iter,"_angles-indexes.csv")),',',Int16)
            diheds = reshape(diheds,length(diheds))
            ts_sim[j] = length(diheds)
            minfs_sim[j] = minFunc(lastQ)
            ts_table[i,j] = length(diheds)
            minfs_table[i,j] = minFunc(lastQ)
            accepted_moves_table[i,j] = length([k for k in diheds if k!=0])
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
    angles = DelimitedFiles.readdlm(string(name,"_angles-values.csv"))
    angles = reshape(angles,length(angles))
    diheds = DelimitedFiles.readdlm(string(name,"_angles-indexes.csv"),',',Int16)
    diheds = reshape(diheds,length(diheds))
    saveTrajectory(string(name,"_trajectory.csv"),Q,angles,diheds)
 end