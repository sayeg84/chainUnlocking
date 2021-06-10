include("types.jl")

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

function saveSimulation(name,P,Q,lastQ,angles,diheds;saveTrajec=true)
    saveChain(string(name,"_P.csv"),P)
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

function saveLtable(name::AbstractString,lvals::AbstractArray)
    open(joinpath(name,"lvals.csv"),"w+") do io
        DelimitedFiles.writedlm(io,lvals,',')
    end
end

function readChain(name::AbstractString)
    chain = DelimitedFiles.readdlm(name,',',Float64)
    return PolygonalNew(chain)
end

function readSingleSimulation(name::AbstractString)
    P = readChain(string(name,"_P.csv"))
    Q = readChain(string(name,"_Q.csv"))
    lastQ = readChain(string(name,"_lastQ.csv"))
    diheds = DelimitedFiles.readdlm(string(name,"_angles-indexes.csv"),',',Int16)
    diheds = reshape(diheds,length(diheds))
    angles = DelimitedFiles.readdlm(string(name,"_angles-values.csv"),',',Float64)
    angles = reshape(angles,length(angles))
    return P,Q,lastQ, diheds, angles
end


function readLSimulation(name::AbstractString)
    ls = DelimitedFiles.readdlm(joinpath(name,"lvals.csv"))
    ls = reshape(ls,length(ls))
    ln = length(ls)
    simuls = [file for file in readdir(name)]
    simuls = [file for file in simuls if file[end-4:end]=="P.csv"]
    # number of total simulations
    lt = length(simuls)
    # filter independent interations
    simuls = [split(file,"_")[1] for file in simuls]
    simuls = unique(simuls)
    simuls = sort(simuls)
    
    ts_mean = zeros(ln)
    ts_error = zeros(ln)
    rmsds_mean = zeros(ln)
    rmsds_error = zeros(ln)
    
    tentativeIters = div(lt,ln,RoundUp)
    ts_table = zeros(ln,tentativeIters)
    rmsds_table = zeros(ln,tentativeIters)
    accepted_moves_table = zeros(ln,tentativeIters)
    for (i,sim) in enumerate(simuls)
        iterations = [file for file in readdir(name) if split(file,"_")[1]==sim]
        iterations = [split(file,"_")[2] for file in iterations]
        iterations = unique(iterations)
        ts_sim = zeros(length(iterations))
        rmsds_sim = zeros(length(iterations))
        for (j,iter) in enumerate(iterations)
            println("reading simulation $(sim),$(iter)")
            P = readChain(joinpath(name,string(sim,"_",iter,"_P.csv")))
            lastQ = readChain(joinpath(name,string(sim,"_",iter,"_lastQ.csv")))
            diheds = DelimitedFiles.readdlm(joinpath(name,string(sim,"_",iter,"_angles-indexes.csv")),',',Int16)
            diheds = reshape(diheds,length(diheds))
            ts_sim[j] = length(diheds)
            rmsds_sim[j] = overlapedRmsd(lastQ,P)
            ts_table[i,j] = length(diheds)
            rmsds_table[i,j] = overlapedRmsd(lastQ,P)
            accepted_moves_table[i,j] = length([k for k in diheds if k!=0])
        end
        ts_mean[i] = Statistics.mean(ts_sim)
        ts_error[i] = Statistics.std(ts_sim)
        rmsds_mean[i] = Statistics.mean(rmsds_sim)
        rmsds_error[i] = Statistics.std(rmsds_sim)
    end
    return ls,ts_mean,ts_error,rmsds_mean,rmsds_error, ts_table, rmsds_table,accepted_moves_table
end




function makeCoordinates(name)
    Q = readChain(string(name,"_Q.csv"))
    angles = DelimitedFiles.readdlm(string(name,"_angles-values.csv"))
    angles = reshape(angles,length(angles))
    diheds = DelimitedFiles.readdlm(string(name,"_angles-indexes.csv"),',',Int16)
    diheds = reshape(diheds,length(diheds))
    saveTrajectory(string(name,"_trajectory.csv"),Q,angles,diheds)
 end