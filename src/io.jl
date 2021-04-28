include("types.jl")

using DelimitedFiles, Statistics

function saveChain(name::AbstractString,P::AbstractChain)
    arr = toArray(P)
    open(name,"w+") do io
        DelimitedFiles.writedlm(io,arr,',')
    end
end

function saveTrajectory(name,Q::PolygonalChain2,angles,diheds)
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
    open(string(name,"lvals.csv"),"w+") do io
        DelimitedFiles.writedlm(io,lvals,',')
    end
end

function readChain(name::AbstractString)
    chain = DelimitedFiles.readdlm(name,',',Float64)
    return PolygonalChain2(chain)
end

function readLSimulation(name::AbstractString)
    ls = DelimitedFiles.readdlm(string(name,"lvals.csv"))
    ls = reshape(ls,length(ls))
    simuls = [file for file in readdir(name)]
    simuls = [file for file in simuls if file[end-4:end]=="P.csv"]
    simuls = [split(file,"_")[1] for file in simuls]
    simuls = unique(simuls)
    simuls = sort(simuls)
    ts_mean = zeros(length(ls))
    ts_error = zeros(length(ls))
    rmsds_mean = zeros(length(ls))
    rmsds_error = zeros(length(ls))
    for (i,sim) in enumerate(simuls)
        iterations = [file for file in readdir(name) if split(file,"_")[1]==sim]
        iterations = [split(file,"_")[2] for file in iterations]
        iterations = unique(iterations)
        sim_ts = zeros(length(iterations))
        sim_rmsds = zeros(length(iterations))
        for (j,iter) in enumerate(iterations)
            println("reading simulation $(sim),$(iter)")
            P = readChain(string(name,sim,"_",iter,"_P.csv"))
            lastQ = readChain(string(name,sim,"_",iter,"_lastQ.csv"))
            diheds = DelimitedFiles.readdlm(string(name,sim,"_",iter,"_angles-values.csv"))
            diheds = reshape(diheds,length(diheds))
            sim_ts[j] = length(diheds)
            sim_rmsds[j] = overlapedRmsd(lastQ,P)
        end
        ts_mean[i] = Statistics.mean(sim_ts)
        ts_error[i] = Statistics.std(sim_ts)
        rmsds_mean[i] = Statistics.mean(sim_rmsds)
        rmsds_error[i] = Statistics.std(sim_rmsds)
    end
    return ls,ts_mean,ts_error,rmsds_mean,rmsds_error
end

