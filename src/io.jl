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
    return funcvals, newQ
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

function readSingleSimulation(name::AbstractString,simul::GDAlgorithm,minFunc)
    source,_ = DelimitedFiles.readdlm(string(name,"_trajectory.csv"),',',header=true)
    funcvals = zeros(size(source,1))
    for i in 1:size(source,1)
        Q = PolygonalChain(source[i,:]) 
        funcvals[i] = minFunc(Q)
    end
    return [PolygonalChain(source[end,:]) ],hcat(funcvals),[size(source,1)]
end

function readSingleSimulation(name::AbstractString,simul::MHAlgorithm,minFunc)
    sources,_ = DelimitedFiles.readdlm(string(name,"_init_Qs.csv"),',',header=true)
    Qs = [PolygonalChain(sources[i,:]) for i in 1:size(sources,1)]
    ang_idxs   = DelimitedFiles.readdlm(string(name,"_ang_idxs.csv"),',',Int16)
    ang_vals   = DelimitedFiles.readdlm(string(name,"_ang_vals.csv"),',',BigFloat)
    m = length(Qs)
    funcvals = zeros(size(ang_idxs))
    accepted_moves = zeros(m)
    for j in 1:m
        println("Reading iteration = $j")
        aux,_ = funcValues(Qs[j],ang_idxs[:,j],ang_vals[:,j],minFunc)
        funcvals[:,j] = aux
        accepted_moves[j] = length([k for k in ang_idxs[:,j] if k!=0])
    end
    sources,_ = DelimitedFiles.readdlm(string(name,"_final_Qs.csv"),',',header=true)
    Qs = [PolygonalChain(sources[i,:]) for i in 1:m]
    #display(toArr
    return  Qs, funcvals, accepted_moves
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
    simul   = getfield(Main,Symbol(metaParams["algorithm"]))()
    minFunc = getfield(Main,Symbol(metaParams["minFunc"]))
    indep_simuls = typeof(simul) <: MHAlgorithm ? metaParams["indep_simuls"] : 1
    # saving `_lastQ.csv` files values
    ts_table             = zeros(ln,indep_simuls)
    minfs_table          = zeros(ln,indep_simuls)
    accepted_moves_table = zeros(ln,indep_simuls)
    for i in 1:ln
        println("Reading lval = $i")
        n1zeros = Int(ceil(log10(ln+1)))
        n1 = lpad(i,n1zeros,'0')
        Qs,funvals,accepted = readSingleSimulation(joinpath(name,n1),simul,minFunc)
        if burnout < 1
            ncut             = Int(ceil(burnout*length(ang_idxs)))
            minfs_table[i,:] = Statistics.mean(fun_vals[ncut:end,:],dims=1)
        else
            minfs_table[i,:] = funvals[end,:]
        end
        accepted_moves_table[i,:] = accepted
        ts_table[i,:] = [length(funvals[:,j]) for j in 1:size(funvals,2)]
    end
    return ls, ts_table, minfs_table, accepted_moves_table
end



function makeCoordinates(name::AbstractString,i::Integer)
    sources,_ = DelimitedFiles.readdlm(string(name,"init_Qs.csv"),',',header=true)
    Qs = [PolygonalChain(sources[i,:]) for i in 1:size(sources,1)]
    ang_idxs   = DelimitedFiles.readdlm(string(name,"ang_idxs.csv"),',',Int16)
    ang_vals   = DelimitedFiles.readdlm(string(name,"ang_vals.csv"),',',BigFloat)
    saveTrajectory(string(name,"_",i,"_trajectory.csv"),Qs[i],ang_vals[:,i],ang_idxs[:,i])
 end

 function generateDrawingData(n::Integer)
    P = curveToChain(polynomialTrefoilCurve,n,-2,2)
    open("../../chain.csv","w+") do io; 
        writedlm(io,to2DArray(P),',') 
    end
    grad = deriv(tangentEnergyFrac,P)
    grad = reshape(grad,(3,n))'
    open("../../normalgrad.csv","w+") do io; 
        writedlm(io,grad,',') 
    end
    gradf = sobolevGradient(P,true)
    gradf = reshape(gradf[1:3*n],(3,n))'
    open("../../sobgrad.csv","w+") do io; 
        writedlm(io,gradf,',') 
    end
end