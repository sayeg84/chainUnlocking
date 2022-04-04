include("algorithms.jl")

abstract type AbstractAlgorithm end

abstract type MHAlgorithm <: AbstractAlgorithm end

struct Annealing <: MHAlgorithm end

struct Genetic <: MHAlgorithm end

abstract type GDAlgorithm <: AbstractAlgorithm end

struct Sobolev <: GDAlgorithm end

struct CurvifiedSobolev <: GDAlgorithm end


function runSingle(P,SimulType::Annealing,savename,parsed_args)
    tempProgram = getfield(Main,Symbol(parsed_args["temp_program"]))()
    minFunc = getfield(Main,Symbol(parsed_args["minFunc"]))
    initQs,finalQs,fun_vals,ang_idxs,ang_vals = multipleSimulatedAnnealing(P,
    minFunc,
    parsed_args["tolerance"],
    parsed_args["max_angle"],
    -1*parsed_args["max_angle"],
    parsed_args["internal"];
    temp_init=parsed_args["temp_init"],
    temp_f=parsed_args["temp_f"],
    iter_per_temp=parsed_args["iter_per_temp"],
    tempProgram=tempProgram,
    population=parsed_args["indep_simuls"],
    max_iter=parsed_args["max_iter"],
    mut_k=parsed_args["mut_k"],
    debug=parsed_args["debug"])
    if !isempty(savename)
        open(string(savename,"fun_vals.csv"),"w+") do io
            DelimitedFiles.writedlm(io,fun_vals,',')
        end
        open(string(savename,"ang_idxs.csv"),"w+") do io
            DelimitedFiles.writedlm(io,ang_idxs,',')
        end
        open(string(savename,"ang_vals.csv"),"w+") do io
            DelimitedFiles.writedlm(io,ang_vals,',')
        end
        saveChains(string(savename,"init_Qs.csv"),initQs)
        saveChains(string(savename,"final_Qs.csv"),finalQs)
    end
    return initQs,finalQs,fun_vals,ang_idxs,ang_vals
end

function runSingle(P::PolygonalChain,SimulType::Genetic,savename::AbstractString,parsed_args)
    selection = getfield(Main,Symbol(parsed_args["selection"]))()
    minFunc = getfield(Main,Symbol(parsed_args["minFunc"]))
    initQs,finalQs,fun_vals,ang_idxs,ang_vals,parents = genetic(P,
    minFunc,
    parsed_args["tolerance"],
    parsed_args["max_angle"],
    -1*parsed_args["max_angle"],
    parsed_args["internal"];
    selection=selection,
    population=parsed_args["indep_simuls"],
    max_iter=parsed_args["max_iter"],
    mut_k=parsed_args["mut_k"],
    debug=parsed_args["debug"])
    if !isempty(savename)
        open(string(savename,"fun_vals.csv"),"w+") do io
            DelimitedFiles.writedlm(io,fun_vals,',')
        end
        open(string(savename,"ang_idxs.csv"),"w+") do io
            DelimitedFiles.writedlm(io,ang_idxs,',')
        end
        open(string(savename,"ang_vals.csv"),"w+") do io
            DelimitedFiles.writedlm(io,ang_vals,',')
        end
        open(string(savename,"parents.csv"),"w+") do io
            DelimitedFiles.writedlm(io,parents,',')
        end
        saveChains(string(savename,"init_Qs.csv"),initQs)
        saveChains(string(savename,"final_Qs.csv"),finalQs)
    end
    return initQs,finalQs,fun_vals,ang_idxs,ang_vals
end

function runSingle(P::PolygonalChain,SimulType::Sobolev,savename::AbstractString,parsed_args)
    P = PolygonalChain([Point(Float64,p) for p in P.vertices])
    trajectory = sobolevGradientDescent(P,
    parsed_args["max_iter"],
    tau=parsed_args["time_step"],
    debug=parsed_args["debug"])
    if !isempty(savename)
        saveTrajectory(string(savename,"trajectory.csv"),P,trajectory)
        saveChain(string(savename,"init_Q.csv"),P)
    end
    finalQ = PolygonalChain(trajectory[end,:])
    return P,finalQ
end

function runSingle(P::PolygonalChain,SimulType::CurvifiedSobolev,savename::AbstractString,parsed_args)
    if !isempty(parsed_args["ndens"])
        ndens = [parse(Int64,s) for s in split(parsed_args["ndens"],",")]
        if length(ndens) != length(P)
            error("ndens $(parsed_args["ndens"]) doesn't have $(length(P)) integers")
        end
    else
        error("Argument `ndens` is required")
    end
    P = PolygonalChain([Point(Float64,p) for p in P.vertices])
    Q = helixify(P,ndens,r=5e-3)
    trajectory = sobolevGradientDescent(Q,ndens,
    parsed_args["max_iter"],
    tau=parsed_args["time_step"],
    debug=parsed_args["debug"])
    if !isempty(savename)
        saveTrajectory(string(savename,"trajectory.csv"),P,trajectory)
        saveChain(string(savename,"init_Q.csv"),P)
    end
    finalQ = PolygonalChain(trajectory[end,:])
    return Q,finalQ
end


