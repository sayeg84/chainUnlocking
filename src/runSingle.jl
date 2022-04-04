include("argumentParsing.jl")
include("io.jl")

function main()
    parsed_args = parse_commandline()
    SimulType = getfield(Main,Symbol(parsed_args["algorithm"]))()
    minFunc = getfield(Main,Symbol(parsed_args["minFunc"]))
    chainFunc = getfield(Main,Symbol(parsed_args["chain"]))
    if !isdir(parsed_args["path"])
        mkdir(parsed_args["path"])
    end
    P = chainFunc(parsed_args["lmin"])
    runSingle(P,SimulType,parsed_args["path"],parsed_args)
end

println("Single simulator")
println()
main()

