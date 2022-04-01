include("io.jl")
include("argumentParsing.jl")
function main()
    parsed_args = parsed_args = parse_commandline()
    P = parametricCurveChain(treefoil,30,0,deg2rad(340))    
    #P = parametricCurveChain(polynomialTrefoil,40,-3,3)    
    #P = knittingNeedle(2.23)
    # chains must be in Float64 type, otherwise errors arise when inverting matrices
    P = PolygonalChain([Point(Float64,p) for p in P.vertices])
    res = sobolevGradientDescent(P,parsed_args["max_iter"],tau=parsed_args["time_step"],debug=parsed_args["debug"])
    if !isdir(parsed_args["path"])
        mkdir(parsed_args["path"])
    end
    name = joinpath(parsed_args["path"],"sobolev.csv")
    saveTrajectory(name,P,res)
end
main()