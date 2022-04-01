include("io.jl")
include("argumentParsing.jl")
function main()
    parsed_args = parsed_args = parse_commandline()
    P = parametricCurveChain(treefoil,30,0,deg2rad(340))    
    minFunc  =  getfield(Main,Symbol(parsed_args["minFunc"]))
    #chainFunc = getfield(Main,Symbol(parsed_args["chain"]))
    #P = parametricCurveChain(polynomialTrefoil,40,-3,3)    
    P = knittingNeedle(2.23)
    # chains must be in BigFloat type, otherwise errors arise when checking intersections
    P = PolygonalChain([Point(BigFloat,p) for p in P.vertices])
    Qs,fun_vals = genetic(P,
    minFunc,
    parsed_args["tolerance"],
    parsed_args["max_angle"],
    -1*parsed_args["max_angle"],
    parsed_args["internal"],;
    max_iter=parsed_args["max_iter"],
    debug=parsed_args["debug"])
    if !isdir(parsed_args["path"])
        mkdir(parsed_args["path"])
    end
    open(joinpath(parsed_args["path"],"fun_vals.csv"),"w+") do io
        DelimitedFiles.writedlm(io,fun_vals,',')
    end
    saveChains(joinpath(parsed_args["path"],"Qs.csv"),Qs)
end
main()