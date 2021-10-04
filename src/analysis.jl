include("algorithms.jl")
include("io.jl")

using Plots, DelimitedFiles, ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--path"
            help = "Folder to save the simulations"
            arg_type = String
            required = true
        "--burnout"
            help = "Fraction of burnout values. Value of 1 saves only the last one"
            arg_type = Float64
            default  = 1.0
    end
    return parse_args(s)
end

const parsed_args = parse_commandline()

function main()    
    ls,ts_mean,ts_error,minfs_mean,minfs_error,ts_table,minfs_table,accepted_moves_table = readLSimulation(parsed_args["path"],parsed_args["burnout"])
    println("saving results")
    open(joinpath(parsed_args["path"],"results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,minfs_mean,minfs_error)
        write(io,"l,t_mean,t_std,minf_mean,minf_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    open(joinpath(parsed_args["path"],"ts_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ts_table,',')
    end
    open(joinpath(parsed_args["path"],"minfs_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,minfs_table,',')
    end
    open(joinpath(parsed_args["path"],"accepted_moves_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,accepted_moves_table,',')
    end
    
    println("Making Plots")
    scatter(ls,ts_mean,yerror=ts_error./2,label=false,title="L analysis",ylabel="steps",xlabel="L")
    savefig(joinpath(parsed_args["path"],"l-t.pdf"))
    scatter(ls,minfs_mean,yerror=minfs_error./2,label=false,title="L analysis",ylabel="minf",xlabel="L")
    savefig(joinpath(parsed_args["path"],"l-minf.pdf"))
    println("Done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end