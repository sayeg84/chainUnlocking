include("io.jl")
using Plots, DelimitedFiles

function main()    
    ls,ts_mean,ts_error,rmsds_mean,rmsds_error,ts_table,rmsds_table,accepted_moves_table = readLSimulation(ARGS[1])
    println("saving results")
    open(joinpath(ARGS[1],"results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,rmsds_mean,rmsds_error)
        write(io,"l,t_mean,t_std,rmsd_mean,rmsd_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    open(joinpath(ARGS[1],"ts_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,ts_table,',')
    end
    open(joinpath(ARGS[1],"rmsds_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,rmsds_table,',')
    end
    open(joinpath(ARGS[1],"accepted_moves_table.csv"),"w+") do io
        DelimitedFiles.writedlm(io,accepted_moves_table,',')
    end
    
    println("Making Plots")
    scatter(ls,ts_mean,yerror=ts_error./2,label=false,title="L analysis",ylabel="steps",xlabel="L")
    savefig(joinpath(ARGS[1],"l-t.pdf"))
    scatter(ls,rmsds_mean,yerror=rmsds_error./2,label=false,title="L analysis",ylabel="rmsd",xlabel="L")
    savefig(joinpath(ARGS[1],"l-rmsd.pdf"))
    println("Done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end