include("io.jl")
using Plots, DelimitedFiles

function main()    
    ls,ts_mean,ts_error,rmsds_mean,rmsds_error = readLSimulation(ARGS[1])
    println("saving results")
    open(string(ARGS[1],"results.csv"),"w+") do io
        table = hcat(ls,ts_mean,ts_error,rmsds_mean,rmsds_error)
        write(io,"l,t_mean,t_std,rmsd_mean,rmsd_std\n")
        DelimitedFiles.writedlm(io,table,',')
    end
    println("Making Plots")
    scatter(ls,ts_mean,yerror=ts_error,label=false,title="L analysis",ylabel="steps",xlabel="L")
    savefig(string(ARGS[1],"l-t.pdf"))
    scatter(ls,rmsds_mean,yerror=rmsds_error,label=false,title="L analysis",ylabel="rmsd",xlabel="L")
    savefig(string(ARGS[1],"l-rmsd.pdf"))
    println("Done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end