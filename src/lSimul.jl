include("io.jl")

function main()    
    ls,ts_mean,ts_error,rmsds_mean,rmsds_error = readLSimulation(ARGS[1])
    println("Making Plots")
    scatter(ls,ts_mean,yerror=ts_error,label=false,title="L analysis",ylabel="steps",xlabel="L")
    savefig(string(ARGS[1],"l-t.pdf"))
    scatter(ls,rmsds_mean,yerror=rmsds_error,label=false,title="L analysis",ylabel="rmsd",xlabel="L")
    savefig(string(ARGS[1],"l-rmsd.pdf"))
    println("Done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    main()
end