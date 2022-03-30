include("io.jl")

using Plots


function makePlot(vals,name)
    fig = Plots.plot(size=(1000,400))
    Plots.plot!(1:length(vals),vals)
    Plots.xlabel!("Sim step")
    Plots.ylabel!("Func value")
    Plots.savefig(string(name,"_min.pdf"))
    Plots.closeall()
end

function specialPlot(vals1,vals2,name)
    fig = Plots.plot(size=(1000,400))
    Plots.plot!(1:length(vals1),vals1,label="recalc")
    Plots.plot!(1:length(vals2),vals2,label="simulated")
    Plots.plot!(1:length(vals2),vals1-vals2,label="diff")
    Plots.xlabel!("Sim step")
    Plots.ylabel!("Func value")
    Plots.savefig(string(name,"_diffs.pdf"))
    Plots.closeall()
end

function main()
    metaParams = readMetaParams(dirname(ARGS[1]))
    minFunc = getfield(Main,Symbol(metaParams["minFunc"]))
    Q, lastQ, ang_idxs, ang_vals = readSingleSimulation(ARGS[1])
    saveTrajectory(string(ARGS[1],"_trajectory.csv"),Q,ang_vals,ang_idxs)
    vals = funcValues(Q,ang_idxs,ang_vals,minFunc)
    simvals = DelimitedFiles.readdlm(string(ARGS[1],"_fun_vals"),',')
    println("Beggining")
    println(vals[begin])
    println(minFunc(Q))
    println("End")
    println(vals[end])
    println(minFunc(lastQ))
    makePlot(vals,ARGS[1])
    specialPlot(vals,simvals,ARGS[1])
    open(string(ARGS[1],"_diffs.csv"),"w+") do io
        DelimitedFiles.writedlm(io,vals-simvals,',')
    end
end
main()