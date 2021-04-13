include("types.jl")

using DelimitedFiles

function saveChain(name::AbstractString,P::AbstractChain)
    arr = toArray(P)
    open(name,"w+") do io
        DelimitedFiles.writedlm(io,arr,',')
    end
end

function saveTrajectory(name,Q::PolygonalChain2,angles_values,angles_indexes)
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
        arr = toArray(newQ)
        coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
        coordinates = join(coordinates,",")
        coordinates = string(coordinates,"\n")
        write(io,coordinates)
        for i in 1:length(angles_values)
            if angles_indexes[i]!= 0 
                dihedralRotate!(newQ,angles_indexes[i],angles_values[i])
                arr = toArray(newQ)
                coordinates = [string(x) for x in reshape(transpose(arr),length(arr))]
                coordinates = join(coordinates,",")
                coordinates = string(coordinates,"\n")
                write(io,coordinates)
            end
        end
    end
end

function saveSimulation(name,P,Q,angles_values,angles_indexes)
    saveChain(string(name,"_P.csv"),P)
    saveChain(string(name,"_Q.csv"),Q)
    saveTrajectory(string(name,"_trajectory.csv"),Q,angles_values,angles_indexes)
    open(string(name,"_angles-indexes.csv"),"w+") do io
        DelimitedFiles.writedlm(io,angles_indexes,',')
    end
    open(string(name,"_angles-values.csv"),"w+") do io
        DelimitedFiles.writedlm(io,angles_values,',')
    end
end