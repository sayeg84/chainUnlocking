include("test.jl")
P = PolygonalChain(6)
n = length(P)
ang_idx = rand(n-1:2*n-3)
phi = rand()*(pi/2) + (-pi/4)
inter_flag , newQ = mutateChain(P,ang_idx,phi,debug=true)