include("test.jl")
P = PolygonalChain(6)
n = length(P)
dihed = rand(n-1:2*n-3)
phi = rand()*(pi/2) + (-pi/4)
inter_flag , newQ = mutateChain(P,dihed,phi,debug=true)