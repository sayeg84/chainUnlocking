include("../src/simulations.jl")
using GLMakie
x = LinRange(0, 10, 100)
y = sin.(x)
colors = repeat([:crimson, :dodgerblue, :slateblue1, :sienna1, :orchid1], 20)

gl_screen = scatter(x, y, color = colors, markersize = 20)
#wait(gl_screen)

#display(current_figure())