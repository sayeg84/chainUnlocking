using ArgParse
using Dates
const date_str = Dates.format(Dates.now(),"YYYY-mm-dd_HH-MM")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--algorithm"
            help = "Algorithm to make simulation"
            arg_type = String
            default = "localRandomSearch"
        "--minFunc"
            help = "Minimizing function"
            arg_type = String
            default = "distToFlat"
        "--chain"
            help = "Chain to simulate"
            arg_type = String
            default = "knittingNeedle"
        "--path"
            help = "Folder to save the simulations"
            arg_type = String
            default = joinpath("..","outputs",date_str)
        "--indep_simuls"
            help = "Number of simulations per l value"
            arg_type = Int
            default = 4
        "--processes"
            help = "Parallel processes to repart the simulations"
            arg_type = Int
            default = 4
        "--lmin"
            help = "minimum value for l"
            arg_type = Float64
            default = 1.4
        "--lmax"
            help = "maximum value for l"
            arg_type = Float64
            default = 1.42
        "--lvals"
            help = "Values for the L interval"
            arg_type = Int
            default = 12
        "--log_l"
            help = "flag to tell if logarithmic space must be filled for the l Values"
            action = :store_true
        "--max_angle"
            help = "maximum angle for rotation"
            arg_type = Float64
            default = pi/40
        "--max_iter"
            help = "maximum number of iterations"
            arg_type = Int
            default = 10_000 
        "--temp_init"
            help = "Initial temperature for simulated annelaing"
            arg_type = Float64
            default = 4.0
        "--temp_f"
            help = "Final temperature for simulated annelaing"
            arg_type = Float64
            default = 1e-2
        "--iter_per_temp"
            help = "Iterations per temperature for simulated annealing"
            arg_type = Int
            default = 20
        "--tolerance"
            help = "Tolerance for minimum value of function"
            arg_type = Float64
            default = 1e-2
        "--internal"
            help= "Flag to indicate if internal angles are also allowed to be changed"
            action = :store_true
        "--debug"
            help= "Flag to indicate if internal angles are also allowed to be changed"
            action = :store_true
    end
    return parse_args(s)
end
