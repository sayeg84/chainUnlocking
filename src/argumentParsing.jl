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
            default = "tangentEnergyFrac"
        "--chain"
            help = "Chain to simulate. Must be function of one real parameter that returns chain or path to csv file with coordinates. It will take last coordinates as chain."
            arg_type = String
            default = "knittingNeedle"
        "--path"
            help = "Folder to save the simulations. Defaults to current date"
            arg_type = String
            default = joinpath("..","outputs",date_str)
        "--indep_simuls"
            help = "Number of independent simulations."
            arg_type = Int
            default = 4
        "--processes"
            help = "Parallel processes to distribute the simulations"
            arg_type = Int
            default = 1
        "--max_iter"
            help = "Maximum number of iterations"
            arg_type = Int
            default = 10_000
        "--tolerance"
            help = "Tolerance for minimum value of function"
            arg_type = Float64
            default = -1000.0
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
        "--internal"
            help= "Flag to indicate if internal angles are also allowed to be changed. Only for MH simulations"
            action = :store_true
        "--max_angle"
            help = "Maximum angle change for DOFs Only for MH simulations"
            arg_type = Float64
            default = pi/40
        "--mut_k"
            help = "Maximum degree of mutation. Only for MH simulations"
            arg_type = Int
            default = 3 
        "--selection"
            help = "Maximum degree of mutation. Only for Genetic simulations"
            arg_type = String
            default = "RouletteWheelSelection" 
        "--selection_k"
            help = "Parameter for selection method, must be smaller than `indep_simuls`. Only for Genetic simulations"
            arg_type = Int
            default = 3 
        "--temp_init"
            help = "Initial temperature. Only valid for Simulated Annealing"
            arg_type = Float64
            default = 4.0
        "--temp_f"
            help = "Final temperature. Only valid for Simulated Annealing"
            arg_type = Float64
            default = 1e-2
        "--iter_per_temp"
            help = "Iterations per temperature. Only valid for Simulated Annealing"
            arg_type = Int
            default = 20
        "--temp_program"
            help = "Temperature program. Only valid for Simulated Annealing"
            arg_type = String
            default = "ExponentialProgram"
        "--time_step"
            help = "Time step size. Only valid for GD simulations"
            arg_type = Float64
            default = 1e-2
        "--ndens"
            help = "Number of new nodes on each link represented as ints separated by commas (\"2,3,4\"). Only valid for CurvifiedSobolev simulation"
            arg_type = String
            default = ""
        "--retake_curve"
            help = "Flag to indicate that curve is already helixified and doesn't need to be done again. Useful for running curvified simulations starting with an already existing chain"
            action = :store_true
        "--debug"
            help= "Flag passed to methods to print debug info on screen"
            action = :store_true
    end
    return parse_args(s)
end
