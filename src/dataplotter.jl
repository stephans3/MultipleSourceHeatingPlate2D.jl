using DelimitedFiles, PlotlyJS
import Dates

# Important: Library Plot.jl shall _not_ be loaded!

timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS");

# Files: file/folder names and paths
basepath        = "results/";
plotsfolder     = "plots/"
input_filename  = "input_history_scenario_";
output_filename = "output_history_scenario_";
temperature2D_filename = "temperature2D_scenario_";


const L = 300.0*10^(-3) # Length of rod
const H = 10.0*10^(-3)  # Height of rod

# Discretization: necessary to import the matrices correctly 
const N₁ = 100; # 50;
const N₂ = 40; #20;
const dx₁ = L/(N₁ - 1);
const dx₂ = H/(N₂ - 1);
const x1span = dx₁/2 : dx₁ : L;
const x2span = (dx₂/2 : dx₂ : H) * 1000;

# Number of scenarios
const N_scenarios = 2;

# Temperature on topside
temperature_on_top = zeros(N₁, N_scenarios);

for i = 1 : N_scenarios

    # Input and Output temporal data
    input_filepath  = basepath*input_filename*string(i)*".txt"
    inputdata       = readdlm(input_filepath, '\t', Float64)
    num_inputsignals = size(inputdata)[2] - 1
    input_averaged   = sum(inputdata[:, 2:end], dims=2) / num_inputsignals;

    output_filepath = basepath*output_filename*string(i)*".txt"
    outputdata      = readdlm(output_filepath, '\t', Float64)
    num_outputsignals   = size(outputdata)[2] - 1
    output_averaged     = sum(outputdata[:, 2:end], dims=2) / num_outputsignals;
    output_variance     = sum((outputdata[:, 2:end] .- output_averaged).^2, dims=2)
    
    input_output_plotname = basepath * plotsfolder  * "input_output_"*string(i)*"_"*string(timestamp) * ".pdf"

    output_plot = Plot(outputdata[:,1], output_averaged, Layout( yaxis = attr(; title="Averaged y", titlefont=attr(; size=16))) )
    input_plot  = Plot(inputdata[:,1], input_averaged, Layout( yaxis = attr(; title="Averaged u", titlefont=attr(; size=16))) )
    output_var_plot = Plot(outputdata[:,1], output_variance, Layout( xaxis = attr(; title="Time [s]", titlefont=attr(; size=20) ) ,yaxis=attr(;title="Variance y", titlefont=attr(; size=16)) ))

    in_out_plot = [output_plot; input_plot; output_var_plot]
    in_out_plot.layout["showlegend"] = false;
    in_out_plot.layout["margin"][:l] = 70;

    PlotlyJS.savefig(in_out_plot, input_output_plotname)
   


    # Temperature distribution: 2 dimensional contour plot 
    temperature2D_filepath  = basepath*temperature2D_filename*string(i)*".txt"
    temperature2D_data      = readdlm(temperature2D_filepath, '\t', Float64 )
    data    = contour(; x=x1span, y=x2span, z=temperature2D_data, colorbar=attr(;title="Temperature [K]", titleside="right", titlefont=attr(;size=20)))
    layout  = Layout(; xaxis = attr(; title="Length [m]", titlefont=attr(; size=20) ) ,yaxis=attr(;title="Heigth [mm]", titlefont=attr(; size=20)) )
    temperature2D_plot = Plot(data, layout)

    temperature2D_plotpath = basepath * plotsfolder * "temperature2D_contour_" * string(i) * "_" *string(timestamp) * ".pdf";
    PlotlyJS.savefig(temperature2D_plot, temperature2D_plotpath)

    temperature_on_top[:,i] = temperature2D_data[:, end];

end 

#Plotting the temperature on top of the plate
trace1 = scatter(; x=x1span, y=temperature_on_top[:,1], mode="lines", name="Scenario 1")
trace2 = scatter(; x=x1span, y=temperature_on_top[:,2], mode="lines", name="Scenario 2")

layout = Layout(; xaxis = attr(; title="Length [m]", titlefont=attr(; size=20), range=[-0.005, 0.301] ) ,yaxis=attr(;title="Temperature [K]", titlefont=attr(; size=20), range=[388, 404]) )
topplot = Plot([trace1, trace2], layout)

topplot_path = basepath * plotsfolder * "top_" * string(timestamp) * ".pdf"

PlotlyJS.savefig(topplot, topplot_path)

