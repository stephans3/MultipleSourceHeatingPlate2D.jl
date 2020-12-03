using LinearAlgebra, SparseArrays, DifferentialEquations

#=
    Scenario 1: Nominal actuators
    Scenario 2: Realistic actuators
=#
const scenario = 1;


# Configuration of physical properties λ and c
function property(temperature :: T , a₀ :: T, a₁ :: T) where T <: Real
    Θ = temperature

    return a₀ + a₁ * Θ
end

# Emitted heat flux on boundary
function boundary_heatflux(boundary_temperature :: T, ambient_temperature ::T, transfer_coefficient :: T, radiation_coefficient :: T) where T <: Real
    θbound = boundary_temperature;
    θamb = ambient_temperature
    h = transfer_coefficient
    k = radiation_coefficient

    return -1*h *(θbound - θamb) -  k *(θbound^4 - θamb^4)
end

# Emitted heat flux on boundary
function boundary_heatflux(boundary_temperature :: Array{T}, ambient_temperature :: Array{T}, transfer_coefficient :: T, radiation_coefficient :: T) where T <: Real
    θbound = boundary_temperature;
    θamb = ambient_temperature
    
    if size(θbound) != size(θamb)
        error("Boundary temperature array and ambient temperature array have different sizes: size(boundary_temperature) = $(size(θbound)) ≠ $(size(θamb)) = size(ambient_temperature)!")
    end
    
    h = transfer_coefficient
    k = radiation_coefficient

    return -1*h *(θbound .- θamb) -  k *(θbound.^4 .- θamb.^4)
end


# Spatial approximation
const L = 300.0*10^(-3) # Length of rod
const H = 10.0*10^(-3)  # Height of rod

const N₁ = 100; # Number of cells in x₁ direction
const N₂ = 40;  # Number of cells in x₂ direction
const Δx₁ = L/N₁; # Discretizaion
const Δx₂ = H/N₂;
const x1span = Δx₁/2 : Δx₁ : L; # Coordinates of cell centers
const x2span = Δx₂/2 : Δx₂ : H;

# Heat transer + radiation
const h = 10.0             # Heat transfer coefficient
const ϵ = 0.6;             # Emissivity
const sb = 5.67*10^(-8);   # Stefan-Boltzmann constant
const k = ϵ*sb;            # Heat radiation coefficient

const ρ = 7800.0;       # Density
const ρinv = inv(ρ);
const Id12 = sparse(I, N₁*N₂, N₁*N₂) # Template matrix for specific heat capacity

const θamb = 300.0; # Ambient temperature


# Initial noisy data
const a₀ = 3;   # Scaling
const a₁ = 10;
const a₂ = 5;
const noise = a₀ * cos.(a₁*2*pi/L * x1span) .* transpose(cos.(a₂ *2*pi/H * x2span))
θinit = ones(N₁, N₂) * 300.0 .+ noise;
θinit = reshape(θinit, N₁*N₂, 1 )


# Input / Output
const Nu = 5; # Number of inputs
const Ny = 5; # Number of outputs

b_in = spzeros(N₁, Nu)          # Actuators' spatial characterization
g_out = spzeros(N₁, Ny)         # Sensors' spatial characterization

const β_in = spzeros(Int64,N₁, Nu)    # Input partition
const γ_out = spzeros(Int64, N₁, Ny   )# Output partition

const xspan = Δx₁/2  : Δx₁ : L;           # Centers of cells
const xc_in = (L/Nu)*((1/2) : 1 : Nu)     # Centers of actuators
const xc_out = (L/Ny)*((1/2) : 1 : Ny)    # Centers of sensors


# Decision of scenario: actuators' spatial characterization
if scenario == 1
    in_coeff = (m = 1.0, M = 0.0, ν = 4.0 )
else
    in_coeff = (m = 1.0, M = 30.0, ν = 4.0 )
end

# Sensors' spatial characterization
out_coeff = (m = 1.0, M = 10.0, ν = 4.0 )

function characterization(x, xcenter, m::Number=1.0, M::Number=1.0, ν::Number=1.0)
    return  m*exp(-1 * norm( M * (x - xcenter) )^ν )
end

# Input characterization
for i = 1 : Nu
    β_in[round(Int, (i-1)*N₁/Nu)+1:round(Int, i*N₁/Nu), i] .= 1;
    b_in[:,i] = characterization.(xspan, xc_in[i], in_coeff.m, in_coeff.M , in_coeff.ν)  
end

b_in = β_in .* b_in;

# Output characterization
for i = 1 : Ny
    γ_out[round(Int, (i-1)*N₁/Ny)+1:round(Int, i*N₁/Ny), i] .= 1;
    g_out[:,i] = characterization.(xspan, xc_in[i], out_coeff.m, out_coeff.M , out_coeff.ν)  
end
  
g_out = γ_out .* g_out;
const g_weight = inv.(sum(transpose(g_out), dims=2))


# Controller
const Kp = ones(Nu) * 10^(4);     # Controller gain
const u_in  = zeros(Nu)           # System input signals
const y_out = zeros(Ny)           # System output signals
const y_ref = ones(Ny) * 400.0;   # Reference temperatures


# History of input signals
global u_hist
u_hist = zeros(Nu,0);

# History of measurements
global y_hist
y_hist = zeros(Ny,0)

# History of time steps
const t_hist = Real[];


# System dynamics: approximated spatial derivatives
const Mx1_lo = spdiagm(-1 => ones(N₁*N₂-1), 0 => -1*ones(N₁*N₂))
const Mx1_hi = spdiagm( 0 => -1*ones(N₁*N₂), 1 => ones(N₁*N₂-1)) 

const Mx2_lo = spdiagm(-1*N₁ => ones( N₁*(N₂-1) ), 0 => -1*ones(N₁*N₂))
const Mx2_hi = spdiagm( 0 => -1*ones(N₁*N₂), N₁ => ones( N₁*(N₂-1) )) 


# Boundary value description
const bound_under = sparse(I,N₁*N₂, N₁)
const bound_top = spzeros(N₁*N₂, N₁)
bound_top[(N₂-1)*N₁+1 : N₁*N₂,1:N₁] = sparse(I, N₁, N₁)

const bound_left = spzeros(N₁*N₂, N₂)
const bound_right = spzeros(N₁*N₂, N₂)

for k = 1 : N₂
    bound_left[1 + (k-1)*N₁, k] = 1
    bound_right[k*N₁, k] = 1
end


# Approximated PDE as large-scale ODE
function heat_eq(dθ,θ,p,t)
    global y_hist
    global u_hist
    push!(t_hist, t);
  
    # Reshape temperature vector to matrix (2-dimensional temperature field)
    θ = reshape(θ, N₁, N₂)

    # Computation of system output  signals
    y_out = g_weight .* (transpose(g_out) * θ[:,end]);
    
    # Controller
    u = Kp .* (y_ref - y_out)
    u = u .* (0.5*(sign.(u) .+ 1.0));
  
    # Saving input and output signals
    y_hist = hcat(y_hist, Matrix(y_out))
    u_hist = hcat(u_hist, u)

    # Averaging temperatures at internal cell boundaries
    θ1m = (θ[1:end-1,:] + θ[2:end,:])/2.0;          # x₁-direction
    θ2m = (θ[:, 1:end-1] + θ[:, 2:end])/2.0;        # x₂-direction
    
    # Specific heat capacity
    c = property.(θ, 330.0, 0.4)
    Cinv = reshape(inv.(c), N₁* N₂,1) .* Id12 

    # Thermal conductivity
    λ₁ = property.( θ1m , 10., 0.1); # x₁-direction
    λ₂ = property.( θ2m , 10., 0.1); # x₂-direction

    # Add addional "virtual" cell at left and right boundary
    λ₁_ext1 = cat( zeros(1,N₂), λ₁, dims=1)
    λ₁_ext2 = cat( λ₁, zeros(1,N₂), dims=1)
    
    # Reshape λ matrix of 2-dimensional temperature field to vector
    λ₁_ext1 = reshape(λ₁_ext1, N₁*N₂, 1 )
    λ₁_ext2 = reshape(λ₁_ext2, N₁*N₂, 1 )
    
    # Add addional "virtual" cell at underside and topside boundary
    λ₂_ext1 = cat( zeros(N₁,1), λ₂, dims=2 ) 
    λ₂_ext2 = cat( λ₂, zeros(N₁,1), dims=2 ) 
    
    # Reshape λ matrix of 2-dimensional temperature field to vector
    λ₂_ext1 = reshape(λ₂_ext1, N₁*N₂, 1 )
    λ₂_ext2 = reshape(λ₂_ext2, N₁*N₂, 1 )

    # Calculate diffusion
    diff_x1 = 1/(Δx₁^2) * (λ₁_ext1 .* Mx1_lo +  λ₁_ext2 .* Mx1_hi)
    diff_x2 = 1/(Δx₂^2) * (λ₂_ext1 .* Mx2_lo +  λ₂_ext2 .* Mx2_hi)

    # Emitted heat flux
    ϕout1 = bound_left * boundary_heatflux.(θ[1,:], θamb, h, k) + bound_right * boundary_heatflux.(θ[end,:], θamb, h, k) # x₁-direction
    ϕout2 = bound_top * boundary_heatflux.(θ[:,end], θamb, h, k) # x₂-direction
    
    # Induced heat flux
    ϕin = bound_under * b_in * u
    
    # Boundary conditions as vectors
    in_out_x1 = 1/Δx₁ * ϕout1
    in_out_x2 = 1/Δx₂ * (ϕout2 + ϕin)

    # Reshape matrix (2D temperature field) to vector
    θ = reshape(θ, N₁*N₂, 1 )

    # Calculate temporal derivative
    dθ .=  ρinv * Cinv * ( (diff_x1 + diff_x2) * θ + in_out_x1 + in_out_x2 )
end

# Sampling time and final time
const Δt = 10^(-3);
const Tfinal = 10.;
const tspan = (0, Tfinal)

# Solve the ODE
prob = ODEProblem( heat_eq, θinit,tspan )
sol = solve(prob,Euler(),dt=Δt,progress=true,save_everystep=false,save_start=true)



# Save and export results
import DelimitedFiles
basepath        = "results/";

input_filename  = basepath * "input_history_scenario_" * string(scenario) * ".txt";
output_filename = basepath *  "output_history_scenario_" * string(scenario) * ".txt";
temperature2D_filename = basepath *  "temperature2D_scenario_" * string(scenario) * ".txt";

# Save signals with time steps (t_hist)
u_io_data = [t_hist transpose(u_hist)]
y_io_data = [t_hist transpose(y_hist)]

# Save temperature as matrix (2D temperature field)
temperature2D = reshape(sol.u[2], N₁, N₂); 

# Export files
DelimitedFiles.writedlm(input_filename, u_io_data, "\t");
DelimitedFiles.writedlm(output_filename, y_io_data, "\t");
DelimitedFiles.writedlm(temperature2D_filename, temperature2D, "\t");

