using DifferentialEquations,Plots,StaticArrays

λᵢ = @SVector [0.0127,0.0317,0.115,0.311,1.4,3.87]
βᵢ = @SVector [0.00247,0.0013845,0.001222,0.0026455,0.000832,0.000169]
β = sum(βᵢ)
Λ = 0.00001
Ρ(t) = 10^-5*(1.5-t)
ρ = @SVector [0.0001,-0.001,Ρ]



struct Params
    ρ
    λᵢ
    βᵢ
    β
    Λ
end


tspan = (0.0,5)
param = Params(ρ[1],λᵢ,βᵢ,β,Λ)
u₀ = zeros(7)
u₀[1] = 1.0
u₀
ΔT = 0.000001


# function PonctualKinetics!(du,u,p,t)
#     Cᵢ = @views u[2:7]
#     du[1] = ((p.ρ - p.β)/p.Λ)*u[1] + sum(p.λᵢ .* Cᵢ)
#     for i in eachindex(p.λᵢ)
#         du[i+1] = p.βᵢ[i]/p.Λ * u[1] - p.λᵢ[i] * Cᵢ[i]
        
#     end
#end
A


function PonctualKinetics!(du,u,p,t)
    Cᵢ = u[2:7]
    du[1] = ((p.ρ - p.β)/p.Λ)*u[1] + sum(p.λᵢ .* Cᵢ)
    for i in eachindex(p.λᵢ)
        du[i+1] = p.βᵢ[i]/p.Λ * u[1] - p.λᵢ[i] * Cᵢ[i]
        
    end

end

prob = ODEProblem(PonctualKinetics!,u₀,tspan,param)

sol = solve(prob,Tsit5(),dt=ΔT)
x = collect(0:5:0.0001)
plot(sol.(x),idxs = 1)
plot(sol[1])
sol.(x)

begin
	function lotka(du,u,p,t)
	  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
	  du[2] = -p[3]*u[2] + p[4]*u[1]*u[2]
	end
	
	p = [1.5,1.0,3.0,1.0]
	prob1 = ODEProblem(lotka,[1.0,1.0],(0.0,10.0),p)
	sol1 = solve(prob1)
	plot(sol1)
end