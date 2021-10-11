##################
# SIR Model
##################

# Modules
using DifferentialEquations
using Plots
using LinearAlgebra

#--- parameters

# N = population size
N = 85000000
# beta = Infection Rate
beta = 0.211
# gamma = multiplicative inverse of Avg. Infectious Period
gamma = 1/12

p = [beta gamma]

#--- ODE problem

function SIR(du,u,p,t)
    beta, gamma = p
    du[1] = - beta*u[2]*u[1]
    du[2] = beta*u[2]*u[1] - gamma*u[2]
    du[3] = gamma*u[2]
end

#--- Initial states of S,I,R

S = N * (999999/1000000)
I = N * (1/1000000) # 1 in a million is infected
R = 0

u0 = [S,I,R]
u0_sum = sum(u0)

u0 = u0/u0_sum #normalizing u0
tspan = (0. , 750. ) #time period

#--- R0 calculation

s = u0[1]

# positive part of derivative of I respect to I
F = beta*s

#negative part of derivative of I respect to I
V = gamma

R0 = F/V
R0_round = round(R0,digits=1)

#--- Solving the problem and making its graph

Problem = ODEProblem(SIR, u0, tspan, p)
solution = solve(Problem)

plot(solution, title = "SIR Model", xticks = 0:tspan[2]/10:tspan[2], label = ["S" "I" "R"])
xlabel!("Days")
annotate!(tspan[2]/2,1.03,"R₀ = $R0_round")


savefig("C:\\Users\\merts\\Google Drive\\3.semester\\Grimm\\fig1")
