##################
# SEIR Model
##################

# Modules
using DifferentialEquations
using Plots
using LinearAlgebra

#--- parameters

# N = population size
N = #85000000
# beta = Infection Rate
beta = #0.211
# gamma =  inverse of Avg. Infectious Period
gamma = #1/12
# epsilon = inverse of Avg. Latent Period
epsilon = #1/5.2

p = [beta gamma epsilon]

#--- ODE problem

function SEIR(du,u,p,t)
    beta, gamma, epsilon = p
    du[1] = -beta*u[3]*u[1]
    du[2] =  beta*u[3]*u[1] - epsilon*u[2]
    du[3] = epsilon*u[2] - gamma*u[3]
    du[4] = gamma*u[3]
end

#--- initial states of S, E, I, R

S = #N * (999997.4/1000000)
E = #N * (1.6/1000000) # 1.6 in a million is exposed
I = #N * (1/1000000) # 1 in a million is infected
R = 0

u0 = [S,E,I,R]
u0_sum = sum(u0)
u0 = u0/u0_sum # normalizing u0
tspan = #(0. , 750. ) # time period

#--- R0 calculation

#parameters for R0

s = u0[1]

# positive part of the derivative of E and I respect to themselves
F = [0. beta*s  #E
    epsilon 0]  #I

# negative part of the derivative of E and I respect to themselves
V = [epsilon 0. #E
    0 gamma]    #I

#Finding spectral radius of F*V_inv
V_inv = inv(V)
eigval = eigvals(F*V_inv)
eigval_abs = zeros(length(eigval))
for i in 1:length(eigval)
    eigval_abs[i] = abs(eigval[i])
end
R0_index = argmax(eigval_abs)

R0 = eigval_abs[R0_index]
R0_round = round(R0,digits=1)


#--- Solving the problem and making its graph

Problem = ODEProblem(SEIR, u0, tspan, p)
solution = solve(Problem)

plot(solution, title = "SEIR Model",xticks = 0:tspan[2]/10:tspan[2], yticks = 0:0.25:1, label = ["S" "E" "I" "R"])
xlabel!("Days")
annotate!(tspan[2]/2,1.03,"R₀=$R0_round")
