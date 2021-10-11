##################
# SEI³Q³RD Model
##################

# Modules
using DifferentialEquations
using Plots
using LinearAlgebra

#--- sigma calculation over time

function sigma_t!(sigma,u_6,u_9,B_by_grp,N,N_len)

    for k in 1:N_len
        if (u_6[k] + u_9[k])*N[k] > B_by_grp[k]
            sigma[k] = (sigma[k]*B_by_grp[k] + (u_6[k] + u_9[k])*N[k] - B_by_grp[k])/((u_6[k] + u_9[k])*N[k])
        end
    end
end

#--- Parameters

# Number of individuals in groups
N = #[59500000 , 25500000]
# Number of individuals in whole environment
N_sum = sum(N)
# Number of groups
N_len = length(N)

# beta_asym =  Infection rate matrix of asymptomatic group
beta_asym = #[.258 .250
            #.250 .242]
# beta_sym = Infection rate matrix of symptomatic group
beta_sym = #[.238 .230
            #.230 .222]
# beta_sev = Infection rate matrix of severe group
beta_sev = #[.078 .070
            #.070 .062]

# gamma_asym = Multiplicative Inverse of Avg. Infectious period of asymptomatic Infectious individuals
gamma_asym = #1/10
# gamma_sym = group specific Inverse of Avg. Infectious period of symptomatic Infectious individuals
gamma_sym = #[1/10 1/15]
# gamma_sev_r = group specific Inverse of Avg. Infectious period of severely symptomatic Infectious individuals with recovery
gamma_sev_r =  #[1/25 1/25]
# gamma_sev_d = group specific Inverse of Avg. Infectious period of severely symptomatic Infectious individuals with death
gamma_sev_d =  #[1/20 1/20]

# epsilon = Multiplicative inverse of Avg. Latent Period
epsilon =  #[1/5.2 1/5.2]

# eta = Fraction of asymptomatic infectious individuals by groups
eta =  #[.200 .002]

# nu = Fraction of severely symptomatic infectious individuals to symptomatic individuals by groups.
nu = #[.050 .516]

#sigma = Lethality rates by group
sigma = #[.009 .090]
# we do not want to lose initial sigma to use on R0 calculation
sigma_init = copy(sigma)

# Available beds
B = #39100
# Bed by groups
B_by_grp = N*B/N_sum

# psi = Percentage of tracing app users by groups
psi = #[2.5/10 2.5/10]

p =(N,N_len,beta_asym,beta_sym,beta_sev,gamma_asym,gamma_sym,gamma_sev_r,gamma_sev_d,epsilon,eta,nu,sigma,B_by_grp,psi)

#--- ODE problem

function SEI3Q3RD!(du,u,p,t)

    N,N_len,beta_asym,beta_sym,beta_sev,gamma_asym,gamma_sym,gamma_sev_r,gamma_sev_d,epsilon,eta,nu,sigma,B_by_grp,psi = p

    u_1 =  u[:,:,1] #S
    u_2 =  u[:,:,2] #E_nt
    u_3 =  u[:,:,3] #E_tr
    u_4 =  u[:,:,4] #I_asym
    u_5 =  u[:,:,5] #I_sym
    u_6 =  u[:,:,6] #I_sev
    u_7 =  u[:,:,7] #Q_asym
    u_8 =  u[:,:,8] #Q_sym
    u_9 =  u[:,:,9] #Q_sev
    u_10 = u[:,:,10] #R
    u_11 = u[:,:,11] #D

    sigma_t!(sigma,u_6,u_9,B_by_grp,N,N_len)

    du[:,:,1] = -(u_4*beta_asym .+ u_5*beta_sym .+ u_6*beta_sev).*u_1
    du[:,:,2] =  (u_4*beta_asym .+ u_5*beta_sym .+ u_6*beta_sev .- (psi.*((psi.*u_5)*beta_sym) .+ psi.*((psi.*u_6)*beta_sev))).*u_1 .- epsilon.*u_2
    du[:,:,3] = (psi.*((psi.*u_5)*beta_sym) .+ psi.*((psi.*u_6)*beta_sev)).*u_1 .-epsilon.*u_3
    @. du[:,:,4] = (eta*epsilon)*u_2 - gamma_asym*u_4
    @. du[:,:,5] = (1 - eta)*(1 - nu)*epsilon*u_2 - gamma_sym*u_5
    @. du[:,:,6] = (1 - eta)*nu*epsilon*u_2 - ((1 - sigma)*gamma_sev_r + sigma*gamma_sev_d)*u_6
    @. du[:,:,7] = eta*epsilon*u_3 - gamma_asym*u_7
    @. du[:,:,8] = (1 - eta)*nu*epsilon*u_3 - gamma_sym*u_8
    @. du[:,:,9] = (1 - eta)*nu*epsilon*u_3 - ((1 - sigma)*gamma_sev_r + sigma*gamma_sev_d)*u_9
    @. du[:,:,10] = gamma_asym*(u_4 + u_7) + gamma_sym*(u_5 + u_8) + (1 - sigma)*gamma_sev_r*(u_6 + u_9)
    @. du[:,:,11] = sigma*gamma_sev_d*u_6 + sigma*gamma_sev_d*u_9

end

#--- initial states of S, E, I_asym, I_sym, I_sev, R, D by groups

S_0 = #N * (999997.4/1000000)
E_nt_0 = #N * (1.6/1000000) * (1 - (psi[1]*N[1]/N_sum + psi[2]*N[2]/N_sum))
E_tr_0 = #N * (1.6/1000000) * (psi[1]*N[1]/N_sum + psi[2]*N[2]/N_sum)
I_asym_0 =  #N * (1/1000000) *(20/100)
I_sym_0 = #N * (1/1000000) * (65/100)
I_sev_0 = #N * (1/1000000) * (15/100)
Q_asym_0 = #[0,0]
Q_sym_0 = #[0,0]
Q_sev_0 = #[0,0]
R_0 = #[0,0]
D_0 = #[0,0]

u0 = zeros(1,N_len,11)
u0[:,:,1] = S_0
u0[:,:,2] = E_nt_0
u0[:,:,3] = E_tr_0
u0[:,:,4] = I_asym_0
u0[:,:,5] = I_sym_0
u0[:,:,6] = I_sev_0
u0[:,:,7] = Q_asym_0
u0[:,:,8] = Q_sym_0
u0[:,:,9] = Q_sev_0
u0[:,:,10] = R_0
u0[:,:,11] = D_0

u0_sum = sum(u0)
u0 = u0/u0_sum #normalizing u0
tspan = #(0. , 750.) #time period

#--- R0 calculation

s = u0[:,:,1]

# positive part of derivative of not traced Exposed groups and Infected groups respect to themselves
F = [0 0 beta_asym[1,1]*s[1] beta_asym[2,1]*s[1] beta_sym[1,1]*s[1] beta_sym[2,1]*s[1] beta_sev[1,1]*s[1] beta_sev[2,1]*s[1] #E_nt_1
    0 0 beta_asym[1,2]*s[2] beta_asym[2,2]*s[2] beta_sym[1,2]*s[2] beta_sym[2,2]*s[2] beta_sev[1,2]*s[2] beta_sev[2,2]*s[2] #E_nt_2
    (eta.*epsilon)[1] 0 0 0 0 0 0 0 #I_asym_1
    0 (eta.*epsilon)[2] 0 0 0 0 0 0 #I_asym_2
    ((1 .+ (eta.*nu)).*epsilon)[1] 0 0 0 0 0 0 0 #I_sym_1
    0 ((1 .+ (eta.*nu)).*epsilon)[2] 0 0 0 0 0 0 #I_sym_2
    (nu.*epsilon)[1] 0 0 0 0 0 (sigma_init.*gamma_sev_r)[1] 0 #I_sev_1
    0 (nu.*epsilon)[2] 0 0 0 0 0 (sigma_init.*gamma_sev_r)[2] #I_sev_2
    ]

# negative part of derivative of not traced Exposed groups and Infected groups respect to themselves
V = [epsilon[1] 0 0 0 0 0 0 0 #E_nt_1
    0 epsilon[2] 0 0 0 0 0 0 #E_nt_2
    0 0 gamma_asym 0 0 0 0 0 #I_asym_1
    0 0 0 gamma_asym 0 0 0 0 #I_asym_2
    ((eta.+nu).*epsilon)[1] 0 0 0 gamma_sym[1] 0 0 0 #I_sym_1
    0 ((eta.+nu).*epsilon)[2] 0 0 0 gamma_sym[2] 0 0 #I_sym_2
    ((eta.*nu).*epsilon)[1] 0 0 0 0 0 (gamma_sev_r .+ (sigma_init.*gamma_sev_d))[1] 0 #I_sev_1
    0 ((eta.*nu).*epsilon)[2] 0 0 0 0 0 (gamma_sev_r .+ (sigma_init.*gamma_sev_d))[2] #I_sev_2
     ]

# Finding spectral radius of F*V_inv
V_inv = inv(V)
eigval = eigvals(F*V_inv)
eigval_abs = zeros(length(eigval))
for i in 1:length(eigval)
    eigval_abs[i] = abs(eigval[i])
end
R0_index = argmax(eigval_abs)

R0 = eigval_abs[R0_index]
R0_round = round(R0,digits=2)

#--- Solving the problem and making its graph

Problem = ODEProblem(SEI3Q3RD!, u0, tspan, p)
solution = solve(Problem)

plot(solution, layout=2, ylim = [0,1], yticks= 0:0.05:1 , title = ["SEI³Q³RD Low-Risk Group" "SEI³Q³RD High-Risk Group"], label = ["S" "S" "E_nt" "E_nt" "E_tr" "E_tr" "I_asym" "I_asym" "I_sym" "I_sym" "I_sev" "I_sev" "Q_asym" "Q_asym" "Q_sym" "Q_sym" "Q_sev" "Q_sev" "R" "R" "D" "D"])
xlabel!("Days")
annotate!(tspan[2]/2,1.03,"R₀=$R0_round")
