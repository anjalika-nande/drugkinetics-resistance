#Defining the global variables
import numpy as np
import time
import scipy.stats as ss


# Pharma model
def pharma_drug(Cmax,th,T,t):
    return Cmax*2**(-np.mod(t,T)/th)

def pharma_efficacy(Cmax,th,T,t,IC50,M,rho):
    eff = 1/(1+(pharma_drug(Cmax,th,T,t)/(rho*IC50))**M)
    return eff

start_time = time.time()

t_final = 300 ; t_eq = 50
Lambda=2e6; d=0.1; beta=5e-10; k=1000; c=10; a=0.5; u=3e-5*1e2; s=0.3; R=(Lambda*beta*k)/(a*d*c); b = beta*k/c; gamma = 1
Cmax = 17; th = 120; T = 120
IC50 = 10; M = 10; rhow = 1; rhor = np.inf

x0=Lambda/d; yw0=0; yr0=0;

fraction_alive_mutants=0.0

number_points = 100
number_iterations = 200

final_count=np.zeros((number_points,number_iterations,3))

t=0.0
interval=0.001

for point in range(number_points):
    t_point = T*point/number_points
    for num in range(number_iterations):
        x=[x0,yw0,yr0]

        while t < t_eq+t_point:

            if x[1] < 0:
                x[1] = 0

            birth_x=Lambda;death_x=d*x[0]+(pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*x[0]*x[1]
            birth_yw=gamma*(1-u/s) + (pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*(1-u)*x[0]*x[1]; death_yw=a*x[1]

            x[0]+=ss.poisson.rvs(birth_x*interval)-ss.poisson.rvs(death_x*interval)
            x[1]+=ss.poisson.rvs(birth_yw*interval)-ss.poisson.rvs(death_yw*interval)

            t+=interval

        while t_eq+t_point <= t < t_eq + t_point + (1/gamma):

            if x[1] < 0:
                x[1] = 0

            if x[2] < 0:
                x[2] = 0

            birth_x=Lambda;death_x=d*x[0]+(pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*x[0]*x[1]+ b*x[0]*x[2]*(1-s)
            birth_yw=gamma*(1-u/s) + (pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*(1-u)*x[0]*x[1]; death_yw=a*x[1]
            birth_yr= gamma*u/s + (pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*(u)*x[0]*x[1]+b*x[0]*x[2]*(1-s); death_yr=a*x[2]

            x[0]+=ss.poisson.rvs(birth_x*interval)-ss.poisson.rvs(death_x*interval)
            x[1]+=ss.poisson.rvs(birth_yw*interval)-ss.poisson.rvs(death_yw*interval)
            x[2]+=ss.poisson.rvs(birth_yr*interval)-ss.poisson.rvs(death_yr*interval)

            if x[1] < 0:
                x[1] = 0

            if x[2] < 0:
                x[2] = 0

            t+=interval

        while t_eq + t_point + (1/gamma) <= t < t_eq + t_point + (1/gamma) + 60:

            if x[1] < 0:
                x[1] = 0

            if x[2] < 0:
                x[2] = 0

            if x[2]==0:
                break

            birth_x=Lambda;death_x=d*x[0]+(pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*x[0]*x[1]+ b*x[0]*x[2]*(1-s)
            birth_yw=gamma*(1-u/s) + (pharma_efficacy(Cmax,th,T,t,IC50,M,rhow))*b*(1-u)*x[0]*x[1]; death_yw=a*x[1]
            birth_yr= b*x[0]*x[2]*(1-s); death_yr=a*x[2]

            x[0]+=ss.poisson.rvs(birth_x*interval)-ss.poisson.rvs(death_x*interval)
            x[1]+=ss.poisson.rvs(birth_yw*interval)-ss.poisson.rvs(death_yw*interval)
            x[2]+=ss.poisson.rvs(birth_yr*interval)-ss.poisson.rvs(death_yr*interval)

            if x[1] < 0:
                x[1] = 0

            if x[2] < 0:
                x[2] = 0

            t+=interval

        final_count[point][num]=x
        t=0.0

for i in range(number_points):
    for j in range(number_iterations):
        if final_count[i][j][2]!=0.0:
            fraction_alive_mutants+=1

print(fraction_alive_mutants/(number_iterations*number_points))
print("T=%i"%T)

print("--- %s seconds ---" % (time.time() - start_time))
