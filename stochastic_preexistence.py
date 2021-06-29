import numpy as np
import scipy.stats as ss
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dosingPeriod",type=int, required=True)
args = parser.parse_args()
dose_interval = args.dosingPeriod

class StochasticSimulation(object):
    """
    """
    def __init__(self, Lambda=2e6, d=0.1, beta=5e-10, k=1000, c=10, a=0.5, u=3e-5, s=0.3):
        """
        """
        self.Lambda = Lambda
        self.d = d
        self.beta = beta
        self.k = k
        self.c = c
        self.a = a
        self.u = u
        self.s = s
        self.R = (Lambda*beta*k)/(a*d*c)

    def pre_treatment(self, t_final, number_iterations):
        """
        """
        x0 = self.Lambda/self.d
        ys0 = 100
        yr0 = 0
        vs0 = 0
        vr0 = 0
        x = [x0,ys0,yr0,vs0,vr0]

        number_iterations = 1
        final_count = np.zeros((number_iterations,len(x)))
        t = 0.0
        interval = 0.001
        plot_times = np.arange(0,t_final+1,1)
        xt = np.zeros((number_iterations,len(plot_times)))
        yst = np.zeros((number_iterations,len(plot_times)))
        yrt = np.zeros((number_iterations,len(plot_times)))
        vst = np.zeros((number_iterations,len(plot_times)))
        vrt = np.zeros((number_iterations,len(plot_times)))
        recorded_time = 0.0
        for num in range(0,number_iterations,1):
            while t<t_final:
                if recorded_time%1000.0==0:
                    index = int(recorded_time/1000.0)
                    xt[num][index] = x[0]
                    yst[num][index] = x[1]
                    yrt[num][index] = x[2]
                    vst[num][index] = x[3]
                    vrt[num][index] = x[4]

                birth_x = self.Lambda;death_x=self.d*x[0]+self.beta*x[0]*x[3]+self.beta*x[0]*x[4]*(1-self.s)
                birth_ys = self.beta*(1-self.u)*x[0]*x[3]; death_ys=self.a*x[1]
                birth_yr = self.beta*self.u*x[0]*x[3]+self.beta*x[0]*x[4]*(1-self.s); death_yr=self.a*x[2]
                birth_vs = self.k*x[1]; death_vs=self.c*x[3]
                birth_vr = self.k*x[2]; death_vr=self.c*x[4]

                x[0] += ss.poisson.rvs(birth_x*interval)-ss.poisson.rvs(death_x*interval)
                x[1] += ss.poisson.rvs(birth_ys*interval)-ss.poisson.rvs(death_ys*interval)
                x[2] += ss.poisson.rvs(birth_yr*interval)-ss.poisson.rvs(death_yr*interval)
                x[3] += ss.poisson.rvs(birth_vs*interval)-ss.poisson.rvs(death_vs*interval)
                x[4] += ss.poisson.rvs(birth_vr*interval)-ss.poisson.rvs(death_vr*interval)

                for i in range(len(x)):
                    if x[i]<0:
                        x[i] = 0
                t += interval
                recorded_time += 1
            final_count[num] = x

        return x
    def epsilon(self, avg, A, T, t):
        """
        """
        return avg-A*np.cos(2*np.pi*t/T)

    def during_treatment(self, x, eavgs, As, eavgr, Ar, T_dose, t_final, number_iterations):
        """
        """

        interval = 0.001
        establishment_probability = []


        T = T_dose
        fraction_alive_mutants = 0.0

        x0 = x[0]; ys0 = x[1]; yr0 = 1; vs0 = x[3]; vr0 = 0


        final_count_treat = np.zeros((number_iterations,len(x)))
        t = 0.0

        for num in range(0,number_iterations,1):
            x_treat = [x0,ys0,yr0,vs0,vr0]
            while t<t_final:

                birth_x = self.Lambda;death_x=self.d*x_treat[0]+(1-self.epsilon(eavgs,As,T,t))*self.beta*x_treat[0]*x_treat[3]+self.beta*x_treat[0]*x_treat[4]*(1-self.s)*(1-self.epsilon(eavgr,Ar,T,t))
                birth_ys = self.beta*(1-self.epsilon(eavgs,As,T,t))*x_treat[0]*x_treat[3]; death_ys = self.a*x_treat[1]
                birth_yr = self.beta*(1-self.epsilon(eavgr,Ar,T,t))*x_treat[0]*x_treat[4]*(1-self.s); death_yr = self.a*x_treat[2]
                birth_vs = self.k*x_treat[1]; death_vs = self.c*x_treat[3]
                birth_vr = self.k*x_treat[2]; death_vr = self.c*x_treat[4]

                x_treat[0]+=ss.poisson.rvs(birth_x*interval)-ss.poisson.rvs(death_x*interval)
                x_treat[1]+=ss.poisson.rvs(birth_ys*interval)-ss.poisson.rvs(death_ys*interval)
                x_treat[2]+=ss.poisson.rvs(birth_yr*interval)-ss.poisson.rvs(death_yr*interval)
                x_treat[3]+=ss.poisson.rvs(birth_vs*interval)-ss.poisson.rvs(death_vs*interval)
                x_treat[4]+=ss.poisson.rvs(birth_vr*interval)-ss.poisson.rvs(death_vr*interval)

                for i in range(len(x_treat)):
                    if x_treat[i]<0:
                        x_treat[i]=0.0
                if t > 5 and x_treat[4]==0.0:
                    break
                t+=interval
            final_count_treat[num]=x_treat
            t=0.0

        for i in range(number_iterations):
            if final_count_treat[i][4]!=0.0:
                fraction_alive_mutants+=1

        establishment_probability.append(fraction_alive_mutants/number_iterations)

        return establishment_probability


sim = StochasticSimulation()
start_time = time.time()
result1 = sim.pre_treatment(200,1)
result2 = sim.during_treatment(result1,0.900,0.090,0.000,0.000,dose_interval,100,20000)
print("--- %s seconds ---" % (time.time() - start_time))
print(result2)
print("T=%i"%dose_interval)
