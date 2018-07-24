import numpy as np
class Muscle:
    
    def __init__(self, Lslack, Lce_o, Fmax, alpha, dt):
        
        self.Lslack = Lslack # tendon slack length
        self.Lce_o  = Lce_o # optimal muscle fiber length
        self.Fmax   = Fmax #maximal isometric DF force
        self.alpha  = alpha # DF muscle fiber pennation angle
        
        
        self.UmaxTendon = 0.04
        self.Umax = 1
        self.width = 0.63 # Max relative length change of CE
        self.FMlen = 1.4 # young adults
        self.Vmax = 10  # young adults
        self.Af = 0.25  #force-velocity shape factor
        
        self.tau_deact = 50e-3 #young adults
        self.tau_act = 15e-3
        self.dt = dt
        
        self.a = 0.01 #inital conditional for ativation
        self.Lnorm_see = 0
        self.Lnorm_ce = 0
        self.Fnorm_tendon = 0
        self.Fnorm_kpe = 0
        self.F0 = 0
        self.Lnorm_cedot = 0
        
    def tendonLength(self, Lm):
        '''
        Compute tendon length

        Inputs:
            Lm = 
            Lce_o = optimal length of the fiber
            Lnorm_ce = normalized contractile element length

        Output:
            Lnorm_see = normalized tendon length   
        '''
        self.Lnorm_see = Lm/self.Lce_o - self.Lnorm_ce*np.cos(self.alpha)

        return self.Lnorm_see
    
    def TendonForce(self):
        '''
        Compute tendon force

        Inputs:
            Lnorm_see = normalized tendon length
            Lslack = slack length of the tendon (non-normalized)
            Lce_o = optimal length of the fiber

        Output:
            Fnorm_tendon = normalized tendon force

        '''
        

        if self.Lnorm_see<self.Lslack/self.Lce_o: 
            self.Fnorm_tendon = 0
        else: 
            self.Fnorm_tendon = ((self.Lnorm_see-self.Lslack/self.Lce_o)/
                            (self.UmaxTendon*self.Lslack/self.Lce_o))**2

        return self.Fnorm_tendon
    
    def ParallelElementForce(self):
        '''
        Compute parallel element force

        Inputs:
            Lnorm_ce = normalized contractile element length

        Output:
            Fnorm_kpe = normalized parallel element force

        '''
        

        if self.Lnorm_ce< 1: 
            self.Fnorm_kpe = 0
        else: 
            self.Fnorm_kpe = ((self.Lnorm_ce-1)/(self.Umax*1))**2 

        return self.Fnorm_kpe
    
    def ForceLengthCurve(self):
        self.F0 = max([0, (1-((self.Lnorm_ce-1)/self.width)**2)])
        
        return self.F0
    
    def ContractileElementDot(self):
    
        '''
        Compute Contractile Element Derivative

        Inputs:
            F0 = Force-Length Curve
            Fce = Contractile element force

        Output:
            Lnorm_cedot = normalized contractile element length derivative

        '''

        

        self.Fnorm_CE = max(0.001,min(self.FMlen*self.a*self.F0 - 0.001, self.Fnorm_CE))

        if self.Fnorm_CE > self.a*self.F0:

            b = ((2 + 2/self.Af)*(self.a*self.F0*self.FMlen - self.Fnorm_CE))/(self.FMlen-1)

        elif self.Fnorm_CE <= self.a*self.F0:

            b = self.a*self.F0 + self.Fnorm_CE/self.Af

        self.Lnorm_cedot = (.25 + .75*self.a)*self.Vmax*((self.Fnorm_CE - self.a*self.F0)/b)

        return self.Lnorm_cedot
    
    def ContractileElementForce(self):
        '''
        Compute Contractile Element force

        Inputs:
            Fnorm_tendon = normalized tendon force
            Fnorm_kpe = normalized parallel element force

        Output:
            Fnorm_CE = normalized contractile element force
        '''
        self.Fnorm_CE = self.Fnorm_tendon/np.cos(self.alpha) - self.Fnorm_kpe
        
        return self.Fnorm_CE
    
    def activation(self, u):
        '''
        Compute activation

        Inputs:
            u = idealized muscle excitation signal, 0 <= u <= 1
            a = muscular activation
            dt = time step

        Output:
            a = muscular activation  
        '''

        
        if u>self.a:
            tau_a = self.tau_act*(0.5+1.5*self.a)
        elif u <= a:
            tau_a = self.tau_deact/(0.5+1.5*self.a)

        #-------
        dadt = (u-self.a)/tau_a # euler

        self.a = self.a + dadt*self.dt
        #-------
        return self.a
    
    def updateMuscle(self, Lm, u):
        self.tendonLength(Lm)

        self.TendonForce() 

        self.ParallelElementForce()     

        #isometric force at Lce from CE force length relationship
        self.ForceLengthCurve()

        self.ContractileElementForce() #Fnorm_CE = ~Fm

        #computing activation
        self.activation(u)

        #calculate CE velocity from Hill's equation    
        self.ContractileElementDot()

        self.Lnorm_ce = self.Lnorm_ce + self.dt*self.Lnorm_cedot