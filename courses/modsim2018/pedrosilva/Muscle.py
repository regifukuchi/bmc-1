import numpy as np
import math

class Muscle:
    
   
    def __init__(self,Lslack,Lce_o,Fmax,alpha,tau_dact,tau_act,dt,Af,F_m_len,v_max):
        
        self.Lslack = Lslack # tendon slack length
        self.Lce_o  = Lce_o # optimal muscle fiber length
        self.Fmax   = Fmax #maximal isometric DF force
        self.alpha  = alpha # DF muscle fiber pennation angle
        
        #Parameters from Nigg & Herzog (2006).
        self.UmaxTendon = 0.04 # SEE strain at Fmax
        self.Umax = 1
        self.width = 0.63 #max relative length change of CE
        
        self.tau_dact = tau_dact
        self.tau_act  = tau_act
        self.dt = dt
        
        self.Af = Af
        self.F_m_len = F_m_len
        self.v_max = v_max
        
        self.a = 0.001 #initial condition for actiation
        self.Lnomr_see = 0
        self.F_tendon_norm = 0
        self.F_kpe_norm = 0
        self.F0 = 0
        self.F_CE_norm = 0
        
        self.Lnorm_slack = self.Lslack/self.Lce_o #normalized slack length
        
    def TendonLength(self,Lm):
        """
        ComputeTendonLength
        
        Inputs:
        
            Lm - Normalized muscle length
        
            Lce_o - Optimal length of the fiber
        
            Lnorm_ce - Nomalized Contract Element Lenght
        
            alpha - Angle between tendon and fiber
        
        Output:
            Lnomr_see - Normalized slack length
        
        """
    
        self.Lnomr_see = Lm/self.Lce_o - self.Lnorm_ce*np.cos(self.alpha)
    
        return self.Lnomr_see
    
    def TendonForce(self):
        """
        ComputeTendonForce
    
        Inputs:
    
            L_SEE_norm - Normalized tendon length
            
            Lnorm_slack - Normalized slack length of the tendon
            
            Lce_o = Optimal length of the fiber
            
        Output:
            F_tendon_norm - Normalized tendon force
            
        """
    
        if self.Lnomr_see < self.Lnorm_slack:
            self.F_tendon_norm = 0
        else:
            self.F_tendon_norm =  ((self.Lnomr_see - self.Lnorm_slack)/(self.UmaxTendon*self.Lnorm_slack))**2
        
        return self.F_tendon_norm
    
    def ParalelElementForce(self):
        """
        Compute Paralel Elemento Force
        
        Inputs:
        
            Lnorm_ce - Normalized Contract Element length
                                 
        Output:
            F_kpe_norm - Normalized paralel element force
    
        """
      
        if self.Lnorm_ce < 1:
            self.F_kpe_norm = 0
        else:
            self.F_kpe_norm = ((self.Lnorm_ce - 1)/(self.Umax*1))**2
        
        return self.F_kpe_norm
    
    def ForceLengthCurve (self):
        """
        Compute Force Length Curve
        
        Inputs:
        
            Lnorm_ce - Normalized Contract Element length
            
        Output:
            F0 - Normalized Force Length Element
        
        """
    
        self.F0 = max(0,(1-((self.Lnorm_ce-1)/self.width)**2))
        
        return self.F0
    
    def ContractileElementForce(self):
        """"
        Incremental Contract Element's length
    
        Inputs:
    
            F_tendon_norm - Normalized tendon force
        
            F_CE_norm - Normalized Contract Element Force
        
            alpha - Angle between tendon and fibers (rad)
        
        Output:
            F_kpe_norm - Normalized paralel element force
      
        """
        
        self.F_CE_norm = self.F_tendon_norm/np.cos(self.alpha) - self.F_kpe_norm
        
        return self.F_CE_norm
    
    def Activation(self,u):
        """
        Compute Activation parameter
        
        Inputs:
        
            a - Muscle activation parameter
            
            u - Muscle excitation signal
            
            dt - time step
            
        Output:
            a - Muscle activation parameter
        
        """
      
        if u > self.a:
            tau_a = self.tau_act*(0.5 + 1.5*self.a)
        else:
            tau_a = self.tau_dact/(0.5 + 1.5*self.a)
        
        dadt = (u-self.a)/tau_a # euler
    
        self.a = self.a + dadt*self.dt
    
        return self.a
    
    def ContractileElementDot(self):
        """"
        Incremental Contract Element's length
        
        Inputs:
        
            F_Lnorm - Normalized Force Length Element
            
            F_kpe_norm - Normalized paralel element force
            
        Output:
            L_CE_dot - Incremental Contract Element's length
             
        """

    
        FCE = min(self.a*self.F0*self.F_m_len-0.001,self.F_CE_norm)
    
        if FCE <= self.a*self.F0:
            b = self.a*self.F0 + FCE/self.Af
            
        else:
            b = (2 + 2/self.Af)*(self.a*self.F0*self.F_m_len - FCE)/(self.F_m_len - 1)
            
        self.Lnorm_cedot  = (0.25 + 0.75*self.a)*self.v_max*((FCE - self.a*self.F0)/b)

        return self.Lnorm_cedot
    
    def updateMuscle(self, Lm, u):
    
        self.TendonLength(Lm)
    
        self.TendonForce()
        
        self.ParalelElementForce()
        
        self.ForceLengthCurve()
    
        self.ContractileElementForce()
    
        self.Activation(u)
                
        self.ContractileElementDot()
    
        self.Lnorm_ce = self.Lnorm_ce + self.dt*self.Lnorm_cedot
        
        return self.F_tendon_norm, self.a, self.Lnorm_ce