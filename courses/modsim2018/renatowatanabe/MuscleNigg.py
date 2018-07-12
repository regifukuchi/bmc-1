class MuscleNigg(object):
    
    def __init__(self, Lslack, Lce_o, Fmax, dt):
        self.Lslack = Lslack
        self.Lce_o = Lce_o
        self.Fmax = Fmax
        
        self.dt = dt
        
        
        self.UmaxTendon = 0.04
        self.UmaxParallel = 100
        self.widthFLcurve = .63
        self.AFvcurve = 0.25
        self.F_m_lenFVcurve = 1.4 #young
        #F_m_len = 1.8 #old
        self.v_maxFVcurve = 10 #young
        #v_max = 8 #old
        self.tau_dact = 0.050 #young #ms
        #tal_dact = 0.060 #old #ms
        self.tau_act = 0.015 #ms
        
        self.LseeNorm = 0
        self.FTendonNorm = 0
        self.FkpeNorm = 0
        self.LceNorm = 0
        self.FLNorm = 0
        self.FCENorm = 0
        self.act = 0
        self.LceNormdot = 0
        
        
        
    def computeTendonForce(self):
        '''
        Compute Tendon Force

        Inputs:

        LseeNorm - Normalized Tendon  length
        Lslack - slack length of the tendon (non-normalized)
        Lce_o - Optimal length of the fiber

        Output:
        FtendonNorm - Normalized tendon force
        '''

        
        if self.LseeNorm<self.Lslack/self.Lce_o: 
            self.FTendonNorm = 0
        else: 
            self.FTendonNorm = ((self.LseeNorm-self.Lslack/self.Lce_o)/(self.UmaxTendon*self.Lslack/self.Lce_o))**2

        
    
    def computeParallelElementForce(self):
        
        if self.LceNorm < 1: 
            self.FkpeNorm = 0
        else: 
            self.FkpeNorm = ((self.LceNorm-1)/(self.UmaxParallel))**2 

           
    def computeForceLengthCurve(self):
        
        self.FLNorm = max([0, (1-((self.LceNorm-1)/self.widthFLcurve)**2)])        
    
    def computeContractileElementDerivative(self):
        
        trueFCE = min((self.F_m_lenFVcurve*self.act*self.FLNorm-1e-4),self.FCENorm)
        if  trueFCE <= self.act*self.FLNorm:
            b = self.act*self.FLNorm + trueFCE  / self.AFvcurve            
        else:
            b = (2 + 2/self.AFvcurve)*(self.act*self.FLNorm*self.F_m_lenFVcurve - trueFCE)/(self.F_m_lenFVcurve - 1)
            
        self.LceNormdot = (0.25 + 0.75*self.act)*self.v_maxFVcurve*(trueFCE   - self.act*self.FLNorm)/b
   
    def computeContractileElementForce(self):
        self.FCENorm = self.FTendonNorm - self.FkpeNorm
      
    def computeTendonLength(self, Lm):
        self.LseeNorm = Lm/self.Lce_o - self.LceNorm
        
    
    def computeActivation(self, n):      
 
        if n > self.act:
            tau_a = self.tau_act*(0.5 + 1.5*self.act)
        else:
            tau_a = self.tau_dact/(0.5 + 1.5*self.act)
    
    
        dActdt = (n-self.act)/tau_a
        
        self.act = self.act + self.dt*dActdt

    def updateMuscle(self, Lm, n):
        self.computeTendonLength(Lm)    
        self.computeTendonForce()        
        self.computeParallelElementForce()      
        self.computeForceLengthCurve()        
        self.computeContractileElementForce()    
        self.computeActivation(n)    
        self.computeContractileElementDerivative()

        self.LceNorm = self.LceNorm + self.dt*self.LceNormdot
