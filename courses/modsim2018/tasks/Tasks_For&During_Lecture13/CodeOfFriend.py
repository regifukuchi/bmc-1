{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Task 12 - Motor Control\n",
    "### Introduction to modeling and simulation of human movement\n",
    "https://github.com/BMClab/bmc/blob/master/courses/ModSim2018.md"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Desiree Miraldo  \n",
    "Renato Watanabe\n",
    "\n",
    "\n",
    "* Based on task for Lecture 11 (+ muscle activation dynamics + pennation angle):\n",
    "\n",
    "Change the derivative of the contractile element length function. The new function must compute the derivative according to the article from Thelen (2003) (Eqs. (1), (2), (6) and (7)):\n",
    "\n",
    "     Thelen D; Adjustment of muscle mechanics model parameters to simulate dynamic contractions in older adults (2003)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "#import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "import math\n",
    "\n",
    "%matplotlib notebook"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Muscle properties"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "Lslack = .223\n",
    "Umax = .04\n",
    "Lce_o = .093 #optmal l\n",
    "width = .63#*Lce_o\n",
    "Fmax = 7400\n",
    "a = 0\n",
    "u = 0.5\n",
    "#b = .25*10#*Lce_o "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Initial conditions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "Lnorm_ce = .087/Lce_o #norm\n",
    "t0 = 0\n",
    "tf = 5\n",
    "h = 1e-3"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "t = np.arange(t0,tf,h)\n",
    "F = np.empty(t.shape)\n",
    "Fkpe = np.empty(t.shape)\n",
    "FiberLen = np.empty(t.shape)\n",
    "TendonLen = np.empty(t.shape)\n",
    "a_dynamics = np.empty(t.shape)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Simulation - Series"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "for i in range (len(t)):\n",
    "    #ramp\n",
    "    if t[i]<=1:\n",
    "        Lm = 0.31\n",
    "    elif t[i]>1 and t[i]<2:\n",
    "        Lm = .31 + .1*(t[i]-1)\n",
    "        #print(Lm)\n",
    "    \n",
    "    #shortening at 4cm/s\n",
    "    Lsee = Lm - Lce\n",
    "    \n",
    "    if Lsee<Lslack: \n",
    "        F[i] = 0\n",
    "    else: \n",
    "        F[i] = Fmax*((Lsee-Lslack)/(Umax*Lslack))**2\n",
    "        \n",
    "        \n",
    "    #isometric force at Lce from CE force length relationship\n",
    "    F0 = max([0, Fmax*(1-((Lce-Lce_o)/width)**2)])\n",
    "    \n",
    "    #calculate CE velocity from Hill's equation\n",
    "    if  F[i]>F0: print('Error: cannot do eccentric contractions')\n",
    "    \n",
    "    Lcedot = -b*(F0-F[i])/(F[i]+a) #vel is negative for shortening\n",
    "    \n",
    "    # --- Euler integration step\n",
    "    Lce += h*Lcedot\n",
    "\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "def TendonForce (Lnorm_see,Lslack, Lce_o):\n",
    "    '''\n",
    "    Compute tendon force\n",
    "\n",
    "    Inputs:\n",
    "        Lnorm_see = normalized tendon length\n",
    "        Lslack = slack length of the tendon (non-normalized)\n",
    "        Lce_o = optimal length of the fiber\n",
    "    \n",
    "    Output:\n",
    "        Fnorm_tendon = normalized tendon force\n",
    "        \n",
    "    '''\n",
    "    Umax = .04\n",
    "    \n",
    "    if Lnorm_see<Lslack/Lce_o: \n",
    "        Fnorm_tendon = 0\n",
    "    else: \n",
    "        Fnorm_tendon = ((Lnorm_see-Lslack/Lce_o)/(Umax*Lslack/Lce_o))**2\n",
    "        \n",
    "    return Fnorm_tendon"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "def ParallelElementForce (Lnorm_ce):\n",
    "    '''\n",
    "    Compute parallel element force\n",
    "    \n",
    "    Inputs:\n",
    "        Lnorm_ce = normalized contractile element length\n",
    "    \n",
    "    Output:\n",
    "        Fnorm_kpe = normalized parallel element force\n",
    "\n",
    "    '''\n",
    "    Umax = 1\n",
    "    \n",
    "    if Lnorm_ce< 1: \n",
    "        Fnorm_kpe = 0\n",
    "    else: \n",
    "        Fnorm_kpe = ((Lnorm_ce-1)/(Umax*1))**2 \n",
    "        \n",
    "    return Fnorm_kpe"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [],
   "source": [
    "def ForceLengthCurve (Lnorm_ce,width):\n",
    "    F0 = max([0, (1-((Lnorm_ce-1)/width)**2)])\n",
    "    return F0"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {
    "collapsed": true
   },
   "source": [
    "def ContractileElementDot(F0, Fnorm_CE, a, b):\n",
    "    \n",
    "    '''\n",
    "    Compute Contractile Element Derivative\n",
    "\n",
    "    Inputs:\n",
    "        F0 = Force-Length Curve\n",
    "        Fce = Contractile element force\n",
    "    \n",
    "    Output:\n",
    "        Lnorm_cedot = normalized contractile element length derivative\n",
    "\n",
    "    '''\n",
    "    \n",
    "    if  Fnorm_CE>F0: print('Error: cannot do eccentric contractions')\n",
    "    \n",
    "    Lnorm_cedot = -b*(F0-Fnorm_CE)/(Fnorm_CE + a) #vel is negative for shortening\n",
    "    \n",
    "    return Lnorm_cedot"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [],
   "source": [
    "def ContractileElementDot(F0, Fnorm_CE, a):\n",
    "    \n",
    "    '''\n",
    "    Compute Contractile Element Derivative\n",
    "\n",
    "    Inputs:\n",
    "        F0 = Force-Length Curve\n",
    "        Fce = Contractile element force\n",
    "    \n",
    "    Output:\n",
    "        Lnorm_cedot = normalized contractile element length derivative\n",
    "\n",
    "    '''\n",
    "    \n",
    "    FMlen = 1.4 # young adults\n",
    "    Vmax = 10  # young adults\n",
    "    Af = 0.25  #force-velocity shape factor\n",
    "    \n",
    "    Fnorm_CE = min(FMlen*a*F0 - 0.001, Fnorm_CE)\n",
    "    \n",
    "    if  Fnorm_CE > a*F0:\n",
    "        \n",
    "        b = ((2 + 2/Af)*(a*F0*FMlen - Fnorm_CE))/(FMlen-1)\n",
    "        \n",
    "    elif Fnorm_CE <= a*F0:\n",
    "        \n",
    "        b = a*F0 + Fnorm_CE/Af\n",
    "    \n",
    "    Lnorm_cedot = (.25 + .75*a)*Vmax*((Fnorm_CE - a*F0)/b)\n",
    "    \n",
    "    return Lnorm_cedot"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [],
   "source": [
    "def ContractileElementForce(Fnorm_tendon,Fnorm_kpe, alpha):\n",
    "    '''\n",
    "    Compute Contractile Element force\n",
    "\n",
    "    Inputs:\n",
    "        Fnorm_tendon = normalized tendon force\n",
    "        Fnorm_kpe = normalized parallel element force\n",
    "    \n",
    "    Output:\n",
    "        Fnorm_CE = normalized contractile element force\n",
    "    '''\n",
    "    Fnorm_CE = Fnorm_tendon/np.cos(alpha) - Fnorm_kpe\n",
    "    return Fnorm_CE"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [],
   "source": [
    "def tendonLength(Lm,Lce_o,Lnorm_ce, alpha):\n",
    "    '''\n",
    "    Compute tendon length\n",
    "    \n",
    "    Inputs:\n",
    "        Lm = \n",
    "        Lce_o = optimal length of the fiber\n",
    "        Lnorm_ce = normalized contractile element length\n",
    "    \n",
    "    Output:\n",
    "        Lnorm_see = normalized tendon length   \n",
    "    '''\n",
    "    Lnorm_see = Lm/Lce_o - Lnorm_ce*np.cos(alpha)\n",
    "    \n",
    "    return Lnorm_see"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [],
   "source": [
    "def activation(a,u,dt):\n",
    "    '''\n",
    "    Compute activation\n",
    "    \n",
    "    Inputs:\n",
    "        u = idealized muscle excitation signal, 0 <= u <= 1\n",
    "        a = muscular activation\n",
    "        dt = time step\n",
    "    \n",
    "    Output:\n",
    "        a = muscular activation  \n",
    "    '''\n",
    "    \n",
    "    tau_deact = 50e-3 #young adults\n",
    "    tau_act = 15e-3\n",
    "    \n",
    "    if u>a:\n",
    "        tau_a = tau_act*(0.5+1.5*a)\n",
    "    elif u <=a:\n",
    "        tau_a = tau_deact/(0.5+1.5*a)\n",
    "    \n",
    "    #-------\n",
    "    dadt = (u-a)/tau_a # euler\n",
    "    \n",
    "    a += dadt*dt\n",
    "    #-------\n",
    "    return a"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Simulation - Parallel"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [],
   "source": [
    "#Normalizing\n",
    "alpha = 25*np.pi/180\n",
    "for i in range (len(t)):\n",
    "    #ramp\n",
    "    if t[i]<=1:\n",
    "        Lm = 0.31\n",
    "    elif t[i]>1 and t[i]<2:\n",
    "        Lm = .31 - .04*(t[i]-1)\n",
    "        #print(Lm)\n",
    "    \n",
    "    #shortening at 4cm/s\n",
    "    u = 0.7 + 0.2*np.sin(np.pi*t[i])\n",
    "    \n",
    "    Lnorm_see = tendonLength(Lm,Lce_o,Lnorm_ce, alpha)\n",
    "\n",
    "    Fnorm_tendon = TendonForce (Lnorm_see,Lslack, Lce_o) \n",
    "    \n",
    "    Fnorm_kpe = ParallelElementForce (Lnorm_ce)     \n",
    "        \n",
    "    #isometric force at Lce from CE force length relationship\n",
    "    F0 = ForceLengthCurve (Lnorm_ce,width)\n",
    "    \n",
    "    Fnorm_CE = ContractileElementForce(Fnorm_tendon,Fnorm_kpe, alpha) #Fnorm_CE = ~Fm\n",
    "    \n",
    "    #computing activation\n",
    "    a = activation(a,u,h)\n",
    "    \n",
    "    #calculate CE velocity from Hill's equation    \n",
    "    Lnorm_cedot = ContractileElementDot(F0, Fnorm_CE,a)\n",
    "    \n",
    "    # --- Euler integration step\n",
    "    Lnorm_ce += h*Lnorm_cedot\n",
    "\n",
    "    \n",
    "    F[i] = Fnorm_tendon*Fmax\n",
    "    Fkpe[i] = Fnorm_kpe*Fmax\n",
    "    FiberLen[i] = Lnorm_ce*Lce_o\n",
    "    TendonLen[i] = Lnorm_see*Lce_o\n",
    "    a_dynamics[i] = a\n",
    "    "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Plots "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "/opt/miniconda3/lib/python3.6/site-packages/matplotlib/axes/_axes.py:545: UserWarning: No labelled objects found. Use label='...' kwarg on individual plots.\n",
      "  warnings.warn(\"No labelled objects found. \"\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<matplotlib.figure.Figure at 0x7fc7a72a3cf8>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig, ax = plt.subplots(1, 1, figsize=(6,6), sharex=True)\n",
    "\n",
    "ax.plot(t,a_dynamics,c='magenta')\n",
    "plt.grid()\n",
    "plt.xlabel('time (s)')\n",
    "plt.ylabel('Activation dynamics')\n",
    "\n",
    "\n",
    "ax.legend()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "/opt/miniconda3/lib/python3.6/site-packages/matplotlib/axes/_axes.py:545: UserWarning: No labelled objects found. Use label='...' kwarg on individual plots.\n",
      "  warnings.warn(\"No labelled objects found. \"\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<matplotlib.figure.Figure at 0x7fc79c391be0>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig, ax = plt.subplots(1, 1, figsize=(6,6), sharex=True)\n",
    "\n",
    "ax.plot(t,F,c='red')\n",
    "plt.grid()\n",
    "plt.xlabel('time (s)')\n",
    "plt.ylabel('Force (N)')\n",
    "\n",
    "\n",
    "ax.legend()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<matplotlib.text.Text at 0x7fc79c1ef550>"
      ]
     },
     "execution_count": 15,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<matplotlib.figure.Figure at 0x7fc7a7258630>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<matplotlib.figure.Figure at 0x7fc79c1c3160>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig, ax = plt.subplots(1, 1, figsize=(6,6), sharex=True)\n",
    "\n",
    "ax.plot(t,FiberLen, label = 'fiber')\n",
    "ax.plot(t,TendonLen, label = 'tendon')\n",
    "plt.grid()\n",
    "plt.xlabel('time (s)')\n",
    "plt.ylabel('Length (m)')\n",
    "ax.legend(loc='best')\n",
    "\n",
    "\n",
    "fig, ax = plt.subplots(1, 3, figsize=(12,4), sharex=True, sharey=True)\n",
    "ax[0].plot(t,FiberLen, label = 'fiber')\n",
    "ax[1].plot(t,TendonLen, label = 'tendon')\n",
    "ax[2].plot(t,FiberLen + TendonLen, label = 'muscle (tendon + fiber)')\n",
    "\n",
    "ax[1].set_xlabel('time (s)')\n",
    "ax[0].set_ylabel('Length (m)')\n",
    "#plt.legend(loc='best')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.3"
  },
  "latex_envs": {
   "LaTeX_envs_menu_present": true,
   "autoclose": false,
   "autocomplete": true,
   "bibliofile": "biblio.bib",
   "cite_by": "apalike",
   "current_citInitial": 1,
   "eqLabelWithNumbers": true,
   "eqNumInitial": 1,
   "hotkeys": {
    "equation": "Ctrl-E",
    "itemize": "Ctrl-I"
   },
   "labels_anchors": false,
   "latex_user_defs": false,
   "report_style_numbering": false,
   "user_envs_cfg": false
  },
  "nbTranslate": {
   "displayLangs": [
    "*"
   ],
   "hotkey": "alt-t",
   "langInMainMenu": true,
   "sourceLang": "en",
   "targetLang": "fr",
   "useGoogleTranslate": true
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}