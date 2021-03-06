import numpy as np

"""def mb(m,tau,p):
    #n = np.asarray([(i**2 * (np.pi * (p-i)**2 + 4*i**2 * np.log(p/i)))/(4 * (p**4 - i**4)) for i in tau])
    n = np.asarray((tau**2 * (np.pi * (p-tau)**2 + 4*tau**2 * np.log(p/tau)))/(4 * (p**4 - tau**4)))
    m_b = m/n
    return m_b"""

def term1(x,p,tau):
    return (2*p*np.sqrt(1/(tau**2+x**2)))/(p**4 - tau**4)

def term2(x,p,tau):
    down = x < p
    up = np.logical_not(down)
    lowx = np.sqrt(1/(x[up]**2-p**2))/(p*(p**2+tau**2))
    highx = 0
    ans = np.zeros(x.shape)
    ans[up] = lowx
    ans[down] = highx
    return ans

def term3(x,p,tau):
    return np.sqrt(1/(x**2+p**2))/(p**3-p*tau**2)

def term5(x,p,tau):
    return 2*np.arctanh(p/np.sqrt(x**2+p**2)) / ((p**3-p*tau**2)*np.sqrt(x**2+p**2))

def term6(x,p,tau):
    return 4*tau*np.arctanh(tau/np.sqrt(x**2+tau**2)) / ((p**4-tau**4)*np.sqrt(x**2+tau**2))

def term4(x,p,tau):
    down = x < p
    up = np.logical_not(down)
    lowx = 2*np.arctan(p/np.sqrt(x[up]**2-p**2)) / ((p**3+p*tau**2)*np.sqrt(x[up]**2 - p**2))
    highx = np.real(2*np.arctan(p/(1j*np.sqrt(p**2-x[down]**2))) / ((p**3+p*tau**2)*1j*np.sqrt(p**2-x[down]**2)))
    ans = np.zeros(x.shape)
    ans[up] = lowx
    ans[down] = highx
    return ans

def conv_tburk(x,y,rs,mb,pval,tauval,sigmac):
    p = pval
    t = tauval
    z = np.sqrt(x**2+y**2)/rs
    term_1 = term1(z,p,t)
    term_2 = term2(z,p,t)
    term_3 = term3(z,p,t)
    term_4 = term4(z,p,t)
    term_5 = term5(z,p,t)
    term_6 = term6(z,p,t)
    sig = (mb/(8*np.pi*sigmac*rs**2))*t**2 * (np.pi*(term_1 - term_2 - term_3) + term_4 - term_5 + term_6)
    return sig

