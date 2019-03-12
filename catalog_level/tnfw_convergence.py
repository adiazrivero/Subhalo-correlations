import numpy as np

"""def rs(msub):
    return 0.1*(msub/10**6)**(1/3)"""

"""def mnfw(m,tau):
    n = np.asarray([(i**2/(i**2 + 1)**2) * ((i**2 - 1) * np.log(i) + i * np.pi - (i**2 + 1)) for i in tau])
    mnfw = m/n
    return mnfw"""

"""def rt_func(m,r3d):
    return (m/10**6)**(1/3) * (r3d/100)**(2/3)"""

"""r_s = rs(masses)
avg_min_rs = np.mean([min(i) for i in r_s])
r3d_arr = np.asarray([np.asarray(i) for i in r_3d])
rt = rt_func(masses,r3d_arr)
max_rt = [max(i) for i in rt]
avg_max_rt = np.mean(max_rt)
t = rt/r_s
m_nfw = mnfw(masses,t)"""

def FNFW(x):
    up = x>1
    down = np.logical_not(up)
    lowx = np.arccosh(1/x[down])/np.sqrt(1 - x[down]**2)
    highx = np.arccos(1/x[up])/np.sqrt(x[up]**2 - 1)
    ans = np.zeros(x.shape)
    ans[down] = lowx
    ans[up] = highx
    return ans

def conv_tnfw(x,y,tau,mnfw,rs,sigmac):
    z = np.sqrt(x**2+y**2)/rs
    Lz = np.log(z/(tau+np.sqrt(tau**2+z**2)))
    Fz = FNFW(z)
    sig = (mnfw/(sigmac*rs**2))*(tau**2/(2*np.pi*(tau**2 + 1)**2))*(((tau**2 + 1)/(z**2 - 1))*(1 - Fz) + 2*Fz - np.pi/np.sqrt(tau**2 + z**2) + ((tau**2 - 1)/(tau*np.sqrt(tau**2 + z**2)))*Lz)
    return sig
