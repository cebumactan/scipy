#from numpy import *
#from scipy import *
import matplotlib.pyplot as myplot
import scipy as sp
import scipy.io
import numpy as np
from scipy.integrate import odeint


myplot.close('all')
FNAME=""
FIGDIR=""

# Load Data eta, Z(:,1), Z(:,2), Z(:,3) ------ eta, p(eta), q(eta), r(eta)
mat = scipy.io.loadmat("pqr.mat")

# data was backward in time. change it in order of forward in time
eta = mat['eta'][::-1]
Z = mat['Z'][::-1]
n = mat['n'][0,0]
A = mat['A'][0,0]
L = mat['L'][0,0]

# preserve original variables
eta_ = eta[:,0]
p_ = Z[:,0]
q_ = Z[:,1]
r_ = Z[:,2]

#figure(1)
#myplot.plot(eta_,p_,eta_,q_,eta_,r_)

# determine translation factor
eta_ = eta_ - 40
xi_ = sp.exp(eta_)

# determine where to cutoff
CUT2 = 4800
CUT = 200
# cutoff the orbit
eta  = eta_[CUT:CUT2]
xi   = xi_[CUT:CUT2]
# xi   = xi - 2*xi(1);

p    = p_[CUT:CUT2]
q    = q_[CUT:CUT2]
r    = r_[CUT:CUT2]

#figure(2)
#myplot.plot(eta,p,eta,q,eta,r)

# constants
m = 0
D = 1+2*A-m-n
a0,a1 = (2+2*A-n)   / D   , (2+2*A)       / D 
b0,b1 = (1+m)       / D   , (1+m+n)       / D 
c0,c1 = (2+2*m)     / D   , (2+2*m+2*n)   / D 
d0,d1 = (-2*A+2*m+n)/ D   , (-2*A+2*m+2*n)/ D 
a = a0 + a1*L
b = b0 + b1*L
c = c0 + c1*L
d = d0 + d1*L

# tilded variables stress, strain, strain rates, vertical velocity, functions of \xi
tildev     = 1/b * ( (p**(-(A-n)))* (r**(n*(1+A))    ))**(1/D)*np.sign(q)*np.abs(q)
tildetheta =       ( (p**(1+n)   )* (r**(n*(1+n))    ))**(1/D)
tildesigma =       ( (p**(-(A-n)))* (r**(n*(1+A))    ))**(1/D)
tildeu     =       ( (p**(1+A)   )* (r**(n*(1+A)+1)  ))**(1/D)

# solve for tildeg
myeta = eta[::-1]
idx = np.size(myeta)-1
mytg_infty = 0

def myodefunc( mytg, myeta, a0, myL, mytu):
    global idx
    dtgdeta = (mytu[idx] - a0*mytg)/myL
    print(idx)
    idx-=1
    return dtgdeta
tildeg = odeint(myodefunc, mytg_infty, myeta, args=(a0, L, tildeu))

# Bar Capital variables : stress, strain, strain rates, vertical velocity, functions of \xi
VBar     = ( (xi)**(-b1   ) ) * (tildev);
ThetaBar = ( (xi)**(-c1   ) ) * (tildetheta);
SigmaBar = ( (xi)**(-d1 ) ) * (tildesigma);
UBar     = ( (xi)**(-b1-1   ) ) * (tildeu);


## solve for Gbar
#myUBar = UBar
#def myodefunc( myG, myxi, a, myL, myU ):
#    print(myxi)
#    dGdxi = (myU[myxi] - a*myG)/(myL*myxi)
#    return dGdxi
#myG_infty = 1
#myxi=xi
#GBar = odeint(myodefunc, myG_infty, myxi, args=(a, L, myUBar))


# Odd extension
xi       = np.concatenate((-xi[::-1]      , xi      ))
VBar     = np.concatenate((-VBar[::-1]    , VBar    ))
# Even extension
ThetaBar = np.concatenate((ThetaBar[::-1] , ThetaBar))
SigmaBar = np.concatenate((SigmaBar[::-1] , SigmaBar))
UBar     = np.concatenate((UBar[::-1]     , UBar    ))


# open figure windows
fth = myplot.figure()
fthp = fth.add_subplot(111)
fu = myplot.figure()
fup = fu.add_subplot(111)
fv = myplot.figure()
fvp = fv.add_subplot(111)
fs = myplot.figure()
fsp = fs.add_subplot(111)

fthp.plot(xi,ThetaBar)
fup.plot(xi,UBar)
fvp.plot(xi,VBar)
fsp.plot(xi,SigmaBar)


poly=['k','b','r','g','m']
mrk=['o','x','^','d','s','*']

tv = [ 0, 1, 3, 7, 12]; js=-1; it=0  
for t in tv:         # for temperature
    
    fac   = (t+1)**L;
    theta = (t+1)**c * ThetaBar
    x     = xi / fac;

    #R = 25;
    R = 50
    fthp.plot(x,theta,'b', linewidth=2);#, markersize=8, marker=mrk[it], markevery=200, label='t='+str(tv[it]));  
    fthp.axis([-R,R,1.0/80,10]);
    fthp.set_ylabel(r"${\theta}(x,t)$",fontsize=30)
    fthp.set_xlabel('$x$',fontsize=25)
    fthp.tick_params(labelsize=15)
    fthp.set_yscale('log')
    print(tv[it])
    fthp.annotate('t='+str(tv[it]),xy=(0,max(theta)),xycoords='data',xytext=(js*120,25),fontsize=15, textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    js*=-1; it+=1;
    figname=FIGDIR+'temperature_log.eps'
    fth.savefig(figname, format='eps')

tv = [ 0, 1, 3, 7, 12]; js=-1; it=0 
#tv = [ 0, 2, 4, 6, 8, 10]; js=-1; it=0   
for t in tv:     # for strain rate
    fac   = (t+1)**L;
    u = (t+1)**(b+L) * UBar
    x     = xi / fac;
    

    R = 10
    #fup.plot(x,u,'b');
    fup.plot(x,u,linewidth=2, markersize=8, marker=mrk[it],markevery=200, label='t='+str(tv[it]));
    fup.axis([-R,R,1e-2,3.5],fontsize=30);
    fup.set_ylabel('$u(x,t)$',fontsize=30)
    fup.set_xlabel('$x$',fontsize=25)
    fup.tick_params(labelsize=15)
    fup.set_yscale('log')
    fup.legend(loc=2)
    fup.annotate('t='+str(tv[it]),xy=(0,max(u)),xycoords='data',xytext=(js*120,25),fontsize=15, textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    js*=-1; it+=1;    
    figname=FIGDIR+'strain_rate_log.eps'
    fu.savefig(figname, format='eps')

tv = [ 0, 1, 3, 7, 12]; js=-1; it=0     
#tv = [ 0, 2, 4, 6, 8, 10]; js=-1; it=0  
for t in tv:       # for v and sigma 
    fac   = (t+1)**L;
    v     = (t+1)**b * VBar
    x     = xi / fac;

    R = 30;
    fvp.plot(x,v,linewidth=2, markersize=8, marker=mrk[it], markevery=200, label='t='+str(tv[it]));  
    #fvp.plot(x,v,'b')
    fvp.axis([-R,R,-2.2,2.2]);
    fvp.set_ylabel('$v(x,t)$',fontsize=30)
    fvp.set_xlabel('$x$',fontsize=25)
    fvp.tick_params(labelsize=15)
    fvp.legend(loc=2)
    it+=1
    #fvp.annotate('t='+str(tv[it]),xy=(js*x[2300+js*it*150],js*v[2300+js*it*150]),xycoords='data',xytext=(js*100,-45),fontsize=15,textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    figname=FIGDIR+'velocity.eps'
    fv.savefig(figname, format='eps')   
    
tv = [ 0, 1, 3, 7, 12]; js=-1; it=0 
#tv = [ 0, 2, 4, 6, 8, 10]; js=-1; it=0  
for t in tv:       # for v and sigma 
    fac   = (t+1)**L;
    sigma = (t+1)**d * SigmaBar
    x     = xi / fac;

    R = 60;

    fsp.plot(x,sigma,'b',linewidth=2);  
    fsp.axis([-R,R,1e-1,20]); 
    fsp.set_ylabel("${\sigma}(x,t)$", fontsize=30)
    fsp.set_xlabel('$x$',fontsize=25)
    fsp.tick_params(labelsize=15)
    fsp.set_yscale('log')
    yl=10; yr=40; ly=(yr-yl)/2.; my=(yr+yl)/2.;
    fsp.annotate('t='+str(tv[it]),xy=(0,min(sigma)),xycoords='data',xytext=(js*120,ly*js+my),fontsize=15,textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    js*=-1; it+=1;
    figname=FIGDIR+'stress_log.eps'
    fs.savefig(figname, format='eps')
#myplot.show()


#figure(1);
#xlabel('$x$','interpreter','latex'); 
#ylabel('${\gamma}(x,t)$','interpreter', 'latex');    
#set(get(gca,'YLabel'),'Rotation',0);
#ylabh = get(gca,'yLabel');
#set(ylabh, 'Units', 'Normalized')
#set(ylabh, 'Position',get(ylabh,'Position').*[2,1,1]);
#hold off;
#figure(2);
#xlabel('$x$','interpreter','latex'); ylabel('${u}(x,t)$','interpreter', 'latex');      
#set(get(gca,'YLabel'),'Rotation',0);
#ylabh = get(gca,'yLabel');
#set(ylabh, 'Units', 'Normalized')
#set(ylabh, 'Position',get(ylabh,'Position').*[2,1,1])
#hold off;
#figure(3);
#xlabel('$x$','interpreter','latex'); ylabel('${v}(x,t)$','interpreter', 'latex');      
#set(get(gca,'YLabel'),'Rotation',0);
#ylabh = get(gca,'yLabel');
#set(ylabh, 'Units', 'Normalized')
#set(ylabh, 'Position',get(ylabh,'Position').*[2,1,1])
#hold off;
#figure(4);
#xlabel('$x$','interpreter','latex'); ylabel('${\sigma}(x,t)$','interpreter', 'latex');  
#set(get(gca,'YLabel'),'Rotation',0);
#ylabh = get(gca,'yLabel');
#set(ylabh, 'Units', 'Normalized')
#set(ylabh, 'Position',get(ylabh,'Position').*[1.7,1,1])
#hold off;

