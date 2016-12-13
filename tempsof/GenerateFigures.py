#from numpy import *
#from scipy import *
import matplotlib.pyplot as myplot
import scipy as sp
import scipy.io
import numpy as np


myplot.close('all')
FNAME=""
FIGDIR=""

# Load Data eta, Z(:,1), Z(:,2), Z(:,3) ------ eta, p(eta), q(eta), r(eta)
mat = scipy.io.loadmat("pqr_2.mat")

# data was backward in time. change it in order of forward in time
eta = mat['eta'][::-1]
Z = mat['Z'][::-1]
n = mat['n'][0,0]
m = mat['m'][0,0]

# preserve original variables
eta_ = eta[:,0];
p_ = Z[:,0]
q_ = Z[:,1]
r_ = Z[:,2]

# determine translation factor
eta_ = eta_ - 25
xi_ = sp.exp(eta_)

# determine where to cutoff
#CUT2 = size(eta,1);
#CUT2 = floor(CUT2/1);
#CUT2 = 4800;
# CUT  = 135;
CUT2 = 4800
CUT = 200
# cutoff the orbit
eta  = eta_[CUT:CUT2]
xi   = xi_[CUT:CUT2]
# xi   = xi - 2*xi(1);

p    = p_[CUT:CUT2]
q    = q_[CUT:CUT2]
r    = r_[CUT:CUT2]

# tilded variables stress, strain, strain rates, vertical velocity, functions of \xi
tildegamma =       (( (p          )* (r**n)                                  ).real)**(1/(2-n));
tildesigma =       (( (p**(-(1-n)))* (r**n)                                  ).real)**(1/(2-n));
tildev     = 1/n * (( (p**(-(1-n)))* (r**n) ))**(1/(2-n))*np.sign(q)*np.abs(q);
tildeu     =       (( (p          )* (r**2)                                  ).real)**(1/(2-n));

# Bar Capital variables : stress, strain, strain rates, vertical velocity, functions of \xi
VBar     = ( (xi)**(-n/(2-n)   ) ) * (tildev);
GammaBar = ( (xi)**(-2/(2-n)   ) ) * (tildegamma);
UBar     = ( (xi)**(-2/(2-n)   ) ) * (tildeu);
SigmaBar = ( (xi)**( 1-n/(2-n) ) ) * (tildesigma);

# Odd extension
xi       = np.concatenate((-xi[::-1]      , xi      ))
VBar     = np.concatenate((-VBar[::-1]    , VBar    ))
# Even extension
GammaBar = np.concatenate((GammaBar[::-1] , GammaBar))
SigmaBar = np.concatenate((SigmaBar[::-1] , SigmaBar))
UBar     = np.concatenate((UBar[::-1]     , UBar    ))


# open figure windows
fg = myplot.figure()
fgp = fg.add_subplot(111)
fu = myplot.figure()
fup = fu.add_subplot(111)
fv = myplot.figure()
fvp = fv.add_subplot(111)
fs = myplot.figure()
fsp = fs.add_subplot(111)

poly=['k','b','r','g','m']
mrk=['o','x','^','d','s']

tv = [0, 0.1, 0.2, 0.3, 0.4, 0.5]; js=-1; it=0
for t in tv:         # for strain
    g0  = 1.
    tau = sp.log(1 + t/g0);
    fac = (1+t/g0)**m;
    x   = xi / fac;

    gamma = (g0*(1+t/g0))     * (fac**(  2/(2-n)   ))   * GammaBar

    #R = 25;
    R = 3.5
    fgp.plot(x,gamma,'b',linewidth=2);  
    fgp.axis([-R,R,1e-2,4]);
    fgp.set_ylabel('${\gamma}(x,t)$',fontsize=30)
    fgp.set_xlabel('$x$',fontsize=25)
    fgp.tick_params(labelsize=15)
    fgp.set_yscale('log')
    fgp.annotate('t='+str(tv[it]),xy=(0,max(gamma)),xycoords='data',xytext=(js*75,25),fontsize=15, textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    js*=-1; it+=1;
    figname=FIGDIR+'strain_log.eps'
    fg.savefig(figname, format='eps')


tv = [0, 0.1, 0.2, 0.3, 0.5];  js=-1; it=0  
for t in tv:     # for strain rate
    g0  = 1.
    tau = sp.log(1 + t/g0);
    fac = (1+t/g0)**m;
    x   = xi / fac;
    u   = (fac**(  2/(2-n)   ))   * UBar;
    
    #R = 25;
    R = 3.5
    #fup.plot(x,u,'b');
    fup.plot(x,u,linewidth=2, markersize=8, marker=mrk[it],markevery=200, label='t='+str(tv[it]));
    fup.axis([-R,R,.8e-1,15],fontsize=30);
    fup.set_ylabel('$u(x,t)$',fontsize=30)
    fup.set_xlabel('$x$',fontsize=25)
    fup.tick_params(labelsize=15)
    fup.set_yscale('log')
    fup.legend(loc=2)
    #fup.annotate('t='+str(tv[it]),xy=(0,max(u)),xycoords='data',xytext=(js*75,25),fontsize=15, textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    js*=-1; it+=1;    
    figname=FIGDIR+'strain_rate_log.eps'
    fu.savefig(figname, format='eps')
  
tv = [ 0, 0.1, 0.2, 0.3, 0.5]; js=-1; it=0  

for t in tv:       # for v and sigma 
    g0  = 1.
    tau = sp.log(1 + t/g0);
    fac = (1+t/g0)**m;
    x   = xi / fac;
    
    sigma = (1/(g0*(1+t/g0))) * (fac**( -1+n/(2-n) ))   * SigmaBar;
    v     =                (fac**(  n/(2-n)   ))   * VBar; 

    R = 25;
    fvp.plot(x,v,linewidth=2, markersize=8, marker=mrk[it], markevery=200, label='t='+str(tv[it]));  
    #fvp.plot(x,v,'b')
    fvp.axis([-R,R,-4,4]);
    fvp.set_ylabel('$v(x,t)$',fontsize=30)
    fvp.set_xlabel('$x$',fontsize=25)
    fvp.tick_params(labelsize=15)
    fvp.legend(loc=2)
    #fvp.annotate('t='+str(tv[it]),xy=(js*x[2300+js*it*150],js*v[2300+js*it*150]),xycoords='data',xytext=(js*100,-45),fontsize=15,textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    figname=FIGDIR+'velocity.eps'
    fv.savefig(figname, format='eps')   
    
    fsp.plot(x,sigma,'b',linewidth=2);  
    fsp.axis([-R,R,5e-1,200]); 
    fsp.set_ylabel('${\sigma}(x,t)$', fontsize=30)
    fsp.set_xlabel('$x$',fontsize=25)
    fsp.tick_params(labelsize=15)
    fsp.set_yscale('log')
    yl=10; yr=40; ly=(yr-yl)/2.; my=(yr+yl)/2.;
    fsp.annotate('t='+str(tv[it]),xy=(0,min(sigma)),xycoords='data',xytext=(js*85,ly*js+my),fontsize=15,textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    js*=-1; it+=1;
    figname=FIGDIR+'stress_log.eps'
    fs.savefig(figname, format='eps')

myplot.show()

#
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
#
