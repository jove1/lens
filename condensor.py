#!/usr/bin/python3
from pylab import *

sdac = float(0x4000)

def saturationf(x,b):
    return where(b != 0, arctan(x*b)/b, x)

def model(p1,p2,p3,p4, b1,b2,b3,b4, bb1,bb2,bb3, b21,b12, dac, cla, full=False):
    "CL1 + CL3 split by CLA + CM + OBJ prefield"
    B1 = b1*saturationf(dac[0]/sdac, bb1)
    B2 = b2*saturationf(dac[1]/sdac, bb2)
    B3 = b3*saturationf(dac[2]/sdac, bb3)
    B4 = b4*dac[3]/sdac
    
    # important cross terms between CL1 and CL3!
    B1, B2 = B1 + b21*B2, B2 + b12*B1

    t = []
    y, dy = 0, 1; t.append(y)
    y = y + p1*dy; t.append(y)
    dy = dy - y*B1*B1  # cl1
    y = y + (p2-p1)*dy; t.append(y)
    dy = dy - y*B2*B2  # cl3
    y = y + (p3-p2)*dy; t.append(y)
    dy = dy - y*B3*B3  # cm
    y = y + (p4-p3)*dy; t.append(y)
    dy = dy - y*B4*B4  # obj. prefield
    y = y + (1-p4)*dy; t.append(y)

    if full:
        return array([0,p1,p2,p3,p4,1]),  array(broadcast_arrays(*t))*cla/2./abs(t[2])
    return t[-1]*cla/2./abs(t[2])

def fit(all_data, params, fit_params, plot_cov=True):

    def fun(vals, log=True):
        params.update(zip(fit_params, vals))
        ret = concatenate([ d[:,3] - sscreen*model(dac=(d[:,0], d[:,1], d[:,2], obj), cla=cla, **params) for d,i,l,sscreen,cla,obj in all_data])
        if log:
            print(sqrt(sum(ret**2)/len(ret)), tuple(vals))
        return ret

    from scipy.optimize import leastsq
    print("-"*80)
    vals, cov, _, mesg, code = leastsq(fun, [params[n] for n in fit_params], full_output=True)
    print("-"*80)
    print(code, mesg.replace("\n  "," "))
    print()

    if code in (1,2,3,4):
        for n,v in zip(fit_params, vals):
            params[n] = v
        res = fun(vals,log=False)
        cov *= sum(res**2)/(len(res)-len(fit_params))
        err = sqrt(diag(cov))
        for n,v,e in zip(fit_params, vals, err):
            params[n] = v
            print(n, "=", v, "+-", e)
        
        if plot_cov:
            m = abs(cov).max() 
            matshow(cov, cmap='seismic', vmin=-m, vmax=m)
            colorbar()
          
    else:
        raise RuntimeError

def plt(all_data, params):
    x = linspace(0, 0xffff, 500)
    for ii in 0, 1, 2:
        figure(ii)
        for d,i,l,sscreen,cla,obj in all_data:
            if i != ii:
                continue
            dac = d[0,:3].tolist() + [obj, ]
            dac[i] = x
            l, = plot(x, sscreen*model(dac=dac, cla=cla, **params), "-", label=l)
            plot(d[:,i], d[:,3], "o", c=l.get_color())

        xlim(0, 0xffff)
        #axhline(-512,c='k')
        #axhline(512,c='k')
        ylim(-550, 550)
        #ylim(-5500, 5500)
        import matplotlib.ticker as ticker
        xax = gca().get_xaxis()
        xax.set_major_locator(ticker.MultipleLocator(0x2000))
        xax.set_major_formatter(ticker.FormatStrFormatter("0x%X"))
        xlabel("DAC")
        ylabel("Spot radius")
        legend()
        grid()

    if 0:
        x /= sdac
        figure()
        plot(x, x,"k-",lw=0.5)
        plot(x, arctan(x*params["bb1"])/params["bb1"], label="CL1")
        plot(x, arctan(x*params["bb2"])/params["bb2"], label="CL3")
        plot(x, arctan(x*params["bb3"])/params["bb3"], label="CM")
        grid()
        legend() 


def interactive(params):
    from matplotlib.widgets import Slider
    from matplotlib.collections import LineCollection

    a0 = 40

    fig = figure(figsize=(8,8))
    cl1_spot={}
    cm_alpha={1: 0xd510, 2: 0x9510, 3: 0x7c50}

    scl1 = Slider(fig.add_axes([0.1,0.26,0.8,0.05]), 'CL1', 0, 0xffff, 0x7630)#, "0x%04x")
    scl3 = Slider(fig.add_axes([0.1,0.19,0.8,0.05]), 'CL3', 0, 0xffff, 0x9400)#, "0x%04x")
    scm  = Slider(fig.add_axes([0.1,0.12,0.8,0.05]), 'CM',  0, 0xffff, cm_alpha[2])#, "0x%04x")
    solp  = Slider(fig.add_axes([0.1,0.05,0.8,0.05]), 'OL pre',  0, 0xffff, std_focus)#, "0x%04x")

    #scm.ax.set_xticks(cm_alpha.values())
    #scm.ax.set_xticklabels(cm_alpha.keys())

    def update(val=None):
        x,y = model(dac=(scl1.val, scl3.val, scm.val, solp.val), full=True, cla=a0, **params)
        l.set_segments([list(zip(x,s*y)) for s in linspace(-1,1,7)])
        fig.canvas.draw_idle()

    scl1.on_changed(update)
    scl3.on_changed(update)
    scm.on_changed(update)
    solp.on_changed(update)
    """
    from jeol.temcon.chat import TEMCON, asyncore
    
    class Lens(TEMCON):
        def found_packet(self, p):
            if p.cmd == 'N128':
                scl1.set_val(p.data.CL1)
                scl3.set_val(p.data.CL3)
                scm.set_val(p.data.CM)
                solp.set_val(p.data.OL_F)
                
    L = Lens('127.0.0.1', 2010)

    def update_timer():
        asyncore.poll(0)
    timer = fig.canvas.new_timer(interval=50)
    timer.add_callback(update_timer)
    timer.start()
    """

    ax = fig.add_axes([0.1,0.4,0.8,0.5])

    l = LineCollection([], linewidths=0.5, colors='k')
    ax.add_collection(l)
    update()

    color = rcParams['axes.prop_cycle'][:3].by_key()['color']
    plot(0,0,"o",c=color[2])
    axvline(params['p1'])
    axvline(params['p2'])
    vlines([params['p2'],params['p2']],[-a0,a0/2.],[-a0/2.,a0], lw=5, colors=color[1])
    axvline(params['p3'])
    axvline(params['p4'])
    axvline(1, c=color[2])
    #axhline(0,c="k",lw=.7)
    show()

sscreen = {}
sscreen[200] = 200*12036e-6 # nominal
sscreen[200] = 1077.7/5/(25400/300.) # measured nominal grid
sscreen[200] = 1077.7/5/(0.89/0.18*252.3/15) # measured measured grid
sscreen[10000] = 10e3*12036e-6 # nominal
sscreen[10000] = 1040.4/20*2160/1000 # calibrated 2160l/mm grid

cla = {1:200, 2:100, 3:40, 4:10}
# measured relative apperture sizes (nominal 200, 100, 40 and 10um):
#a = array([(892+871)/2., (437+426)/2., (177+170)/2., (54+52)/2.])
#a/a[0]*200 == array([200.        ,  97.90130459,  39.36471923,  12.02495746])

std_focus = 0x841e

all_data = [
    # LOWMAG
    ("lens_cl1_x200_cla3.txt",            0, "CL1",            sscreen[200], cla[3], 0),
    ("lens_cl1_cm0x4000_x200_cla3.txt",   0, "CL1 CM=0x4000",  sscreen[200], cla[3], 0),
    ("lens_cl1_cm0x8000_x200_cla3.txt",   0, "CL1 CM=0x8000",  sscreen[200], cla[3], 0),
    ("lens_cl1_cl30x4000_x200_cla3.txt",  0, "CL1 CL3=0x4000", sscreen[200], cla[3], 0),
    ("lens_cl1_cl30x8000_x200_cla3.txt",  0, "CL1 CL3=0x8000", sscreen[200], cla[3], 0),

    ("lens_cl3_x200_cla3.txt",            1, "CL3",            sscreen[200], cla[3], 0),
    ("lens_cl3_cl10x4a00_x200_cla3.txt",  1, "CL3 CL1=0x4a00", sscreen[200], cla[3], 0),
    ("lens_cl3_cl10x8000_x200_cla3.txt",  1, "CL3 CL1=0x8000", sscreen[200], cla[3], 0),
    ("lens_cl3_cm0x8000_x200_cla3.txt",   1, "CL3 CM=0x8000",  sscreen[200], cla[3], 0),

    ("lens_cm_x200_cla3.txt",             2, "CM",             sscreen[200], cla[3], 0),
    ("lens_cm_cl10x5000_x200_cla3.txt",   2, "CM CL1=0x5000",  sscreen[200], cla[3], 0),
    ("lens_cm_cl10x8000_x200_cla3.txt",   2, "CM CL1=0x8000",  sscreen[200], cla[3], 0),
    ("lens_cm_cl30x8000_x200_cla3.txt",   2, "CM CL3=0x8000",  sscreen[200], cla[3], 0),

    # MAG std. focus
#    ("lens_cl3_x10k_spot1_alpha1_cla1_x.txt", 1, "CL3 spot1 alpha1", sscreen[10000], cla[1], std_focus),
    ("lens_cl3_x10k_spot1_alpha2_cla2.txt", 1, "CL3 spot1 alpha2", sscreen[10000], cla[2], std_focus),
    ("lens_cl3_x10k_spot1_alpha3_cla2.txt", 1, "CL3 spot1 alpha3", sscreen[10000], cla[2], std_focus),

    ("lens_cl1_x10k_cla2.txt",              0, "CL1 10k",      sscreen[10000], cla[2], std_focus),
    ("lens_cl3_x10k_cla2.txt",              1, "CL3 10k",      sscreen[10000], cla[2], std_focus),
#    ("lens_cm_x10k_cla2_x.txt",             2, "CM 10k",       sscreen[10000], cla[2], std_focus),

    ("lens_cm_x10k_cla2.txt",               2, "CM 10k",       sscreen[10000], cla[2], std_focus),
    ("lens_cl3_x10k_spot1_alpha1_cla1.txt", 1, "CL3 spot1 alpha1", sscreen[10000], cla[1], std_focus),


]

if __name__ == "__main__":
    import sys
    params = dict(
            p1=0.4, b1=2, bb1=-0.05,
            p2=0.6, b2=1.5, bb2=-0.007,
            p3=0.9, b3=2, bb3=-0.007,
            p4=0.95, b4=sqrt(10),
            b21=0.07,
            b12=0.07
        )
    
    params = {'b21': 0.049701559347947726, 'b4': 6.953730569906942, 'b1': 2.371673206397684, 'b2': 1.480912999526893, 'b3': 1.9299922411622203, 'p2': 0.6070576187160656, 'p3': 0.875011836415697, 'p1': 0.33511630827997235, 'p4': 0.9949764172436991, 'b12': 0.05368547179114984, 'bb3': 0.13174385003940706, 'bb2': 0.15890324293197633, 'bb1': 0.7578559645183915}



    if "fit" in sys.argv[1:]:
        all_data = [ (loadtxt(t[0]),)+t[1:] for t in all_data]

        #fit(all_data[:13], params, "p1,p2,p3,b1,b2,b3,bb1,bb2,bb3,b12,b21".split(","), plot_cov=0)
        #fit(all_data[:-2], params, "p1,p2,p3,p4,b1,b2,b3,b4,bb1,bb2,bb3,b12,b21".split(","), plot_cov=0)
        print(params)

        plt(all_data, params)
        interactive(params)
        #show()
    
    else:
        interactive(params)
