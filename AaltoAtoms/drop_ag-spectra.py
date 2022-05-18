import win32com.client
stm=win32com.client.Dispatch("pstmafm.stmafmrem")
stm.stmbeep()
import pstmafm
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--init_dist', type=float, required=False) # Angstroms
parser.add_argument('--max_dist', type=float, required=False)
args = parser.parse_args()


zpiezoconst = float(stm.getparam('Zpiezoconst')) #Angstroms/V
stm.setparam('FBLogIset', 500) # pA # beware that this number depends on gain 
stm.setparam('BiasVolt.[mV]', 1000) # 
stm.setparam('Vertmangain', 6) #
stm.setp('VERTMAN.SPECLENGTH.SEC', 3)
tipheight_b = stm.signals1data(2, 0.1, 1)
print('\t tip height before VM: %1.5lf angstroms' %tipheight_b)

#vertical manipulation
stm.vert_setextzv(1)

#first line is # of points and offset in DSP clock cycles
# DSP_Clock is typically 50kHz
DSP_clock_freq = float(stm.getparam('DSP_Clock')) # (50000) Hz
print(f'DSP clock frequency {DSP_clock_freq}')
timedelay = 2.5/2501.0 # seconds between vertical ramp points 
tdelay_cycles = int(DSP_clock_freq*timedelay) # DSP cycles 

print("tdelay_cycles: %lf" %tdelay_cycles)
tipcrashdist = 5.1 # Angstroms 
tcd_mv = tipcrashdist/zpiezoconst*1000 # mV


def save_vert_drop_file(tcd_mv: float, tdelay_cycles: float, ramp_points: int=2501, fname: str = r'C:\Users\phys-asp-lab\Dropbox\personal-folders\Abe\drop Ag script\Ag_drop.VERT'):
    # maximum # of points: 100000
    bias_mV = 0
        
    try:
        with open(fname, 'w') as f:
            # header 
            f.write('[VZDATA]\n')
            
            # number of points, delay time between points in DSP cycles
            f.write('%d %d\n' %(ramp_points, int(tdelay_cycles)))
            
            # V (mV), z (mV)
            f.write('%d %d\n' %(bias_mV, 0))
            
            # linear ramp to crash 
            for i in range(1000):
                f.write('%d %d\n' %(bias_mV, int(tcd_mv/1000.0*i)))
            
            # stay at crash for 0.5s 
            for i in range(500):
                f.write('%d %d\n' %(bias_mV, tcd_mv))
            
            # linear ramp back to no tip height change 
            for i in range(1000):
                f.write('%d %d\n' %(bias_mV, tcd_mv/1000.0*(1000-i)))
            
            f.write('%d %d\n' %(bias_mV, 0))
            
            f.close()
    except Exception as e:
        print(e)
        print('couldnt open Ag_drop.VERT, try closing the file')
        exit(0)
        
save_vert_drop_file(tcd_mv,tdelay_cycles)

stm.vert_loadextzvfile(r"C:\Users\phys-asp-lab\Dropbox\personal-folders\Abe\drop Ag script\Ag_drop.VERT")
print('loaded vert file')

stm.setparam('VertFBMode', 4) # I(V,z)

print('performing VM')
stm.btn_vertspec(128,0)
# stm.vert_setextzv(0)

# # get the z and I data for this VM
vertdata_z = stm.vertdata(2,0) #(channel, units)
vertdata_I = stm.vertdata(3,0) 
try:
    with open(r'C:\Users\phys-asp-lab\Dropbox\personal-folders\Abe\drop Ag script\Ag_drop_%lf_angstroms.dat' %(tipcrashdist), 'w') as f:
        f.write('z\tI\n')
        for n, d in enumerate(vertdata_z):
            f.write('%lf\t%lf\n' %(d[0], vertdata_I[n][0]))
        f.close()
    plt.plot(vertdata_z)
    plt.show()
except Exception as e:
    print(e)
    print('could not open log file')
    exit(0)
stm.waitms(2000) # s

tipheight_a = stm.signals1data(2, 0.1, 1)
print('\t tip height after VM: %1.5lf angstroms' %tipheight_a)

diff_pm = abs(tipheight_b - tipheight_a)*100

print("tip apparent height changed %1.2lf pm" %(diff_pm))

if diff_pm > 20: # pm 
    print('i think something was dropped!')
