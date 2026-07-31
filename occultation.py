#!/usr/bin/python3

# Short Python code to determine the relative occultation depths from both reflected and emitted lights, in parts per million.
# The code requires user input of stellar effective temperature, stellar radius, planetary radius, orbital semi-major axis, a wavelength range, and a few albedo choices.
# A plot is created, turning on the relevant option will save this figure. All parts of the script can be tuned as needed.
# 
# by Adam Stevenson, inspired by code from Amaury Triaud
# last modified 31/7/2026
# Queries and improvements to a.t.stevenson@bham.ac.uk or adam.stevenson29@gmail.com


import numpy as np
import matplotlib.pyplot as plt
import sys


# ----------------------  USER INPUTS ------------------------------

teff = 3000
rp = 9.58
rstar=0.23
sma=0.026
wl = np.arange(0.6e-6, 5.30e-6, 1e-9) # e.g. NIRSPEC PRISM

# define if albedo is known or not
# if unknown, the figure will have three panels. if known, there will only be a single albedo choice and one panel in the figure.
albedo_known = False

if albedo_known == False:
    # enter three exploratory values
    albedo=[0.05,0.3,0.6]
else:
   # enter the actual albedo value
   albedo = [0.1] #leave as a 1-item list or suffer the consequences

targ_name = 'Target'

# to save the figure created
save_figure = True


# ---------------------- END OF INPUTS ------------------------------


def main():
    if len(albedo)==3:
        fig,ax=plt.subplots(nrows=1,ncols=3,figsize=(12,5))

        occult1=reflected_emitted(teff,rp,rstar,sma,wl,albedo[0])
        ax[0].plot(wl*1e6,occult1[2],label='Reflected',c='dodgerblue',linewidth=3)
        ax[0].plot(wl*1e6,occult1[3],label='Emitted',c='darkred',linewidth=3)
        ax[0].plot(wl*1e6,occult1[2]+occult1[3],label='Total',c='k',linewidth=2,linestyle='--')
        ax[0].set_yscale('log')
        ax[0].legend(loc='best')

        occult2=reflected_emitted(teff,rp,rstar,sma,wl,albedo[1])
        ax[1].plot(wl*1e6,occult2[2],label='Reflected',c='dodgerblue',linewidth=3)
        ax[1].plot(wl*1e6,occult2[3],label='Emitted',c='darkred',linewidth=3)
        ax[1].plot(wl*1e6,occult2[2]+occult2[3],label='Total',c='k',linewidth=2,linestyle='--')
        ax[1].set_yscale('log')
        ax[1].legend(loc='best')

        occult3=reflected_emitted(teff,rp,rstar,sma,wl,albedo[2])
        ax[2].plot(wl*1e6,occult3[2],label='Reflected',c='dodgerblue',linewidth=3)
        ax[2].plot(wl*1e6,occult3[3],label='Emitted',c='darkred',linewidth=3)
        ax[2].plot(wl*1e6,occult3[2]+occult3[3],label='Total',c='k',linewidth=2,linestyle='--')
        ax[2].set_yscale('log')
        ax[2].legend(loc='best')


        # These are suitable for the current target. Feel free to modify!
        ax[0].set_ylim(1e0,1e3)
        ax[1].set_ylim(1e0,1e3)
        ax[2].set_ylim(1e0,1e3)


        ax[0].set_title(f'Albedo={albedo[0]}')
        ax[1].set_title(f'Albedo={albedo[1]}')
        ax[2].set_title(f'Albedo={albedo[2]}')

        fig.suptitle(f'{targ_name}',fontsize=15)
        fig.supxlabel(r'Wavelength [$\mu$m]',fontsize=15)
        fig.supylabel('Occultation depth [ppm]',fontsize=15)
        fig.tight_layout();
        if save_figure==True:
            plt.savefig(f'occulation_depths_{targ_name}.pdf',bbox_inches='tight')
        plt.show()

    elif len(albedo)==1:
        fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(4,5))

        occult=reflected_emitted(teff,rp,rstar,sma,wl,albedo[0])
        ax.plot(wl*1e6,occult[2],label='Reflected',c='dodgerblue',linewidth=3)
        ax.plot(wl*1e6,occult[3],label='Emitted',c='darkred',linewidth=3)
        ax.plot(wl*1e6,occult[2]+occult[3],label='Total',c='k',linewidth=2,linestyle='--')
        ax.set_yscale('log')
        ax.legend(loc='best')

        # These are suitable for the current target. Feel free to modify!
        ax.set_ylim(1e0,1e3)

        ax.set_title(f'Albedo={albedo[0]}')

        fig.suptitle(f'{targ_name}',fontsize=15)
        fig.supxlabel(r'Wavelength [$\mu$m]',fontsize=15)
        fig.supylabel('Occultation depth [ppm]',fontsize=15)
        fig.tight_layout()
        if save_figure==True:
            plt.savefig(f'occulation_depths_{targ_name}.pdf',bbox_inches='tight')
        plt.show()
    else:
        # 1x3 plot is fine, could iterate over whatever number of albedo inputs but too many options is overkill. 
        print('Wrong number of albedo inputs, cancelling for now')
        sys.exit()

# ----------------------  END OF MAIN ------------------------------


#-------------------------------------------------------------------
#-------------------------------------------------------------------



# ----------------------  FUNCTIONS --------------------------------



def reflected_emitted(teff,rp,rstar,sma,wave_array=np.arange(0.6e-6, 5.30e-6, 1e-9) ,albedo=0.05):
    """

    Calculates the black body curves of reflected and emitted light. Also returns the occultation depths for these two components.

    Inputs:
    teff = effective temperature of host star
    rp = radius of planet, in units of Earth radius
    rstar = radius of star, in units of solar radius
    sma = semi major axis of orbit, in AU
    wave_array = wavelength range of interest, in metres, e.g. np.arange(0.5e-6, 10e-6, 1e-9) or similar
    albedo = estimate of albedo for planet. Default set to typical HJ ball-park

    """

    # SI conversions, from astropy.constants
    Rearth = 6378100.0
    Rsun = 695700000.0
    au = 149597870700.0
    
    rp_a = (rp*Rearth)/(sma*au)
    rs_a = (rstar*Rsun)/(sma*au)
    depth = ((rp*Rearth)/(rstar*Rsun))**2

    teq = eq_temp(teff,rs_a,albedo) 

    star_BB = planck_wave(teff,wave_array)
    reflection_BB = albedo * rp_a**2 * star_BB
    emission_BB = depth * planck_wave(teq, wave_array)

    reflection_depth = (reflection_BB/star_BB) * 1e6 # multiplication factor to convert to ppm
    emission_depth = (depth * planck_wave(teq, wave_array)/star_BB) * 1e6

    return reflection_BB, emission_BB, reflection_depth, emission_depth
    

def eq_temp(teff,rs_a,albedo):
  """
  
  calculates the equilibrium temperature of an exoplanet, taking albedo into account

  Inputs:
  teff = effective temperature of host star
  rs_a = the ratio between the stellar radius and the orbital semi-major axis of the planet
  albedo = the albedo

  """

  return teff * (0.5*rs_a)**0.5 * (1-albedo)**0.25

def planck_wave(temp, wave_array):
    """
    Calculate the spectral radiance Planck function for a given wavelength range
    
    Inputs:
    temp = temperature of the object, be it star or planet
    wave_array = wavelength grid over which to sample the planck function, e.g. np.arange(0.5e-6, 10e-6, 1e-9)

    """

    # constants from astropy.constants
    h = 6.62607015e-34
    c = 299792458.0
    k_B = 1.380649e-23

    a = 2.0 * h * c**2
    b = h * c/(wave_array * k_B * temp)
    B = a / ((wave_array**5) * (np.exp(b) - 1.0))
    return B

# ----------------------  END OF FUNCTIONS ------------------------------



if __name__ == "__main__":
    main()


