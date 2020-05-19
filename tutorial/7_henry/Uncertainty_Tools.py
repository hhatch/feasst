import numpy as np
from scipy import stats

def deltaS_Uncertainty_Canonical(Kstar, Kstar_var, beta, input_dict, alpha, nthreads=1):
    #Code assumes the same number of samples for each Kstar[i]
    #alpha = confidence level
    samples = int(input_dict["trials"]) * nthreads
    
    stdev_Ki_mean = []
    
    #Uncertainty in individual coefficients
    for s2Ki in Kstar_var:
        # Standard deviation of the mean
        stdev_Ki_mean.append(np.sqrt(s2Ki/float(samples-1)))
    stdev_Ki_mean = np.array(stdev_Ki_mean)
    coverage_factor = stats.t.ppf(1-(1.-alpha)/2., samples-1)
    CI_Ki = stdev_Ki_mean*coverage_factor
    
    #Uncertainty in derived quantities
    #  Direct computation of DeltaS
    deltaS = np.log(Kstar[0]) - beta*Kstar[0]/Kstar[1]
    #  Estimate variance via linearized uncertainty propagation
    deltaS_var = ((1./Kstar[0] + beta*Kstar[1]/(Kstar[0]**2))**2)*Kstar_var[0] + ((beta/Kstar[0])**2)*Kstar_var[1]
    #  Estimate degrees of freedom via Welch-Satterthwaite Equation
    nu_deltaS = (deltaS_var**2)/(
        ((1./Kstar[0] + beta*Kstar[1]/(Kstar[0]**2))**4)*(Kstar_var[0]**2)/float(samples) 
        + ((beta/Kstar[0])**4)*(Kstar_var[1]**2)/float(samples)
    )
    nu_deltaS = int(nu_deltaS)
    stdev_deltaS_mean = np.sqrt(deltaS_var/float(nu_deltaS))
    coverage_factor = stats.t.ppf(1-(1.-alpha)/2., nu_deltaS-1)
    CI_deltaS= stdev_deltaS_mean*coverage_factor

    return CI_Ki, CI_deltaS
