import numpy as np
import matplotlib.pyplot as plt


# generates gaussian function based on abscissa interval x, average mu and standard deviation sig
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


# generates equispaced gaussian functions as list of numpy arrays based on range of wavelengths lmbd
# and number of functions n_gauss
def f_gauss(lmbd, n_gauss):
    # define the standard deviation of the gaussian
    fwhm = (lmbd[lmbd.size-1] - lmbd[0]) / n_gauss
    # define full width half maximum
    sigma = fwhm / np.sqrt(2 * np.log(2))
    # create list to store data
    gauss_list = []
    for i in range(n_gauss):
        # define mean of distribution
        mu = lmbd[0] + ((2*i+1)*fwhm)/2
        # save distribution in list
        gauss_list.append(gaussian(lmbd, mu, sigma))
    return gauss_list


n = 15
waves = np.arange(400, 785, 5)
out = f_gauss(waves, n)
# plot the functions
for j in out:
    plt.plot(j)
plt.ylabel('gaussians')
plt.show()
