import numpy as np
import os
import time
import nrrd
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
    for m in range(n_gauss):
        # define mean of distribution
        mu = lmbd[0] + ((2*m+1)*fwhm)/2
        # save distribution in list
        gauss_list.append(gaussian(lmbd, mu, sigma))
    return gauss_list


# define range of wavelengths under analysis
waves = np.arange(400, 785, 5)

# get time
time_start = time.time()

# cycle on images
for n_img in range(9, 61):
    # put a 0 in front of the number if it has no decimals
    if n_img < 10:
        n_img = '0{}'.format(n_img)

    # import nrrd file
    data = nrrd.read('./T{}.nrrd'.format(n_img))[0]  # Load multi-spectral image from nrrd file

    # # import multi-spectral files in list of arrays
    # for k in range(1, len(waves) + 1):
    #     if k == 1:
    #         data = np.loadtxt('../SpecTex/T{}/{}.txt'.format(n_img, k))
    #         # add third dimension in order to then append other spectral data
    #         data = np.expand_dims(data, axis=2)
    #     else:
    #         # dummy variable used to import the data and expand the dimensions
    #         dummy = np.loadtxt('../SpecTex/T{}/{}.txt'.format(n_img, k))
    #         dummy = np.expand_dims(dummy, axis=2)
    #         data = np.append(data, dummy, axis=2)
    # a = (data == data1)
    print('Finished loading data from image no.{}'.format(n_img))

    # cycle on number of channels
    for n in range(1, 9):
        # generate the gaussian sensitivity functions
        sensitivity = f_gauss(waves, n)
        # # plot the functions
        # for j in out:
        #     plt.plot(j)
        # plt.ylabel('gaussians')
        # plt.show()

        # generate empty simulated array, change third dimension so that it corresponds to number of filters
        sim_shape = list(data.shape)
        sim_shape[2] = n
        sim_shape = tuple(sim_shape)
        simulation = np.zeros(sim_shape)
        # introduce list to save results
        filtered = []
        # filter each pixel of data
        for i in range(len(sensitivity)):
            # weight each pixel of data on ith sensitivity and sum over third dimension
            filtered.append(np.dot(data, sensitivity[i]))
        # convert list to array
        filtered = np.asarray(filtered, dtype=np.float64)
        # reshape in order to get right dimensions
        filtered = filtered.reshape((data.shape[0], data.shape[1], n))
        # generate directory, if needed
        directory = './T{}'.format(n_img)
        if not os.path.exists(directory):
            os.makedirs(directory)
        # save filtered image, ith channel
        # np.savetxt('{}/{}.txt'.format(directory, i+1), filtered, delimiter=' ')
        # Create the image as a nrrd-file
        nrrd.write('{}/{}_{}.nrrd'.format(directory, n_img, n), filtered, {u'type': 'double'})

    # get how long does one iteration take
    if n_img == '01':
        print('One iteration takes approximately {}'.format(time.time() - time_start))
