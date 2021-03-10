# ResacNet
Downscaling of ocean fields by fusion of heterogeneous observations using Deep Learning algorithms
We present a method to increase low resolution satellite data fields by merging them with high
resolution data provided by additional sensors. The merging was achieved using deep learning
algorithms, in the present study Convolution Neural Networks (CNN). Experiments were conducted on
simulated ocean data in the Gulf Stream region given by the NATL60 model outputs. We merged low
resolution simulated altimeter data (100x100 km) with high resolution SST data ( 2.x2.km). We were
able to retrieve sea surface ocean currents at high resolution with a good accuracy. We then added
noise to the low resolution altimeter data. Sensitivity experiments showed that taking the SST into
account greatly increases the accuracy of the velocity retrieval. The velocity information embedded in
the transport equation modeling the SST advection is taken into account by the CNN. Incorporating in
the algorithm conservative properties such as energy or enstrophy increased the accuracy of the
method.
