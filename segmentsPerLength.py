import PyPN
import numpy as np



for diameter in [0.1,0.5,1.,2., 3., 4.]:

    # axon definitions
    unmyelinatedParameters = {'fiberD': diameter}

    # set all properties of the bundle
    bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': 10000,  # um Axon length
                        'randomDirectionComponent': 0.0,
                        # 'bundleGuide': bundleGuide,

                        'numberOfAxons': 1,  # Number of axons in the bundle
                        'pMyel': 0,  # Percentage of myelinated fiber type A
                        'pUnmyel': 1,  # Percentage of unmyelinated fiber type C
                        'paramsMyel': unmyelinatedParameters,  # parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                        'tStop': 10,
                        'timeRes': 0.0025,

                        'saveI':True,
                        'saveV': False,
                        # 'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                        'numberOfSavedSegments': 50,
                        # number of segments of which the membrane potential is saved to disk
                        }

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**bundleParameters)

    print 'Micrometer bundle per segment: %3.1f' % (10000./float(bundle.axons[0].get_number_of_segs()))