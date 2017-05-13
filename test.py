import PNPy
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# --------------------------------- DEFINITION ------------------------------
# ---------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=100
dt=0.0025

# ----------------------------- axon params ---------------------------

myelinatedParameters = {'fiberD': 2}
unmyelinatedParameters = {'fiberD': 2}

segmentLengthAxon = 15
rdc = 0.2 # random direction component

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
bundleLength = 40000
nAxons = 5

# bundle guide
bundleGuide = PNPy.createGeometry.get_bundle_guide_straight(bundleLength, segmentLengthAxon)

# ------------------------ intracellular stimulation params -----------------

# parameters of signals for stimulation
rectangularSignalParamsIntra = {'amplitude': 50., #50,  # Pulse amplitude (mA)
                                'frequency': 20.,  # Frequency of the pulse (kHz)
                                'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                                'stimDur': 0.05,  # Stimulus duration (ms)
                                'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                                'delay': 0.,  # ms
                                # 'invert': True,
                                # 'timeRes': timeRes,
                                }

intraParameters = {'stimulusSignal': PNPy.signalGeneration.rectangular(**rectangularSignalParamsIntra)}

# ------------------------- extracellular stimulation params -----------------

rectangularSignalParamsExtra = {'amplitude': 3000, # Pulse amplitude (nA)
                                'frequency': 1,  # Frequency of the pulse (kHz)
                                'dutyCycle': 0.5, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                                'stimDur': 1.,  # Stimulus duration (ms)
                                'waveform': 'MONOPHASIC', # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                                'delay': 0.,  # ms
                                # 'invert': True,
                                # 'timeRes': timeRes,
                                }

elecPosStim = PNPy.createGeometry.circular_electrode(bundleGuide, positionAlongBundle=12500, radius=235,
                                                     numberOfPoles=2, poleDistance=1000)
extPotMechStim = PNPy.Extracellular.precomputedFEM(bundleGuide) # , 'oil190Inner50Endoneurium')

extraParameters = {'stimulusSignal': PNPy.signalGeneration.rectangular(**rectangularSignalParamsExtra),
                   'electrodePositions': elecPosStim,
                   'extPotMech': extPotMechStim}

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 220,
                          'positionAlongBundle': bundleLength*0.5,
                          'numberOfPoles': 1,
                          'poleDistance': 1000,
                          }

electrodePoints = PNPy.createGeometry.circular_electrode(**recordingParametersNew)

extracellularMechs = []
extracellularMechs.append(PNPy.Extracellular.homogeneous(sigma=1))
extracellularMechs.append(PNPy.Extracellular.precomputedFEM(bundleGuide))
extracellularMechs.append(PNPy.Extracellular.analytic(bundleGuide))

# ------------------------------------------------------------------------------
# --------------------------- PNPy object instantiation  -----------------------
# ------------------------------------------------------------------------------

# set all properties of the bundle
bundleParameters = {'radius': 180,  #um Radius of the bundle (match carefully to extracellular mechanism)
                    'randomDirectionComponent': rdc,
                    'bundleGuide': bundleGuide,

                    'numberOfAxons': nAxons,  # Number of axons in the bundle
                    'pMyel': .5,  # Percentage of myelinated fiber type A
                    'pUnmyel': .5,  # Percentage of unmyelinated fiber type C
                    'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                    'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                    'tStop': tStop,
                    'timeRes': dt,

                    # 'saveI':True,
                    'saveV': False,

                    # 'numberOfSavedSegments': 50,  # number of segments of which the membrane potential is saved to disk
                    }

# create the bundle with all properties of axons
bundle = PNPy.Bundle(**bundleParameters)

# spiking through a single electrical stimulation
bundle.add_excitation_mechanism(PNPy.StimIntra(**intraParameters))
bundle.add_excitation_mechanism(PNPy.StimField(**extraParameters))

# add recording electrodes. One for each extracellular medium
for extracellularMech in extracellularMechs:
    bundle.add_recording_mechanism(PNPy.RecordingMechanism(electrodePoints, extracellularMech))

# ------------------------------------------------------------------------------
# -------------------------------- PNPy calculation  ---------------------------
# ------------------------------------------------------------------------------

# run the simulation
bundle.simulate()

PNPy.save_bundle(bundle)
print 'bundle saved.'

# ------------------------------------------------------------------------------
# -------------------------------- Result Plotting  ----------------------------
# ------------------------------------------------------------------------------

t, SFAPs = bundle.get_SFAPs_from_file(0)
plt.plot(t, SFAPs)
plt.show()

bundle = None
