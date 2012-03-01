from scipy import interpolate

from MSTdata import *


# Here we specify the location of the necessary calibration files:

calib_dir = path.expanduser('~/p/mst/chers/fitting/calib_files')

# Digitizer info
v_sat = 4094.9        # Saturation of a 12-bit digitzer

# Spectrometer info
wenter = 0.40         # Entrance slit half-width (square) in Angstroms
wexit = 0.18          # Exit fiber radius in Angstroms
nxpts = 100           # Number of points to use in the trasnfer functions

class IDSII_analysis(MSTdata):

    def get_DNB_signals(self):
        
        """Read in the current, voltage, etc. from the Diagnostic
        Neutral Beam and determine when the beam was fired."""

        signal_names = ['CHERS_DNB_Voltage',
                        'CHERS_DNB_Current',
                        'CHERS_ARC_Current',
                        'CHERS_Det_Current',
                        'CHERS_Grid3_Voltage'
                        ]
        subtree = 'mraw_rs'
        self.get(signal_names, subtree)
        

        # Fix some of the calibration values and units
        DNB = self.data[subtree]
        DNB['CHERS_DNB_Voltage']['signal'] *= 10.0
        DNB['CHERS_DNB_Voltage']['units'] = 'kV'
        DNB['CHERS_DNB_Current']['signal'] *= 2.0
        DNB['CHERS_DNB_Current']['units'] = 'A'        
        DNB['CHERS_ARC_Current']['signal'] *= 100.0
        DNB['CHERS_ARC_Current']['units'] = 'A'
        DNB['CHERS_Grid3_Voltage']['signal'] *= 100.0
        DNB['CHERS_Grid3_Voltage']['units'] = 'V'

        # First, let's condition the signal so that the DNB current is
        # truely zero before the beam comes on. I'm assuming here that
        # we aren't firing the beam before the plasma.
        DNB_current = DNB['CHERS_DNB_Current']
        pre_plasma = DNB_current['time'] < 0.0
        DNB_current['signal'] -= DNB_current['signal'][pre_plasma].mean()

        # Now let's define the beam-on and beam-off times and store
        # them.
        DNB_on = DNB_current['signal'] > 0.5
        DNB['start'] = DNB_current['time'][DNB_on][0]
        DNB['stop'] = DNB_current['time'][DNB_on][-1]

    def get_idsii_raw(self):

        """Read in the raw IDSII spectrometer data."""
        subtree = 'mraw_ids'
        
        # The spectrometer has two fiber views projected onto channels
        # 1-16 and 17-32. Here we will gather both views into 2D
        # arrays.
        nodes_fiber1 = [ r'idsii_' + str(i) for i in range(1,17) ]
        nodes_fiber2 = [ r'idsii_' + str(i) for i in range(17,33) ]
        nodes = nodes_fiber1 + nodes_fiber2

        # First gather all of the signals
        self.get(nodes, subtree)

        # Now create references that collect the 
        n = len(self.data[subtree][nodes_fiber1[0]])
        fiber1_data = zeros([16,n])
        fiber2_data = zeros([16,n])

        # Create two 2D arrays for more efficient data analysis
        raw_data = self.data['mraw_ids']
        null_time = raw_data['idsii_1']['time'] < 0.0

        spec_emiss_1 = zeros([16, len(raw_data['idsii_1']['time'])])
        spec_emiss_2 = zeros_like(spec_emiss_1)
        for i in range(1,17):
            spec_emiss_1[i-1,:] = \
                raw_data['idsii_' + str(i)]['signal'] - \
                raw_data['idsii_' + str(i)]['signal'][null_time].mean()
            spec_emiss_2[i-1,:] = \
                raw_data['idsii_' + str(i+16)]['signal'] - \
                raw_data['idsii_' + str(i+16)]['signal'][null_time].mean()

        #if not self.data.has_key('idsii_analysis'):
        #    self.data['idsii_analysis'] = {}
        self.data['idsii_analysis']['spec_emiss_1'] = spec_emiss_1
        self.data['idsii_analysis']['spec_emiss_2'] = spec_emiss_2
        self.data['idsii_analysis']['time'] = raw_data['idsii_1']['time']

    def wavelength_calibration(self, plot=False):

        """This routine provides the wavelength for a given angle of
        the spectrometer grating calibrated by observed line emission
        from Cd, H, He, and Hg lamps.

        The calibration curve and data are ploted if the keyword
        'plot' is set to True."""

        analysis = self.data['idsii_analysis']

        # Calibration lines from Cd lamp
        cal_wavelength = array([3403.652, 3261.055, 3610.508, 4678.149, 
                                2836.900, 2880.767, 2980.620, 
                       #        228.8022d, 
        # Calibration line from H lamp             
                                4861.330, 4340.470,
        # Calibration line from He lamp             
                                4713.146, 4471.479, 
        # Calibration line from Hg lamp             
                                2967.280, 2536.520, 
        # Calibration line from IDS-II measurement (not sure on this one--MDN)
                                2270.1
                                ])

        cal_grating_angle = array([37.644, 39.525, 34.824, 17.398,
                                   44.877, 44.339, 43.098, 
                                   #        51.380, 
                                   13.492, 23.645, 
                                   16.691, 21.335, 
                                   43.265, 48.480, 
                                   51.580
                                   ])
        i = cal_grating_angle.argsort()

        # Next, let's create an interpolation of the data so that we
        # can determine the wavelength for a given grating angle.
        spec_cal = interpolate.interp1d(cal_grating_angle[i], 
                                        cal_wavelength[i], 
                                        kind='quadratic')
        lam_0 = spec_cal(analysis['grating_angle'])
        
        analysis['grating_calibration'] = spec_cal
        analysis['lam_0'] = lam_0
        
        if plot == True:

            angles = arange(14, 51, 0.5)
            p.plot(angles, spec_cal(angles))
            p.xlabel('Spectrometer Grating Angle [deg]')
            p.ylabel('Wavelength \u212B')

        # Next, let's find the corresponding calibration file for the
        # individual channel centroids and assign wavelengths to each
        # spectrometer channel.
        cent_fname = { 3433: 'cent_600_mar09_new.dat'}[int(lam_0)]
        with open(path.join(calib_dir, cent_fname)) as f_cent:
            centroids = f_cent.read()
        lambda_12 = array(centroids.split()).astype('float')
        analysis['wavelength_1'] = lambda_12[0:16]
        analysis['wavelength_2'] = lambda_12[16:32]

    def specify_emission_line_model(self):

        # We can determine the impurity model to use based on what
        # spectral line the spectrometer is viewing.
        line_models = { 3601: ('Al+4', 26.9815, 1),
                        3586: ('Al+3', 26.9815, 1),
                        3433: ('C+5',  12.0107, 1),
                        3367: ('N+4',  14.007,  1),
                        3210: ('Al+11',26.9815, 1),
                        3106: ('Al+13',26.9815, 1),
                        3064: ('O+5',  15.9994, 1),
                        3047: ('O+4',  15.9994, 1),
                        2981: ('B+5',  10.811,  1),
                        2974: ('O+8',  15.9994, 1),
                        }
        if line_models.has_key(int(lam_0)):
            (imp_ion, m_ion, order) = line_models[int(lam_0)]
                    

    def PMT_normalizations(self):

        # IDS PMT normalization
        # ==> set PMT value to nearest 50 V setting betwen 400 and 800 V
        pmt_voltage_1 = ((( long(abs(pmt_voltage)/50.+0.5)*50 )>400)<800)
        # ==> read in appropriate normalization file

        if pmt_volt == 600: 
            filename=path.join(calib_dir,'norm_600_jun04.dat')
        else: filename = path.join(calib_dir,'norm_apr04.dat')

        with open(filename) as f:
            pmt_norm = f.read()
        pmt_norm = array(pmt_norm.split()).astype('float')
        self.data['ids_analysis']['norm'] = pmt_norm

    def transfer_funcion():
        
        """Transfer functions that describe the wavelength sensitivity
        of each spectrometer channel."""


        # We will need the normalization from the PMTs here:
        analysis = self.data['ids_analysis']
        if not analysis.has_key('norm'):
            self.PMT_normalizations()
        pmt_norm = self.data['ids_analysis']['norm']

        # Create a prototypical transfer function as a combination of
        # tanh functions
        
        xwave1 = arange(nxpts)/(nxpts-1.0) * ( analysis['lambda_1'][15] +
                                               6.0*wenter + 2.0*
                                               wexit) - wexit - 3.0*wenter

        xwave2 = arange(nxpts)/(nxpts-1.0) * ( analysis['lambda_2'][15] +
                                               6.0*wenter + 2.0*
                                               wexit) - wexit - 3.0*wenter
        
