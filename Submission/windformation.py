import numpy as np
import pandas as pd

dv = 1 #amount of space each point represents in m**3
dr = dv**(1/3)

xmax = 20 #max x values in m
ymax = 20 #max y values in m
zmax = 20 #max z values in m

space = np.zeros((xmax, ymax, zmax))

def sun_e(space, eavg=1.08e7, estd=1):
        '''
            :Arguments:
                    *space*
                                numpy.array 3-D array of zeros
                                        *eavg*
                                                    int average value of energy from sun in J/m**2/day;
                                                                default = 1.08e7 (average over Chicago)
                                                                        *estd*
                                                                                    int standard deviation of energy from sun from eavg;
                                                                                                default = 1
                                                                                                    :Returns:
                                                                                                            *energies*
                                                                                                                        numpy.array values of energy from sun in J/m**2/day using
                                                                                                                                    edist with mean eavg
                                                                                                                                        '''
                                                                                                                                            energies = np.zeros(space.shape)
                                                                                                                                                for i in range(np.array(space.shape)[0]):
                                                                                                                                                            for j in range(np.array(space.shape)[1]):
                                                                                                                                                                            e = np.random.normal(loc=eavg, scale=estd)
                                                                                                                                                                                        energies[i, j, -1] = e

                                                                                                                                                                                            return energies

def cloud(space, loc=(10,10,10), size=(10,1,1), hcloud=1, hmed=.95, hair=.69):
        '''
            :Arguments:
                    *space*
                                numpy.array 3-D array of zeros
                                        *loc*
                                                    tuple location of center of cloud;
                                                                default = (10,10,10)
                                                                        *size*
                                                                                    tuple size of cloud;
                                                                                                default = (10,1,1)
                                                                                                        *hcloud*
                                                                                                                    int humidity of cloud;
                                                                                                                                default = 1
                                                                                                                                        *hair*
                                                                                                                                                    int humidity of air;
                                                                                                                                                                default = 0
                                                                                                                                                                    :Returns:
                                                                                                                                                                            *humidities*
                                                                                                                                                                                        numpy.array values of humidity
                                                                                                                                                                                            '''
                                                                                                                                                                                                assert all([i >= 0 for i in np.array(loc)-np.array(size)/2]) and \
                                                                                                                                                                                                            all([i >= 0 for i in np.array(space.shape)-np.array(loc)+np.array(size)/2]), \
                                                                                                                                                                                                                'Cloud must be within space.'

                                                                                                                                                                                                                    humidities = space

                                                                                                                                                                                                                        humidities[:][:][:] = hair

                                                                                                                                                                                                                            n = np.array(size) % 2

                                                                                                                                                                                                                                i = np.array(loc)[0]
                                                                                                                                                                                                                                    j = np.array(loc)[1]
                                                                                                                                                                                                                                        k = np.array(loc)[2]

                                                                                                                                                                                                                                            a0 = int(i-np.array(size)[0]/2)
                                                                                                                                                                                                                                                a1 = int(i+np.array(size)[0]/2) + n[0]
                                                                                                                                                                                                                                                    b0 = int(j-np.array(size)[1]/2)
                                                                                                                                                                                                                                                        b1 = int(j+np.array(size)[1]/2) + n[1]
                                                                                                                                                                                                                                                            c0 = int(k-np.array(size)[2]/2)
                                                                                                                                                                                                                                                                c1 = int(k+np.array(size)[2]/2) + n[2]

                                                                                                                                                                                                                                                                    humidities[a0-1:a1+1, b0-1:b1+1, c0-1:c1+1] = hmed
                                                                                                                                                                                                                                                                        humidities[a0:a1, b0:b1, c0:c1] = hcloud

                                                                                                                                                                                                                                                                            return humidities

                                                                                                                                                                                                                                                                        def tornado(space, loc=10, lpress=10, rpress=1):
                                                                                                                                                                                                                                                                                '''
                                                                                                                                                                                                                                                                                    :Arguments:
                                                                                                                                                                                                                                                                                            *space*
                                                                                                                                                                                                                                                                                                        numpy.array 3-D array of zeros
                                                                                                                                                                                                                                                                                                                *loc*
                                                                                                                                                                                                                                                                                                                            int location of pressure difference;
                                                                                                                                                                                                                                                                                                                                        default = 10

                                                                                                                                                                                                                                                                                                                                            :Returns:
                                                                                                                                                                                                                                                                                                                                                        *pressures*
                                                                                                                                                                                                                                                                                                                                                                    numpy.array values of pressure
                                                                                                                                                                                                                                                                                                                                                                        '''
                                                                                                                                                                                                                                                                                                                                                                            assert 0 < loc and loc < np.array(space.shape)[0], 'Pressure difference must be within space.'

                                                                                                                                                                                                                                                                                                                                                                                pressures = altitude_pressures(space)

                                                                                                                                                                                                                                                                                                                                                                                    pressures[loc:, :, :] *= rpress

                                                                                                                                                                                                                                                                                                                                                                                        pressures[:loc, :, :] *= lpress

                                                                                                                                                                                                                                                                                                                                                                                            return pressures


                                                                                                                                                                                                                                            #wind

                                                                                                                                                                                                                                        def energytopressure(energies, dV):
                                                                                                                                                                                                                                                ps = energies/dV
                                                                                                                                                                                                                                                    return ps
                                                                                                                                                                                                                                                def altitude_pressures(space):
                                                                                                                                                                                                                                                        pressures = space
                                                                                                                                                                                                                                                            for h in range(np.array(space.shape)[2]):
                                                                                                                                                                                                                                                                        p = 101325*np.e**(-0.000119806*h)
                                                                                                                                                                                                                                                                                pressures[:, :, h] = p
                                                                                                                                                                                                                                                                                    return pressures

