import numpy as np
import pandas as pd

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
    '''
        :Arguments:
            *energies*
                numpy.array values of energy from sun in J/m**2/day using
                edist with mean eavg
            *dV*
                volume of each point in the space (m**3)
        :Returns:
            *ps*
                numpy.array values of pressures after heating
    '''
    ps = energies/dV
    return ps

def altitude_pressures(space):
    '''
        :Arguments:
            *space*
                numpy.array 3-D array of zeros
        :Returns:
            *pressures*
                numpy.array values of pressure
    '''
    pressures = space
    for h in range(np.array(space.shape)[2]):
        p = 101325*np.e**(-0.000119806*h)
        pressures[:, :, h] = p
    return pressures
def calc_J_air(pressures, D, dV=dv, T=329):
    '''
        :Arguments:
            *pressures*
                numpy.array values of pressure
            *D*
                Diffusion coefficient (m**2/s)
            *dV*
                volume of each point in the space (m**3)
            *T*
                temperature (Kelvins)
        :Returns:
            *js*
                numpy array of the differences of pressures between adjcent points in pressures
    '''
    R = 8.314598
    js = np.zeros([np.array(pressures.shape)[0], np.array(pressures.shape)[1], np.array(pressures.shape)[2], 9])
    for i in range(np.array(pressures.shape)[0]):
        for j in range(np.array(pressures.shape)[1]):
            for k in range(np.array(pressures.shape)[2]):
                p = pressures[i, j, k]

                n = p*dV/(R*T)
                phi = n/dV

                try:
                    assert i-1 != -1
                    pleft = pressures[i-1, j, k]
                    n = pleft*dV/(R*T)
                    phileft = n/dV
                    dphileft = phi - phileft
                except:
                    phileft = np.nan
                    dphileft = 0

                try:
                    pright = pressures[i+1, j, k]
                    n = pright*dv/(R*T)
                    phiright = n/dv
                    dphiright = phi - phiright
                except:
                    phiright = np.nan
                    dphiright = 0

                try:
                    assert j-1 != -1
                    pfront = pressures[i, j-1, k]
                    n = pfront*dV/(R*T)
                    phifront = n/dV
                    dphifront = phi - phifront
                except:
                    phifront = np.nan
                    dphifront = 0

                try:
                    pback = pressures[i, j+1, k]
                    n = pback*dV/(R*T)
                    phiback = n/dV
                    dphiback = phi - phiback
                except:
                    phiback = np.nan
                    dphiback = 0

                try:
                    assert k-1 != -1
                    pdown = pressures[i, j, k-1]
                    n = pdown*dV/(R*T)
                    phidown = n/dV
                    dphidown = phi - phidown
                except:
                    phidown = np.nan
                    dphidown = 0

                try:
                    pup = pressures[i, j, k+1]
                    n = pup*dV/(R*T)
                    phiup = n/dV
                    dphiup = phi - phiup
                except:
                    phiup = np.nan
                    dphiup = 0

                gradphix = (dphileft - dphiright)/2 #phiright-phileft
                gradphiy = (dphiback - dphifront)/2 #phifront-phiback
                gradphiz = (dphidown - dphiup)/2 #phiup-phidown

                gradphi = np.array([dphileft, dphiright, dphiback, dphifront, dphidown, dphiup,
                                    gradphix, gradphiy, gradphix])/(dV)**(1/3)
                J = -gradphi*D
                js[i, j, k, :] = J
    return js

def pressure_diffusion(pressures, D=1.76e-5, dV=dv, T=329, dt=1):
    '''
        :Arguments:
            *presures*
                numpy.array values of pressure
             *D*
                Diffusion coefficient (m**2/s)
            *dV*
                volume of each point in the space (m**3)
            *T*
                temperature (Kelvins)
            *dt*
                time step
        :Returns:
            *pressures*
                numpy.array values of pressure after taking js into account
    '''
    js = calc_J_air(pressures, D)
    R = 8.314598
    Rarray = R*np.ones(pressures.shape)
    ns = np.divide(pressures*dV, Rarray*T)

    for i in range(np.array(pressures.shape)[0]):
        for j in range(np.array(pressures.shape)[1]):
            for k in range(np.array(pressures.shape)[2]):
                a = js[i, j, k, 0]*D*dt*np.random.normal(loc=1, scale=.2)
                b = js[i, j, k, 1]*D*dt*np.random.normal(loc=1, scale=.2)
                c = js[i, j, k, 2]*D*dt*np.random.normal(loc=1, scale=.2)
                d = js[i, j, k, 3]*D*dt*np.random.normal(loc=1, scale=.2)
                e = js[i, j, k, 4]*D*dt*np.random.normal(loc=1, scale=.2)
                f = js[i, j, k, 5]*D*dt*np.random.normal(loc=1, scale=.2)

                if a > 0:
                    ns[i, j, k] += a
                    ns[i-1, j, k] -= a

                if b > 0:
                    ns[i, j, k] += b
                    ns[i+1, j, k] -= b

                if c > 0:
                    ns[i, j, k] += c
                    ns[i, j+1, k] -= c

                if d > 0:
                    ns[i, j, k] += d
                    ns[i, j-1, k] -= d

                if e > 0:
                    ns[i, j, k] += e
                    ns[i, j, k-1] -= e

                if f > 0:
                    ns[i, j, k] += f
                    ns[i, j, k+1] -= f

    pressures = ns*R*T/dV
    return pressures

def integrate_from_sun(space, D=1.76e-5, dV=dv, tmax=20,
                       dt=1):
    '''
        :Arguements:
            *space*
                numpy.array 3-D array of zeros
            *D*
               Diffusion coefficient (m**2/s)
            *dV*
                volume of each point in the space (m**3)
            *tmax*
                maximum time the function will run
            *dt*
                time step
        :Returns:
            *pt*
                pressure after solar energy is added
    '''
    pressures = np.zeros(space.shape)
    times = np.arange(0, tmax, dt)
    pt = np.zeros([len(times), np.array(space.shape)[0],
                   np.array(space.shape)[1], (np.array(space.shape)[2])])
    ps = altitude_pressures(space)
    for t in range(len(times)):
        pt[t] = ps
        energies = sun_e(space)
        ps += energytopressure(energies, dV)
        ps = pressure_diffusion(ps, D=D, dt=dt)
    return pt

def calc_J_water(humidities, dV=dv, D=2.82e-5, T=329):
    '''
        :Argeuments:
            *humidities*
                numpy.array of humdity values
            *dv*
                volume of each point in the space (m**3)
            *D*
               Diffusion coefficient for water in air (m**2/s)
            *T*
                temperature (kelvins)
        :Returns:
            *js*
                numpy array of the differences of humdity between adjcent points
    '''
    js = np.zeros([np.array(humidities.shape)[0], np.array(humidities.shape)[1], np.array(humidities.shape)[2], 9])
    for i in range(np.array(humidities.shape)[0]):
        for j in range(np.array(humidities.shape)[1]):
            for k in range(np.array(humidities.shape)[2]):
                h = humidities[i, j, k]
                phi = h

                try:
                    assert i-1 != -1
                    hleft = humidities[i-1, j, k]
                    phileft = hleft
                    dphileft = phi - phileft
                except:
                    phileft = np.nan
                    dphileft = 0

                try:
                    hright = humidities[i+1, j, k]
                    phiright = hright
                    dphiright = phi - phiright
                except:
                    phiright = np.nan
                    dphiright = 0

                try:
                    assert j-1 != -1
                    hfront = humidities[i, j-1, k]
                    phifront = hfront
                    dphifront = phi - phifront
                except:
                    phifront = np.nan
                    dphifront = 0

                try:
                    hback = humidities[i, j+1, k]
                    phiback = hback
                    dphiback = phi - phiback
                except:
                    phiback = np.nan
                    dphiback = 0

                try:
                    assert k-1 != -1
                    hdown = humidities[i, j, k-1]
                    phidown = hdown
                    dphidown = phi - phidown
                except:
                    phidown = np.nan
                    dphidown = 0

                try:
                    hup = humidities[i, j, k+1]
                    phiup = hup
                    dphiup = phi - phiup
                except:
                    phiup = np.nan
                    dphiup = 0

                gradphix = (dphileft - dphiright)/2 #phiright-phileft
                gradphiy = (dphiback - dphifront)/2 #phifront-phiback
                gradphiz = (dphidown - dphiup)/2 #phiup-phidown

                gradphi = np.array([dphileft, dphiright, dphiback, dphifront, dphidown, dphiup,
                                    gradphix, gradphiy, gradphix])/(dV)**(1/3)
                J = -gradphi*D
                js[i, j, k, :] = J
    return js

def water_diffusion(humidities, D=2.82e-5, dV=dv, dt=1):
    '''
        :Argeuments:
            *humidities*
                numpy.array of humdity values
            *D*
               Diffusion coefficient for water in air (m**2/s)
            *dv*
                volume of each point in the space (m**3)
            *dt*
                time step
        :Returns:
            *humidities*
                numpy.array of humdity values after t seconds of diffusion
    '''
    js = calc_J_water(humidities, D=D)
    ns = humidities*dV

    for i in range(np.array(humidities.shape)[0]):
        for j in range(np.array(humidities.shape)[1]):
            for k in range(np.array(humidities.shape)[2]):
                a = js[i, j, k, 0]*dt*np.random.normal(loc=1, scale=.2)
                b = js[i, j, k, 1]*dt*np.random.normal(loc=1, scale=.2)
                c = js[i, j, k, 2]*dt*np.random.normal(loc=1, scale=.2)
                d = js[i, j, k, 3]*dt*np.random.normal(loc=1, scale=.2)
                e = js[i, j, k, 4]*dt*np.random.normal(loc=1, scale=.2)
                f = js[i, j, k, 5]*dt*np.random.normal(loc=1, scale=.2)

                if a > 0:
                    ns[i, j, k] += a
                    ns[i-1, j, k] -= a

                if b > 0:
                    ns[i, j, k] += b
                    ns[i+1, j, k] -= b

                if c > 0:
                    ns[i, j, k] += c
                    ns[i, j+1, k] -= c

                if d > 0:
                    ns[i, j, k] += d
                    ns[i, j-1, k] -= d

                if e > 0:
                    ns[i, j, k] += e
                    ns[i, j, k-1] -= e

                if f > 0:
                    ns[i, j, k] += f
                    ns[i, j, k+1] -= f

    humidities = ns/dV
    return humidities

def integrate_from_cloud(space, D=2.82e-5, loc=(10, 10, 10), size=(5, 1, 1), dV=dv, tmax=20,
                       dt=1):
    '''
        :Argeuments:
            *space*
                numpy.array 3-D array of zeros
            *D*
               Diffusion coefficient for water in air (m**2/s)
            *loc*
                int location of pressure difference;
                default = 10
            *size*
                size of the Cloud
            *dv*
                volume of each point in the space (m**3)
            *tmax*
                maximum time the function will run
            *dt*
                time step
        :Returns:
            *ht*
                humdity values after diffusion
    '''
    times = np.arange(0, tmax, dt)
    ht = np.zeros([len(times), np.array(space.shape)[0],
                   np.array(space.shape)[1], (np.array(space.shape)[2])])
    hs = cloud(space, loc=loc, size=size)
    for t in range(len(times)):
        ht[t] = hs
        hs = water_diffusion(hs, D=D, dt=dt)
    return ht
