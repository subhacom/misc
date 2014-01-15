# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#!/usr/bin/env python
"""Simulate the elastic bouncing of some molecule inside a
compartment."""

import numpy as np
import scipy.spatial.distance as dist
import mayavi
from mayavi import mlab

# np.random.seed(1)
# We emulate movement of molecules starting with two kinds A and B
# which react to create C.  A + B -> C

particleA = {'icount': 1000, # initial count
             'color': (1,0,0),   # color for display
             'name': 'A'}
particleB  = {'icount': 500,
              'color': (0,0,1),
              'name': 'B'}
particleC  = {'icount': 0,
              'color': (0,1,0),
              'name': 'C'}

particleA['pos'] = np.zeros((particleA['icount'], 3))
particleB['pos'] = np.zeros((particleB['icount'], 3))
# Since 1 molecule each of A and B create 1 molecule of C, the maximum
# no. of C is minimum of A and B
particleC['pos'] = np.zeros((min(particleB['icount'], particleA['icount']), 3))

# Scalars for color of points <- Not used
# particleA['sc'] = np.ones(particleA['icount']) * 255
# particleB['sc'] = np.ones(particleB['icount'])
# particleC['sc'] = np.zeros((min(particleB['icount'], particleA['icount']), 3))

# Scalars for size of sphere glyphs representing molecules. Setting
# this to 0 will cause the molecule to vanish from view
particleA['ss'] = np.ones(particleA['icount'])
particleB['ss'] = np.ones(particleB['icount'])
particleC['ss'] = np.zeros(min(particleB['icount'], particleA['icount'])) 
print particleC['ss']

# Cylindrical compartment inside which all this happens
compartment = {'length': 10e-6,
               'radius': 1e-6}

# If A and B come within bindingRadius, they turn into C
bindingRadius = 0.1e-6

particleScale = 1e-7
particleSpeed = 0.05e-6
particleResolution=8

steps = 100

# Enable off-screen rendering
# mlab.options.offscreen = True

#comp = mlab.plot3d([0, compartment['length']], [0, 0], [0, 0], tube_radius=compartment['radius'], tube_sides=100, opacity=0.4)
# Display the compartment using a cylinder glyph
comp = mlab.points3d([0], [0], [0], mode='cylinder', opacity=0.4)
comp.glyph.glyph_source.glyph_source.radius = compartment['radius'] + particleScale
comp.glyph.glyph_source.glyph_source.height = compartment['length'] + 2 * particleScale
comp.glyph.glyph_source.glyph_source.center = [0, 0, 0]
comp.glyph.glyph_source.glyph_source.resolution = 20
# The cylinder should start and 0,0,0 and extend to 1,0,0
comp.glyph.glyph_source.glyph_position = 'tail' 
comp.actor.actor.orientation = [1, 0, 0]

# Assign position to A molecules - x position is one sided Gaussian
particleA['pos'][:,0] = np.abs(np.random.normal(0.0, compartment['length']/5, particleA['pos'].shape[0]))
# y is Gaussian
particleA['pos'][:,1] = np.random.normal(0.0, compartment['radius']/10, particleA['pos'].shape[0])
# z is uniform positive
particleA['pos'][:,2] = np.random.rand(particleA['pos'].shape[0])*compartment['radius']

# B starts from  the other end of the cylinder
particleB['pos'][:,0] = np.abs(compartment['length'] - 
                               np.abs(np.random.normal(0.0,
                                                       compartment['length']/5,
                                                       particleB['pos'].shape[0])))
particleB['pos'][:,1] = np.random.normal(0.0,
                                         compartment['radius']/10,
                                         particleB['pos'].shape[0])
particleB['pos'][:,2] = np.random.normal(0.0,
                                         compartment['radius']/10,
                                         particleB['pos'].shape[0])

#pts = mlab.plot3d(particleA['pos'][:,0], particleA['pos'][:,1], particleA['pos'][:,2]) # This is just cool sphere

# Now we create the 3D points
apts = mlab.points3d(particleA['pos'][:,0], particleA['pos'][:,1], particleA['pos'][:,2], 
                     particleA['ss'], scale_factor=particleScale, color=particleA['color'], resolution=particleResolution, opacity=0.7)
bpts = mlab.points3d(particleB['pos'][:,0], particleB['pos'][:,1], particleB['pos'][:,2], 
                     particleB['ss'], scale_factor=particleScale, color=particleB['color'], resolution=particleResolution, opacity=0.7)
cpts = mlab.points3d(particleC['pos'][:,0], particleC['pos'][:,1], particleC['pos'][:,2],
                     particleC['ss'], scale_factor=particleScale, color=particleC['color'], resolution=particleResolution, opacity=0.7)

frameno = 0
f = mlab.gcf() #mlab.figure(1, bgcolor=(0,0,0), fgcolor=(1,1,1), size=(1600, 1200))
# for ii in dir(f.scene): print ii
# f.scene.set_size((1600,1200))
f.scene.background = (1.0, 1.0, 1.0)
# Animate
@mlab.animate(delay=10, ui=False) 
def anim():
    global frameno
    # for ii in range(steps):
    #     print 'step',ii
    while True:
        # f.scene.render()
        # scene.reset_zoom()
        # mlab.savefig('/data/subha/molmotion/particle%05d.png' % (frameno), size=(1600,1400))
        # f.scene.save_png('/data/subha/molmotion/particle%05d.png' % (frameno))
        frameno += 1
        particleA['pos'] += np.random.normal(0, 1, particleA['pos'].shape) * particleSpeed
        # Replace positions outside the X boundary by nearest boundary position
        x1 = np.where((particleA['pos'][:,0] < 0), 0, particleA['pos'][:,0])
        x2 = np.where(x1 > compartment['length'], compartment['length'], x1)
        particleA['pos'][:,0] = x2
        # Fix y-z boundary
        yzlen = np.sqrt(particleA['pos'][:,1]**2 + particleA['pos'][:,2]**2)    # Distance from origin in yz plane
        bad = np.nonzero(yzlen > compartment['radius'])
        # if y0,z0 is the position in yz plane, and angle of the vector is
        # `a` with y axis, then y1 = R * cos(a) = R * y/r and z1 = R *
        # sin(a) = R * z/r, where R is radius of cylinder, r is length of
        # the original position vector.
        particleA['pos'][bad, 1] *= compartment['radius'] / yzlen[bad]
        particleA['pos'][bad, 2] *= compartment['radius'] / yzlen[bad]
        # We want to ignore particles which have been excluded - exclusion
        # is indicated by 'ss'=0
        deadA = np.nonzero(particleA['ss'] <= 0)[0]
        particleA['pos'][deadA,:] = -1e9 # put large value to remove chance of interaction with B

        particleB['pos'] += np.random.normal(0, 1, particleB['pos'].shape) * particleSpeed
        x1 = np.where(particleB['pos'][:,0] < 0, 0, particleB['pos'][:,0])
        x2 = np.where(x1 > compartment['length'], compartment['length'], x1)
        particleB['pos'][:,0] = x2
        yzlen = np.sqrt(particleB['pos'][:,1]**2 + particleB['pos'][:,2]**2)
        bad = np.nonzero(yzlen > compartment['radius'])
        particleB['pos'][bad,1] *= compartment['radius'] / yzlen[bad]
        particleB['pos'][bad,2] *= compartment['radius'] / yzlen[bad]
        deadB = np.nonzero(particleB['ss'] <= 0)[0]
        particleB['pos'][deadB,:] = 1e9
        particleC['pos'] += np.random.normal(0, 1, particleC['pos'].shape) * particleSpeed
        x1 = np.where(particleC['pos'][:,0] < 0, 0, particleC['pos'][:,0])
        x2 = np.where(x1 > compartment['length'], compartment['length'], x1)
        particleC['pos'][:,0] = x2
        yzlen = np.sqrt(particleC['pos'][:,1]**2 + particleC['pos'][:,2]**2)
        bad = np.nonzero(yzlen > compartment['radius'])
        particleC['pos'][bad,1] *= compartment['radius'] / yzlen[bad]
        particleC['pos'][bad,2] *= compartment['radius'] / yzlen[bad]
        # Now calculate the pairwise distances between A and B molecules
        abdist = dist.cdist(particleA['pos'], particleB['pos'], 'euclidean')
        closeA, closeB = np.nonzero(abdist < bindingRadius)
        deadC = np.nonzero(particleC['ss'] <= 0.0)[0]
        if len(closeA)  > 0:
            # print 'Step', ii
            # Multiple A can overlap with one B and vice versa, only the
            # smaller number participates in reaction
            if closeA.shape[0] < closeB.shape[0]:
                fewReacParticle = particleA
                fewReacIndex = np.unique(closeA)
            else:
                fewReacParticle = particleB
                fewReacIndex = np.unique(closeB)
            # print 'Number of close molecules:', closeA.shape[0], closeB.shape[0], np.unique(closeA).shape[0], np.unique(closeB).shape[0]
            # print 'deadC starting:', deadC[0], 'to length', len(deadC)
            # print 'deadA:', deadA.shape[0]
            # print 'deadB:', deadB.shape[0]
            # print 'Reacted:', fewReacIndex.shape[0]
            # Fill C positions sequentially - number of entries only
            # increases from 0
            particleC['pos'][deadC[0]:deadC[0]+fewReacIndex.shape[0],:] = fewReacParticle['pos'][fewReacIndex,:]
            particleA['ss'][closeA] = 0.0
            particleB['ss'][closeB] = 0.0
            particleC['ss'][deadC[0]:deadC[0]+fewReacIndex.shape[0]] = 1.0
            # print 'not-deadC:', np.nonzero(particleC['ss'] > 0)[0]
        #    print abdist[close]
        #    print 'A',particleA['ss']
        #    print 'B',particleB['ss']
        #    print 'C',particleC['ss']
        apts.mlab_source.set(x=particleA['pos'][:,0], y=particleA['pos'][:,1], z=particleA['pos'][:,2], scalars=particleA['ss'])
        bpts.mlab_source.set(x=particleB['pos'][:,0], y=particleB['pos'][:,1], z=particleB['pos'][:,2], scalars=particleB['ss'])
        cpts.mlab_source.set(x=particleC['pos'][:,0], y=particleC['pos'][:,1], z=particleC['pos'][:,2], scalars=particleC['ss'])
        yield


a = anim()
mlab.show()
