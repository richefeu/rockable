#!/usr/bin/env python
##################################################################################################
# Sphere deposition in box
##################################################################################################
#
#Author: jm vanson
#date 01/2020

import os
from math import pi
from rollpy.PreProcessor import *
from libpy.material import UO2

xs = 0.0

diametre = [0.1, 1.]
rayon = [diam/2. for diam in diametre]
forme = ['Octahedron', 'Prisme']

lx, ly, lz = [2.5, 2.5, 11. ]
nb_particles = 50
nbsteps = 240000
dt = 1.e-5
scaling = 1.
algo = 'WP'
coeff_frottement = 0.2

#nbsteps = 100.
#dt = 1.E-10

tstart = 0                            # temps de debut
confstart = 0                         # sauvegarde reprise
fname_shape = 'Shape.shp'     
fname_simu = 'input.txt'

NbInterVerlet = 100                     # Nbre de pas de temps entre chaque maj de la liste de voisins
freq_write_res = int(nbsteps-1)     # frequence ecriture fichiers resultats 

thickness_plane = lx/100.              # Planes thickness ratio
kn = 10000

# instantiation des classes simulation et DrivingSystem
simu = Simulation()
drv = DrivingConf()
simu.periodicity = [lx, ly, 0.]

#########################################################################################
# Define elementary shapes
#########################################################################################
# les elementary shape sont definis en diametre
ePlanXY = ElementaryPlan(name='PlanXY', lx=1.2*lx, ly=1.2*ly, lz=thickness_plane) # on cree les corps elementaires
ePlanXZ = ElementaryPlan(name='PlanXZ', lx=lx, ly=thickness_plane, lz=lz)
ePlanYZ = ElementaryPlan(name='PlanYZ', lx=thickness_plane, ly=ly, lz=lz)
eSphere = ElementarySphere(name='Sphere', radius=0.5)
ePrisme = ElementaryEquilateralPrisme(name='Prisme', size = 1., radius=0.1*diametre[0]/diametre[1])#rayon[0]*0.1)
ePoly = ElementaryRegularPolyhedron(name='Octahedron', nb_vert=6, size=1., radius=0.1)#0.1)
shapeConf = ShapeConf([ePlanXY, ePlanXZ, ePlanYZ, ePoly, ePrisme, eSphere]) # on cree la config des corps elementaires
eBodies = [ePoly, ePrisme]
 
#########################################################################################
# Define simulation conf basics
#########################################################################################

simu.shape_file = fname_shape
simu.nb_particles = nb_particles
simu.tstart = tstart
simu.tmax = dt*nbsteps
simu.dt = dt
simu.nb_inter_verlet = NbInterVerlet
simu.nb_inter_conf = freq_write_res

##########################################################################################
#      Create contact laws
##########################################################################################
# lois de contact plan / sphere
simu.forceLaw = 'Hertz'
CL_plan_sphere = HertzLaw(groupId1=0, groupId2=1, knContact=kn, en2Contact=-1., ktContact=kn, muContact=coeff_frottement, krContact=kn, murContact=0.)
simu.list_force_laws.append(CL_plan_sphere)

# Loi de contact sphere / sphere
#CL_sphere_sphere = AvalanchesLaw(groupId1=1, groupId2=1, knContact=kn, en2Contact=0.01, ktContact=kn,\
#                               muContact=0.2, krContact=kn, murContact=0.0)
CL_sphere_sphere = HertzLaw(groupId1=1, groupId2=1, knContact=kn, en2Contact=-1., ktContact=kn, muContact=coeff_frottement, krContact=kn, murContact=0.)
simu.list_force_laws.append(CL_sphere_sphere)

#########################################################################################
# Create materials
#########################################################################################
matSphere = MatRigid(id=1, density=UO2.density()*1E-6)
simu.list_materials.append(matSphere)
matPlan = MatRigid(id=0, density=UO2.density()*1E-6)
simu.list_materials.append(matPlan)

###############################################################################################################
#   Define planes of the box
###############################################################################################################
idbody = 0
nbplane = 0

idbody += 1
nbplane +=1
bottomPlane = Body(id=idbody, elementaryBody=ePlanXY, idGroup=0, idCluster=0, homothety=1., position=[0., 0., -lz/2.])
simu.list_bodies.append(bottomPlane)

idbody += 1
nbplane +=1
backPlane = Body(id=idbody, elementaryBody=ePlanXZ, idGroup=0, idCluster=0,  homothety=1., position=[0., -ly/2., 0.])
simu.list_bodies.append(backPlane)

idbody += 1
nbplane +=1
frontPlane = Body(id=idbody, elementaryBody=ePlanXZ, idGroup=0, idCluster=0,  homothety=1., position=[0., ly/2., 0.])
simu.list_bodies.append(frontPlane)

idbody += 1
nbplane +=1
leftPlane = Body(id=idbody, elementaryBody=ePlanYZ, idGroup=0, idCluster=0,  homothety=1., position=[-lx/2., 0., 0.])
simu.list_bodies.append(leftPlane)

idbody += 1
nbplane +=1
rightPlane = Body(id=idbody, elementaryBody=ePlanYZ, idGroup=0, idCluster=0,  homothety=1., position=[lx/2., 0., 0.])
simu.list_bodies.append(rightPlane)

idbody += 1
nbplane +=1
topPlane = Body(id=idbody, elementaryBody=ePlanXY, idGroup=0, idCluster=0,  homothety=1., position=[0., 0., +lz/2.])
simu.list_bodies.append(topPlane)

#drv_force_top_plane = ExternalDriving(bodyid=idbody, type='_z_For_', value=-1.E-5)
#drv.append(drv_force_top_plane)

#########################################################################################
# read and compute granulometry
#########################################################################################

#nb, diamii = granulo_bimodal(radius1=diametre[0], radius2=diametre[1], x1=1. - 2./nb_particles, nb=nb_particles, type='Nombre' )
nb, diamii = granulo_bimodal(radius1=diametre[0], radius2=diametre[1], x1=xs, nb=nb_particles, type='Volume' )
nb_particles = nb

try:
   fig, ax = plt.subplots(1, sharex=True)
   val, bins, truc = ax.hist(radii, bins=30)#, normed=True)
   ax.set_xlabel('Radius')
   plt.savefig('Granulometry_number.pdf')
   np.savetxt('Granulometry_number.dat', np.array(([bins[:-1], val])).transpose())
   volii = 4./3.*pi*np.array((radii))**3.
   fig2, ax2 = plt.subplots(1, sharex=True)
   val, bins, truc = ax2.hist(radii, bins=30, weights=volii)#,  normed=True)
   ax2.set_xlabel('Radius')
   plt.savefig('Granulometry_volume.pdf')
   np.savetxt('Granulometry_volume.dat', np.array(([bins[:-1], val])).transpose())
except :
   print('Visualisation de la granulometrie non disponible!')   

####################################################
#computing theoretical box size
####################################################


vol_particles = np.sum(4./3.*pi*(diamii*scaling/2.)**3.)
compacite_theo = 0.5
lx2lz_ratio = 5.
lx_th = (vol_particles/compacite_theo/lx2lz_ratio)**(1./3.)

print('\nTaille de boite | actuelle   | conseillee |')
print('lx              | '+f'{lx:10.2f}'+' | ' +f'{lx_th:10.2f}'+' |   lx = ly ' )
print('ly              | '+f'{ly:10.2f}'+' | ' +f'{lx_th:10.2f}'+' |   lz = '+f'{lx2lz_ratio:4.1f}'+'lx ' )
print('lz              | '+f'{lz:10.2f}'+' | ' +f'{lx_th*lx2lz_ratio:10.2f}'+' |   compacite theorique = '+f'{compacite_theo:4.2f} \n' )


################################################################################################################
# Particles deposit in a box: define position of the spheres
################################################################################################################

diamii_temp =  np.sort(diamii)[::-1]
#nb_comp_particles, coor, radii2, compacity = Depot_aleatoireJM(radii_temp*scaling, [lx/2., ly/2., lz/2.], lx - 2*thickness_plane, ly-2*thickness_plane, lz-2*thickness_plane, 2000)
#nb_comp_particles, coor, radii2, ptype, compacity = Depot_PotentialEnergyGrid(radii_temp*scaling, [lx/2., ly/2., lz/2.], lx - 2*thickness_plane, ly-2*thickness_plane, lz-2*thickness_plane, part_type=[], max_compacity=1., discretisation_ratio=50.)
nb_comp_particles, coor, radii2 = Depot_aleatoire_RSA(diamii_temp*scaling/2., lx- 2*thickness_plane, ly- 2*thickness_plane, lz- 2*thickness_plane, algorithm=algo, shape="Cube", exclusionDistance=np.min(diamii_temp)*0.0025)
diamii = np.array(radii2)/scaling*2.

#simu.nb_particles = nb_comp_particles

########################################################################################################################
# Add spheres to the simulation configuration
########################################################################################################################

for i in range(nb_comp_particles):
   #idbody += 1
   #body = Body(id=idbody, elementaryBody='Sphere', idGroup=1, idCluster=0, homothety=radii[i], position=[coor[3*i]-lx/2., coor[3*i+1]-ly/2., coor[3*i+2]-lz/2.])
   if abs(diamii[i] - diametre[0]) < 1.E-10: 
      idbody += 1
      bodyelementaryBody = ePoly
      body = Body(id=idbody, elementaryBody=bodyelementaryBody, idGroup=1, idCluster=0, homothety=diamii[i], position=[coor[3*i]-lx/2.+thickness_plane, coor[3*i+1]-ly/2+thickness_plane, coor[3*i+2]-lz/2.+thickness_plane])
      body.rotateRandom()
      simu.list_bodies.append(body)
   elif abs(diamii[i] - diametre[1]) < 1.E-10:
      idbody += 1
      bodyelementaryBody = ePrisme
      body = Body(id=idbody, elementaryBody=bodyelementaryBody, idGroup=1, idCluster=0, homothety=diamii[i], position=[coor[3*i]-lx/2.+thickness_plane, coor[3*i+1]-ly/2+thickness_plane, coor[3*i+2]-lz/2.+thickness_plane])
      body.rotateRandom()
      simu.list_bodies.append(body)
   else:
      print('ERROR: rayon non autorise = '+str(diamii[i]))
      print('rayons non autorises = '+str(diametre))
      quit()
      

##############################################################################################################
#   MAJ Simulation parameters
##############################################################################################################
simu.n_driven = nbplane
simu.nb_particles = idbody
simu.D_verlet = np.max(diamii)/2.*0.5
simu.d_verlet = np.max(diamii)/2.*0.05
#simu.D_verlet = 1.
#simu.d_verlet = 0.1
#simu.gravity = [0., 0., 0.]

##############################################################################################################
# Write input file for rockable
##############################################################################################################
shapeConf.write(fname_shape)
simu.write(fname_simu)
drv.write()
