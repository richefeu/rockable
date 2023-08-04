#!/usr/bin/env python
##################################################################################################
# Sphere deposition in box
##################################################################################################
#
#Author: jm vanson
#date 01/2020

import os
from rollpy.PreProcessor import *

nb_particles = 2       # on fixe le nombre de particules a generer
lx = 4.              # taille de la boite
ly = 4.
lz = 4.

tstart = 0                            # temps de debut
nbsteps = 150000                         # nbsteps
dt = 1.e-5                             # pas de temps
confstart = 0                         # sauvegarde reprise
fname_shape = 'Shape.shp'     
fname_simu = 'input.txt'

NbInterVerlet = 10                     # Nbre de pas de temps entre chaque maj de la liste de voisins
freq_write_res = int(nbsteps-1)     # frequence ecriture fichiers resultats 

thickness_plane = 0.025               # Planes thickness ratio
kn = 1E5
coeff_frottement = 0.2

# instantiation des classes simulation et DrivingSystem

simu = pre_simu.Simulation()
drv = DrivingConf()
simu.periodicity = [lx, ly, 0.]

#########################################################################################
# Define elementary shapes
#########################################################################################

ePlanXY = ElementaryPlan(name='PlanXY', lx=1., ly=1., lz=0.05) # on cree les corps elementaires
ePlanXZ = ElementaryPlan(name='PlanXZ', lx=1., ly=0.05, lz=1.)
ePlanYZ = ElementaryPlan(name='PlanYZ', lx=0.05, ly=1., lz=1.)
#eSphere = ElementarySphere(name='Sphere', radius=1.)
ePoly = ElementaryEquilateralPrisme(name='Prisme')
shapeConf = ShapeConf([ePlanXY, ePlanXZ, ePlanYZ, ePoly]) # on cree la config des corps elementaires


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
# lois de contact plan / sphere
simu.forceLaw = 'Hertz'
CL_plan_sphere = HertzLaw(groupId1=0, groupId2=1, knContact=kn, en2Contact=0.0001, ktContact=kn, muContact=coeff_frottement, krContact=kn, murContact=0.)
simu.list_force_laws.append(CL_plan_sphere)

# Loi de contact sphere / sphere
#CL_sphere_sphere = AvalanchesLaw(groupId1=1, groupId2=1, knContact=kn, en2Contact=0.01, ktContact=kn,\
#                               muContact=0.2, krContact=kn, murContact=0.0)
CL_sphere_sphere = HertzLaw(groupId1=1, groupId2=1, knContact=kn, en2Contact=0.0001, ktContact=kn, muContact=coeff_frottement, krContact=kn, murContact=0.)
simu.list_force_laws.append(CL_sphere_sphere)

#########################################################################################
# Create materials
#########################################################################################
matSphere = MatRigid(id=1, density=19300)
simu.list_materials.append(matSphere)
matPlan = MatRigid(id=0, density=19300)
simu.list_materials.append(matPlan)
                                                                    

###############################################################################################################
#   Define planes of the box
###############################################################################################################
idbody = -1
nbplane = 0

idbody += 1
nbplane += 1
bottomPlane = Body(id=idbody, elementaryBody=ePlanXY, idGroup=0, idCluster=0, homothety=lx, position=[0., 0., -lz/2.-lx*thickness_plane])
simu.list_bodies.append(bottomPlane)

#idbody +=1
#nbplane += 1
#backPlane = Body(id=idbody, elementaryName='PlanXZ', idGroup=0, idCluster=0,  homothety=lx, position=[0., -ly/2.-lx*thickness_ratio, -lz/2.-lx*thickness_ratio+lx/2.])
#simu.list_bodies.append(backPlane)
#
#idbody += 1
#nbplane += 1
#frontPlane = Body(id=idbody, elementaryName='PlanXZ', idGroup=0, idCluster=0,  homothety=lx, position=[0., ly/2.+lx*thickness_ratio, -lz/2.-lx*thickness_ratio+lx/2.])
#simu.list_bodies.append(frontPlane)
#
#idbody += 1
#nbplane += 1
#leftPlane = Body(id=idbody, elementaryName='PlanYZ', idGroup=0, idCluster=0,  homothety=lx, position=[-lx/2.-lx*thickness_ratio, 0.0, -lz/2.-lx*thickness_ratio+lx/2.])
#simu.list_bodies.append(leftPlane)
#
#idbody += 1
#nbplane += 1
#rightPlane = Body(id=idbody, elementaryName='PlanYZ', idGroup=0, idCluster=0,  homothety=lx, position=[lx/2.+lx*thickness_ratio, 0.0, -lz/2.-lx*thickness_ratio+lx/2.])
#simu.list_bodies.append(rightPlane)
#
#idbody += 1
#nbplane += 1
#topPlane = Body(id=idbody, elementaryName='PlanXY', idGroup=0, idCluster=0,  homothety=lx, position=[0., 0., +lz/2.+lx*thickness_ratio])
#simu.list_bodies.append(topPlane)
#
#drv_force_top_plane = ExternalDriving(bodyid=idbody, type='_z_For_', value=-1.E-5)
#drv.append(drv_force_top_plane)

#########################################################################################
# read and compute granulometry
#########################################################################################

nb, radii = granulo_monomodal(radius=0.5, nb=nb_particles)
#nb, radii = granulo_bimodal(radius1=1., radius2=2., x1=0.5, nb=100, type='Volume')
#nb, radii = granulo_random(radius_min=-1., radius_max=1., nb=1)
#nb, radii = granulo_from_file(fname='granulometry.txt', nb=nb_particles, nbins=30, tronque_min=5, tronque_max=3, visu=0)

################################################################################################################
# Particles deposit in a box: define position of the spheres
################################################################################################################

#nb, coor, radii, compacity = Depot_aleatoireJM(radii, [lx/2., ly/2., lz/2.], lx, ly, lz, 2000)
scaling = 1.
radii_temp =  np.sort(radii)[::-1]
#nb_comp_particles, coor, radii2, compacity = Depot_aleatoireJM(radii_temp*scaling, [lx/2., ly/2., lz/2.], lx - 2*thickness_plane, ly-2*thickness_plane, lz-2*thickness_plane, 2000)
#nb_comp_particles, coor, radii2, ptype, compacity = Depot_PotentialEnergyGrid(radii_temp*scaling, [lx/2., ly/2., lz/2.], lx - 2*thickness_plane, ly-2*thickness_plane, lz-2*thickness_plane, part_type=[], max_compacity=1., discretisation_ratio=50.)
nb_comp_particles, coor, radii2 = Depot_aleatoire_RSA(radii_temp*scaling, lx- 2*thickness_plane, ly- 2*thickness_plane, lz- 2*thickness_plane, shape="Cube")
radii = np.array(radii2)/scaling
simu.nb_particles = nb

########################################################################################################################
# Add spheres to the simulation configuration
########################################################################################################################

for i in range(simu.nb_particles):
   idbody += 1
   body = Body(id=idbody, elementaryBody=ePoly, idGroup=1, idCluster=0, homothety=radii[i], position=[coor[3*i]-lx/2., coor[3*i+1]-ly/2., coor[3*i+2]-lz/2.])
   body.rotateRandom()
   simu.list_bodies.append(body)
   
##############################################################################################################
#   MAJ Simulation parameters
##############################################################################################################
simu.n_driven = nbplane
simu.nb_particles += nbplane
simu.D_verlet = np.min(radii)*0.1
simu.d_verlet = np.min(radii)*0.05

##############################################################################################################
# Write input file for rockable
##############################################################################################################
shapeConf.write(fname_shape)
simu.write(fname_simu)
drv.write()