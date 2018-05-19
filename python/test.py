# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 17:14:28 2018

@author: jc
"""

import train as nn
#traing with train data
#nn.trainNN(8,0.0001,30000)
#Prediction with test data
#nn.predict()

#Test cases
#nn.generateDataSet(e=2000,E=78.5,nanoparticle_charge=-60,counterion_valency=1, total_gridpoints=1082)

#nn.generateDataSet(e=200,E=78.5,nanoparticle_charge=-60,counterion_valency=2, total_gridpoints=1082)

#nn.generateDataSet(e=20,E=78.5,nanoparticle_charge=-60,counterion_valency=2, total_gridpoints=1082)

#nn.generateDataSet(e=2,E=78.5,nanoparticle_charge=-60,counterion_valency=1, total_gridpoints=1082)

#nn.generateDataSet(e=200,E=2,nanoparticle_charge=-60,counterion_valency=1, total_gridpoints=1082)

'''
nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=1, total_gridpoints=132)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=1, total_gridpoints=752)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=1, total_gridpoints=1272)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=2, total_gridpoints=132)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=2, total_gridpoints=752)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=2, total_gridpoints=1272)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=3, total_gridpoints=132)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=3, total_gridpoints=752)

nn.generateDataSet(e=2,E=30,nanoparticle_charge=-20,counterion_valency=3, total_gridpoints=1272)

'''

paralist = nn.generateDataSet(e=2,E=78.5,nanoparticle_charge=-30,counterion_valency=3, total_gridpoints=752)
print(paralist[5])
print(paralist[6])