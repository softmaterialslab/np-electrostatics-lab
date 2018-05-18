# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:20:23 2018

@author: jc
"""
import numpy as np
import pandas as pd
from sklearn.utils import shuffle

#tanh function and its derivative
def tanhF(x,deriv=False):
	if(deriv==True):
	    return (1-(tanhF(x)*tanhF(x)))
	return np.tanh(x)

#logistic function and its derivative
def logisticF(x,deriv=False):
	if(deriv==True):
	    return (logisticF(x)*(1-logisticF(x)))
	return (1/(1+np.exp(-x)))


#Accuracy Calulation
def get_accuracy(Y,yPred):
    m = np.size(Y,0)
    n = np.size(Y,1)
    count=0
    for i in range(n):
        a= Y.item(0, i)
        b= yPred.item(0, i)
        if(a==b):
            count=count+1
    acc=count/n*100;
    print(acc)

#multilayerNeuralNetwork
def multilayerNeuralNetwork(Z,y,middlelayerSize=8,learningRate=0.05,iterations=30000):
    
    #Change the number of N in the middle layer 
    m = np.size(Z,0)
    n = np.size(Z,1)
    p = np.size(y,0)
    q = np.size(y,1)  
    b = np.ones((n,), dtype=int)
    X1 = np.concatenate((Z, b[None,:]), axis=0)

    np.random.seed(1)
    
    # randomly initialize weights with mean 0
    W1 = 2*np.random.random((middlelayerSize,m+1)) - 1
    W2 = 2*np.random.random((p,middlelayerSize+1)) - 1
    
    for i in range(iterations): 
        #forward propagation in 1st layer
        Z1 = np.dot(W1,X1)
        X2= tanhF(Z1)
        #Adding bias for 2nd layer
        X2 = np.concatenate((X2, b[None,:]), axis=0)
        #forward propagation in 2nd layer
        Z2 = np.dot(W2,X2)
        YHat = logisticF(Z2)
        #backpropagation in 2nd layer
        dB2 = (YHat-y)*logisticF(Z2,True)
        dW2 = np.dot(X2,dB2.T)
        #backpropagation in 1st layer
        nW2 = np.delete(W2, middlelayerSize, 1)
        dB1 = np.dot(nW2.T,dB2)*tanhF(Z1,True)
        dW1 = np.dot(X1,dB1.T)
        #learning parameters
        W1 = W1-(learningRate*dW1.T);
        W2 = W2-(learningRate*dW2.T);
        #Error calulated as cross entropy:
        error=y*np.log(YHat)
        cross_entropy=np.sum(-np.sum(error, axis=1),axis=0)
        Error_string = f"Iteration :  {i} cross_entropy as an error :  {cross_entropy}"
        print(Error_string)
        
    np.savez("T_weights.npz", W1=W1, W2=W2)


#Caluclate the Result
def getResult_Binary(Z):
    
    global YHatR
     
    data = np.load("bin/ml_data/T_weights.npz")
    W1 = data['W1']
    W2 = data['W2']
    
    m = np.size(Z,0)
    n = np.size(Z,1)
 
    b = np.ones((n,), dtype=int)
    X1 = np.concatenate((Z, b[None,:]), axis=0)

    #forward propagation in 1st layer
    Z1 = np.dot(W1,X1)
    X2= tanhF(Z1)
    #Adding bias for 2nd layer
    X2 = np.concatenate((X2, b[None,:]), axis=0)
    #forward propagation in 2nd layer
    Z2 = np.dot(W2,X2)
    YHat = logisticF(Z2)
    YHatR = np.around(YHat, decimals=0).astype(int)
    return YHatR
    #print(YHatR)

#Caluclate the Result
def getResult(Z):
    
    global YHat
     
    data = np.load("bin/ml_data/T_weights.npz")
    W1 = data['W1']
    W2 = data['W2']
    
    m = np.size(Z,0)
    n = np.size(Z,1)
 
    b = np.ones((n,), dtype=int)
    X1 = np.concatenate((Z, b[None,:]), axis=0)

    #forward propagation in 1st layer
    Z1 = np.dot(W1,X1)
    X2= tanhF(Z1)
    #Adding bias for 2nd layer
    X2 = np.concatenate((X2, b[None,:]), axis=0)
    #forward propagation in 2nd layer
    Z2 = np.dot(W2,X2)
    YHat = logisticF(Z2)
    
    #YHat_normed = (YHat - YHat.min(axis=1)) / (YHat.max(axis=1)-YHat.min(axis=1))
 
    return YHat

#Train multilayerNeuralNetwork
def getTrainDataSet(): 
    rawdata = pd.read_csv('bin/ml_data/data.csv', names=['e','E','nanoparticle_charge','counterion_valency','total_gridpoints','cpmd_fake_mass','fake_temperature','R','R_v','FD','output','Time'])
    rawdata=rawdata.drop(rawdata.index[0])
    rawdata = shuffle(rawdata)
    traing_data_raw=rawdata.drop(['R','R_v','FD','output','Time'], axis=1).as_matrix().astype('float64')
    output_raw = rawdata['output'].copy().as_matrix().astype('int')
    x = np.array(traing_data_raw)
    y = np.array(output_raw[:,None])
    return ( x, y )


#Train multilayerNeuralNetwork
def trainNN(middlelayerSize=8,learningRate=0.0001,iterations=30000): 
    rawdata = pd.read_csv('bin/ml_data/data.csv', names=['e','E','nanoparticle_charge','counterion_valency','total_gridpoints','cpmd_fake_mass','fake_temperature','R','R_v','FD','output','Time'])
    rawdata=rawdata.drop(rawdata.index[0])
    rawdata = shuffle(rawdata)
    traing_data_raw = rawdata.drop(['R','R_v','FD','output','Time'], axis=1).as_matrix().astype('float64')
    output_raw = rawdata['output'].copy().as_matrix().astype('int')
    x = np.array(traing_data_raw)
    y = np.array(output_raw[:,None])
    
    x_normed = x / x.max(axis=0)
    #min-max normalization
    #x_normed = (x - x.min(axis=0)) / (x.max(axis=0)-x.min(axis=0))
    
    x=x_normed
     
    multilayerNeuralNetwork(x.T,y.T,middlelayerSize,learningRate,iterations)
    
#predict accuaracy for test data file
def predict():    
    rawdata = pd.read_csv('bin/ml_data/test_data.csv', names=['e','E','nanoparticle_charge','counterion_valency','total_gridpoints','cpmd_fake_mass','fake_temperature','R','R_v','FD','output','Time'])
    rawdata=rawdata.drop(rawdata.index[0])

    #print(rawdata)

    traing_data_raw=rawdata.drop(['R','R_v', 'FD','output','Time'], axis=1).as_matrix().astype('float64')
    output_raw = rawdata['output'].copy().as_matrix().astype('int')
    x = np.array(traing_data_raw)
    y = np.array(output_raw[:,None])

    m = np.size(x,0)
    trainData = getTrainDataSet()
    
    allData = np.concatenate((x, trainData[0]), axis=0) 
    yData = np.concatenate((y, trainData[1]), axis=0)

    x_normed = allData / allData.max(axis=0)
    #print(x_normed)

    result= getResult_Binary(x_normed.T)
    
    #result2= getResult(x_normed.T)
    #print(result2)

    get_accuracy(yData.T,result)

    
    
#Generate dataset for a given 'e','E','nanoparticle_charge','counterion_valency','total_gridpoints' values
def generateDataSet(e=2,E=78.5,nanoparticle_charge=-60,counterion_valency=1, total_gridpoints=752):
    
    cpmd_fake_mass=0
    cpmd_fake_mass_step=1
    mass_size=100
    fake_temperature=0.000
    fake_temperature_step=0.001
    tempsize=100
    
    data = np.empty((0,7), float)
    
    for i in range(mass_size):
        cpmd_fake_mass = cpmd_fake_mass + cpmd_fake_mass_step;
        fake_temperatureTemp=fake_temperature
        #print(cpmd_fake_mass)
        for j in range(tempsize):
            fake_temperatureTemp = fake_temperatureTemp + fake_temperature_step;
            data = np.append(data, np.array([[e,E,nanoparticle_charge,counterion_valency,total_gridpoints,cpmd_fake_mass,fake_temperatureTemp]]), axis=0)
            
    
    #combining user data with test data
    m = np.size(data,0)
    x = getTrainDataSet()[0]    
    allData = np.concatenate((data, x), axis=0) 
    
    #np.random.shuffle(data)
    
    data_normalized = allData / allData.max(axis=0)
    #min-max normalization
    #data_normalized = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0))
    
    res=getResult(data_normalized.T)
    
    res= res[:,0:m]

    #print("Finding best parameters with NN")
    #print("selected index : " + str(res.argmax(axis=1)))
    #print("confidence score : " + str(res.max(axis=1)))
    #print("Best para set : ")
    #print(data[res.argmax(axis=1),:].flatten())
    
    '''
    print("worst para set")
    print(res.argmin(axis=1))
    print(res.min(axis=1))
    print(data[res.argmin(axis=1),:])
    #print(res)
    #print(data_normed)
    '''
    
    return data[res.argmax(axis=1),:].flatten()

