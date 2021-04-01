import numpy as np


def actFunc1(z):
    
    result = z**2
    
    return result

def actFunc2(z):
    
    result = z**3 
    
    return result

def generate_initial_values(nrInputValues,nrNodesLayer):
    
    # generate random weight and intercept is zero
    w = np.zeros((nrNodesLayer,nrInputValues)) + 0.2
    b = np.zeros((nrNodesLayer,1))

    return w, b

def propagate(w1,w2,b1,b2,X,Y):
    
    # total number of training images
    m = X.shape[1]
    
    # forward propagation first layer (hidden layer)
    # linear combination
    z1 = np.dot(w1,X) + b1
    
    # activation function
    activated1 = actFunc1(z1)
    
    
    # forward propagation output layer
    # linear combination
    z2 = np.dot(w2,activated1) + b2
    
    # activation function
    activated2 = actFunc2(z2)
    
    # calculate cost
    cost = - np.sum(y*np.log(activated2),axis=1)
    
    return activated1, activated2, z1, z2, cost


def backProp(Y,w2,z2,b2,a2,w1,z1,b1,a1,x):
    
    # perform back propagation lowest layer
    summW2 = 0 
    summB2 = 0
    for i in range(Y.shape[1]):
        # reshape Y
        yy = np.array([Y[:,i] for j in range(w2.shape[1])]).T
        zz = np.array([z2 for j in range(w2.shape[1])]).T
        aa = np.array([a2 for j in range(w2.shape[1])]).T
       
        productW2 = yy * 3*aa / zz
      
        summW2 = summW2 + productW2
        
    print(summW2.shape)
    
    dW2 = -np.sum(summW2,axis=0)
    dB2 = -np.sum(Y-a2,axis=0)
    
    
    # perform backward propagation hidden layer
    summW1 = 0
    summB1 = 0
    
    for i in range(Y.shape[1]):
        
        yy = np.array([Y[:,i] for j in range(w2.shape[1])]).T
        zz = np.array([z2 for j in range(w2.shape[1])]).T
       
        productW2 = yy * 3*w2 / zz
        
        
        xx = np.array([x[:,i].T for j in range(a1.shape[0])])
        zz = np.array([z1[:,i].T for j in range(x.shape[0])]).T
        
        productW3 = 2 * xx * zz
        productB1 = 2 * zz
        
        product = np.dot(productW2.T,productW3.T)
        productB = np.dot(productW2.T,productB1.T)
    
        
        summW1 += product
        summB1 += productB
    
    
    
    #print(summW1.shape)
    dW1 = -np.sum(summW1,axis=0).T
    dB1 = -np.sum(summB1,axis=0).T
    
    print(dW1.shape)
        
    return dW2, dB2, dW1, dB1








def propagateBackward(A2,A1,Y,z2,z1):
    
    # propagate backward from cost to second layer activated value
    # returns a matrix of size (possibilities x nrPictures)
    try:
        right = np.divide(1- Y, 1 - A2)
    except RuntimeWarning:
        right = 0
    
    try:
        left = np.divide(Y,A2)
    except RuntimeWarning:
        left = 0
        
    dA =  - (left - right)
    print(dA.shape)
    
    # sigmoid backward in second layer
    s = 1/(1+np.exp(-z2))
    dZ = dA * s * (1-s)
    
    # linear backward layer two
    m = A2.shape[1]
    dW2 = 1/m * np.dot(dZ,A2.T)
    db2 = 1/m * np.sum(dZ,axis=1,keepdims=True)
    
    # backward in layer one
    s1 = 1/(1+np.exp(-z1))
    dZ1 = dA * s1 * (1-s1)
    
    # linear backward layer one
    m1 = A1.shape[1]
    dW1 = 1/m1 * np.dot(dZ1,A1.T)
    db1 = 1/m1 *np.sum(dZ1,axis=1,keepdims=True)
    
    return dW1, dW2, db1, db2

def updateWeights(w1,w2,b1,b2,dW1,dW2,db1,db2,learnRate):
    
    w1 = w1 - learnRate * dW1
    w2 = w2 - learnRate * dW2
    b1 = b1 - learnRate * db1
    b2 = b2 - learnRate * db2
    
    return w1, w2, b1, b2
    

learnRate = 1
nodesInHidden = 10
iterations = 2
learnrate = 1

x = np.array([[0.5,0.4,0.9],
              [0.5,0.4,0.8],
              [0.5,0.4,0.7]])
y = np.array([[1,0,1],
              [0,1,0]])

w1,b1 = generate_initial_values(x.shape[0], nodesInHidden)
w2,b2 = generate_initial_values(nodesInHidden, y.shape[0])


for i in range(iterations):
    
    a1, a2, z1, z2, cost = propagate(w1, w2, b1, b2, x, y)
    dW2, dB2, dW1, dB1 = backProp(y, w2, z2, b2, a2,w1,z1,b1,a1,x)
    
    print(w1.shape)
    print(dW1.shape)
    # correct
    #w1 = w1 - learnrate * dW1
    w2 = w2 - learnrate * dW2
    b2 = b2 - learnrate * dB2
    
#W1, dW2, db1, db2 = propagateBackward(a2,a1,y,z2,z1)

#1, w2, b1, b2 = updateWeights(w1, w2, b1, b2, dW1, dW2, db1, db2, learnRate)   
    


