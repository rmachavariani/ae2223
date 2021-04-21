import numpy as np
import preprocessingv2 as v2


def actFunc1(z):

    # quadratic activation function
    result = z**2
    
    return result

def actFunc2(z):

    # cubic activation function
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

    dW1 = -np.sum(summW1,axis=0).T
    dB1 = -np.sum(summB1,axis=0).T
        
    return dW2, dB2, dW1, dB1


def makeEstimate(x,w1,w2,b1,b2,nrClasses):

    # hidden layer
    z1 = np.dot(w1,x) + b1
    a1 = actFunc1(z1)

    # output layer
    z2 = np.dot(w2,a1) + b2
    a2 = actFunc2(z2)

    maxIndex = np.argmax(a2)

    yGuess = np.zeros((nrClasses,1))
    yGuess[maxIndex] = 1

    return yGuess

    
# parameters
learnRate = 1
nodesInHidden = 5
iterations = 2
learnrate = 0.5

# get the picture data
x,y,names = v2.create_arrays(500,300)

# generate initial values
w1,b1 = generate_initial_values(x.shape[0], nodesInHidden)
w2,b2 = generate_initial_values(nodesInHidden, y.shape[0])

# iterate the weights
for i in range(iterations):

    # perform forward propagation
    a1, a2, z1, z2, cost = propagate(w1, w2, b1, b2, x, y)

    # perform backward propagation
    dW2, dB2, dW1, dB1 = backProp(y, w2, z2, b2, a2,w1,z1,b1,a1,x)

    # update weights and intercept
    w1 = w1 - learnrate * dW1
    b2 = w1 - learnrate * dB1
    w2 = w2 - learnrate * dW2
    b2 = b2 - learnrate * dB2