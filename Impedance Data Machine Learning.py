import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn
import keras

plt.style.use('dark_background')

training_list = ['0.00%', '25.4%', '34.8%', '44.1%', '62.7%', 
'81.4%', '90.7%', '100%']

test_list = ['53.4%', '72.0%']

# Importing Data - Producing Image like Representation from Time/Frequency Domain Data
def get_data(name_list):
    directory = r'C:\Users\migue\Desktop\Part III\IP\Low Frequency\New Machine Learning Data'
    tag = 'Data_'
    x = []
    for i in range(len(name_list)):
        file = tag + name_list[i] + '.csv'
        path = os.path.join(directory, file)
        df = pd.read_csv(path)
        img = np.tile(df['Gamma (/gamma)'], (len(df['Time Constant (/tau)']), 1))
        x.append(img)
    x = np.array(x)
    return x

# Scaling Image Data 
def data_scaler(x, scaler):
    x_copy = x.copy().astype('float64') 
    for i, X in enumerate(x_copy[:,:,:]): 
        x_copy[i,:,:] = scaler.fit_transform(X.reshape(510*510,1)).reshape(510,510) 
    return x_copy 

x_train = get_data(training_list)
x_test = get_data(test_list)

y_train = [0.00, 25.4, 34.8, 44.1, 62.7, 81.4, 90.7, 100]
y_test = [53.4, 72.0]

scaler = sklearn.preprocessing.MinMaxScaler(feature_range=(0,1))

X_train = data_scaler(x_train, scaler)
X_test = data_scaler(x_test, scaler)

X_train = X_train.reshape(-1,510,510,1)
X_test = X_test.reshape(-1,510,510,1)

output_scaler = sklearn.preprocessing.StandardScaler()

y_train = np.array(y_train)
y_test = np.array(y_test)

output_scaler.fit(y_train.reshape(-1,1))

Y_train = output_scaler.transform(y_train.reshape(-1,1))
Y_test = output_scaler.transform(y_test.reshape(-1,1))

# Model (CNN)
def cnn(X_train, X_test, Y_train, Y_test):
    cnn = keras.models.Sequential()
    cnn.add(keras.layers.Conv2D(64, kernel_size=(3,3), activation='relu', input_shape=X_train.shape[1:]))
    cnn.add(keras.layers.MaxPooling2D(pool_size=(2,2)))
    cnn.add(keras.layers.Conv2D(128, kernel_size=(3,3), activation='relu'))
    cnn.add(keras.layers.MaxPooling2D(pool_size=(2,2)))
    cnn.add(keras.layers.Flatten())
    cnn.add(keras.layers.Dense(1))
    cnn.compile(loss='mse', optimizer='adam', metrics=['mse'])
    history = cnn.fit(X_train, Y_train, validation_data = (X_test,Y_test), epochs=100, batch_size=64) 
    return cnn, history

cnn, history = cnn(X_train, X_test, Y_train, Y_test)

# Training Behaviour
plt.figure()
plt.semilogy(history.history['loss'], color='red')
plt.semilogy(history.history['val_loss'], '--', color='red')
plt.title('Optimisation History')
plt.ylabel('Loss Function (MSE)')
plt.xlabel('Epoch')
plt.legend(['Training Data', 'Test Data'], loc='upper right')
plt.show()