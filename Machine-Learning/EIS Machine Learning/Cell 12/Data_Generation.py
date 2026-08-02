import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import keras
import os
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf

plt.style.use('dark_background')

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)

directory = r'D:\EIS Machine Learning\Cell 12'

# train_list = ['00','10','20','40','50','60','80','90','100']
# test_list = ['30']
# train_list = ['DRT 00.csv','DRT 10.csv','DRT 20.csv','DRT 40.csv','DRT 60.csv','DRT 80.csv','DRT 90.csv']
# test_list = ['DRT 30.csv','DRT 50.csv']

# final_test = ['100']
# final_test = ['DRT 70.csv']

# Cell 1
train_list = ['0','5','10','15','20','25','30','40','50','60','65','75','80','85','90','95','100']
test_list = ['35','45','55']
final_test = ['70']

# Import Data
def get_data(name):
    path = directory + " \ ".strip() + "Cell 1 EIS at " + name  + "%" + ".txt"

    data = pd.read_csv(path, sep="\t")

    # path = os.path.join(directory,name)

    # data = pd.read_csv(path)

    data = data.to_numpy()

    freq = data[:,0]
    ZI = data[:,1]
    # ZI = data[:,0]
    ZII = data[:,2]
    # ZII = data[:,1]
    # freq = None
    # mask = -ZII > 0.0001

    # freq = freq[mask]
    # ZI = ZI[mask]
    # ZII = ZII[mask]

    return freq, ZI, ZII

# Create Shift-Invariant Transforms and Image Representation
def shift_data(ZI,ZII):
    shift_ZI = np.random.uniform(-0.03,0.03,10000)
    shift_ZII = np.random.uniform(-0.5,0.5,10000)
    images = []
    nyquist = np.array([ZI,ZII])
    images.append(nyquist)

    for a, b in zip(shift_ZI,shift_ZII):
        ZI_new = ZI + a
        ZII_new = ZII + b

        image = np.array([ZI_new,ZII_new])
        images.append(image)

    return images

# Create Dataset
x_train = []
x_test = []
final_test_data = []

for train in train_list:
    freq_train, ZI_train, ZII_train = get_data(train)
    images_train = shift_data(ZI_train,ZII_train)
    x_train.append(images_train)

for test in test_list:
    freq_test, ZI_test, ZII_test = get_data(test)
    images_test = shift_data(ZI_test,ZII_test)
    x_test.append(images_test)

for c in final_test:
    freq_c, ZI_c, ZII_c = get_data(c)
    images_c = shift_data(ZI_c,ZII_c)
    final_test_data.append(images_c)

x_train = np.stack(x_train)
x_test = np.stack(x_test)
x_final = np.stack(final_test_data)

# Save this before reshaping
n_train_images = x_train.shape[1]
n_test_images = x_test.shape[1]

x_train_flat = x_train.reshape(x_train.shape[0] * x_train.shape[1], -1)
x_test_flat = x_test.reshape(x_test.shape[0] * x_test.shape[1], -1)
final_test_flat = x_final.reshape(x_final.shape[0] * x_final.shape[1], -1)

# Fit ONLY on training data
scaler = MinMaxScaler(feature_range=(0, 1))
scaler.fit(x_train_flat)

# Transform both datasets using the same scaler
X_train = scaler.transform(x_train_flat)
X_test = scaler.transform(x_test_flat)
X_final = scaler.transform(final_test_flat)

# Reshape for CNN
X_train = X_train.reshape(-1, 51, 2, 1)
X_test = X_test.reshape(-1, 51, 2, 1)
X_final = X_final.reshape(-1, 51, 2, 1)
# X_train = X_train.reshape(-1, 510, 2, 1)
# X_test = X_test.reshape(-1, 510, 2, 1)
# X_final = X_final.reshape(-1, 510, 2, 1)

# Y_train = np.repeat(
#     [0.0, 10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 90.0, 100.0],
#     n_train_images
# )

# Y_test = np.repeat(
#     [30.0],
#     n_test_images
# )

Y_train = np.repeat(
    [0.0,5.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,65.0,75.0,80.0,85.0,90.0,95.0,100.0],
    n_train_images
)

Y_test = np.repeat(
    [35.0,45.0,55.0],
    n_test_images
)

print("X_train:", X_train.shape)
print("Y_train:", Y_train.shape)

print("X_test:", X_test.shape)
print("Y_test:", Y_test.shape)

#CNN
cnn = keras.models.Sequential()

cnn.add(
    keras.layers.Conv2D(
        32,
        kernel_size=(3, 2),
        activation='relu',
        input_shape=X_train.shape[1:]
    )
)

cnn.add(keras.layers.MaxPooling2D(pool_size=(2, 1)))

cnn.add(
    keras.layers.Conv2D(
        64,
        kernel_size=(3, 1),
        activation='relu'
    )
)

cnn.add(keras.layers.MaxPooling2D(pool_size=(2, 1)))

cnn.add(keras.layers.Flatten())

cnn.add(keras.layers.Dense(64, activation='relu'))

# Regression output
cnn.add(keras.layers.Dense(1, activation='linear'))

cnn.compile(
    loss='mse',
    optimizer='adam',
    metrics=['mae']
)

try:
    history = cnn.fit(
        X_train,
        Y_train,
        validation_data=(X_test, Y_test),
        epochs=7000,
        batch_size=2048
    )

    plt.plot(history.history['mae'])
    plt.plot(history.history['val_mae'])

    plt.title('Mean Absolute Error')
    plt.ylabel('MAE')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'])

    plt.show()

except KeyboardInterrupt:
    print("Training interrupted. Saving model...")
    cnn.save("interrupted_model.keras")

# Save model
# cnn.save('Best_model.keras')
# cnn.save('EIS-SoC_large.keras')
cnn.save('Cell1.keras')

# cnn = keras.models.load_model("Best_model.keras")

cnn.summary()

# model = keras.models.load_model('EIS-SoC.keras')
# model = keras.models.load_model('EIS-SoC_large.keras')

# final_tests = ['00','10','20','30','40','50','60','70','80','90','100']
final_tests = ['0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100']

for c in final_tests:
    final_test_data = []

    freq_c, ZI_c, ZII_c = get_data(c)
    images_c = shift_data(ZI_c,ZII_c)
    final_test_data.append(images_c)
    x_final = np.stack(final_test_data)
    final_test_flat = x_final.reshape(x_final.shape[0] * x_final.shape[1], -1)
    X_final = scaler.transform(final_test_flat)
    X_final = X_final.reshape(-1, 51, 2, 1)

    pred = cnn.predict(X_final)
    print(pred.mean())

# pred = model.predict(X_final)
# pred = cnn.predict(X_final)
# print(pred.mean())