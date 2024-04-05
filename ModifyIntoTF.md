# Modify my data into tensorflow

## INTRODUCTION
My main code is writen in fortron, so the output data are all `.dat` files.
What I want to do is modify these `.dat` files into tensorflow format, which can be used and analysised in python!
## Origincal output.dat file
Here's the contents of the output.dat file, which is the total energy of the many-body system.
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/262d3ea1-af79-4d49-a553-6cd9f3e5b914)



## codes 
Below is my code that can translate input.dat into tensorflow:
```
import numpy as np
import tensorflow as tf

import matplotlib
matplotlib.use('TkAgg')  # Set the backend

import matplotlib.pyplot as plt

# Your plotting code here


# Step 1: Read the .dat file
data = np.loadtxt('input.dat')

# Step 2: Preprocess the data if needed
# For example, you might need to split your data into input features and labels

# Assuming your data has input features in the first columns and labels in the last column
X = data[:, :-1]  # Input features
y = data[:, -1]   # Labels

# Normalize the input features if needed
# X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)

# Step 3: Use TensorFlow
# Here's a simple example of building a neural network using TensorFlow
model = tf.keras.Sequential([
    tf.keras.layers.Dense(64, activation='relu', input_shape=(X.shape[1],)),
    tf.keras.layers.Dense(32, activation='relu'),
    tf.keras.layers.Dense(1)  # Assuming you have a regression task with one output
])

model.compile(optimizer='adam', loss='mse', metrics=['mae'])

# Train the model
model.fit(X, y, epochs=10, validation_split=0.2)

# Evaluate the model
loss, mae = model.evaluate(X, y)
print("Test Loss:", loss)
print("Test MAE:", mae)

# Make predictions
predictions = model.predict(X)


import numpy as np
import matplotlib.pyplot as plt

# Assuming you have some data loaded and preprocessed
# For example, X and y as input features and labels

# Plotting the data
plt.figure(figsize=(8, 6))
plt.scatter(X, y, color='blue', label='Data points')
plt.xlabel('Input features')
plt.ylabel('Labels')
plt.title('Data Distribution')
plt.legend()
plt.grid(True)
plt.show()

plt.savefig('plot.png')

```
And the last section of this code, it plot the tensorflow:
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/162abc49-1c77-47c4-9529-d91c79ddbd80)
while below is the orinal data format plotted by xmgrace:
![image](https://github.com/ubsuny/ExactDiagonalisation-CP2P2024/assets/50903294/1ec6dffc-55e2-4be5-9b83-fac9fcc75dd8)


## Next Step
Since I already successfully modify my data into tensorflow and plot it, what I want to do next is to analysis the other output data from main code; I would need to do Wigner-grid transform and Ftransform on the data, which is in tensorflow format.
