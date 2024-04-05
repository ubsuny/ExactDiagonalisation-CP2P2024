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
