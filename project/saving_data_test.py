import numpy as np

# Create a 3D numpy array
my_array = np.random.rand(3, 4, 5)

# Save the array to a binary file using numpy's save() function
np.save('my_array.npy', my_array)

print("success")

# Use that checks if it's the same data and returns true check
# Load the array from the file
loaded_array = np.load('second.npy')
print(loaded_array[1])

# Check if the loaded array is the same as the original
print(np.array_equal(my_array, loaded_array))