import numpy as np

x            = np.array([3,5,7,1,9,8,6,6])
y            = np.array([2,1,5,10,100,6])

index        = np.argsort(x)
sorted_x     = x[index]
sorted_index = np.searchsorted(sorted_x, y)

yindex       = np.take(index, sorted_index, mode="clip")
mask         = x[yindex] != y

result       = np.ma.array(yindex, mask=mask)

print(result)
