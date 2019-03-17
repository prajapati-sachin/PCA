import numpy as np


X = map(float, input().split(' '))
X = np.array(X)
# X = X.dot(X.transpose())

print(X.shape)

# u, s, v = np.linalg.svd(X)
# print(u)
# print(s.shape)