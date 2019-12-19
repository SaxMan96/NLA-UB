import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def pca_covariance(x, n=None):
    if n is None:
        n = x.shape[0]
    y = x.T / np.sqrt(n - 1)
    U, S, V = np.linalg.svd(y, full_matrices=False)
    return S**2, V


def pca_correlation(x, n=None, transpose=False):
    if n is None:
        n = x.shape[0]
    std = np.std(x, axis=0)
    if transpose:
        x = x.T
    x = x / std
    return pca_covariance(x, n)


def load_data_example():
    return pd.read_csv("example.dat", header=None, sep=" ", dtype=float).values


def load_data_RCsGoff():
    data = pd.read_csv("RCsGoff.csv")
    X = data.iloc[:, 1:].values.T
    X = X - np.mean(X, axis=0)
    names = data.columns[1:].values
    return X.T, names 
    
print("-"*20,"example.dat","-"*20)
print("pca correlation")

X = load_data_example().T
s2, Vt = pca_correlation(X)

variance = s2 / np.sum(s2)
# plt.title("Variance accumulated in each of the principal components")
# plt.scatter(range(len(variance)), variance)
# plt.show()

std = np.std(Vt, axis=0)
# plt.title("Standard deviation of each of the principal components")
# plt.scatter(range(len(Vt)), std)
# plt.show()

print("\tvariance: ", *variance, sep="\n\t\t")
print("\tstd: ", *std, sep="\n\t\t")

pca_space = np.dot(Vt, X)
# plt.title("PCA Space")
# plt.scatter(pca_space[0], pca_space[1])
# plt.show()


print("pca covariance")

X = load_data_example().T
s2, Vt = pca_covariance(X, 58581)

variance = s2 / np.sum(s2)
# plt.title("Variance accumulated in each of the principal components")
# plt.scatter(range(len(variance)), variance)
# plt.show()

std = np.std(Vt, axis=0)
# plt.title("Standard deviation of each of the principal components")
# plt.scatter(range(len(std)), std)
# plt.show()

print("\tvariance: ", *variance, sep="\n\t\t")
print("\tstd: ", *std, sep="\n\t\t")

pca_space = np.dot(Vt, X)
# plt.title("PCA Space")
# plt.scatter(pca_space[0], pca_space[1])
# plt.show()

print("-"*20, "data_RCsGoff", "-"*20)

X, sample_names = load_data_RCsGoff()
s2, Vt = pca_covariance(X, 58581)

variance = s2 / np.sum(s2)
# plt.title("Variance accumulated in each of the principal components")
# plt.scatter(range(len(variance)), variance)
# plt.show()

pca_space = np.dot(Vt, X)
# plt.title("PCA Space")
# plt.scatter(pca_space[0], pca_space[1])
# plt.show()

scree_plot_treshold = 2000000
# plt.plot(s2, "ro-")
# plt.plot([scree_plot_treshold] * len(s2), "b-")
# plt.show()

print("Number of principal components needed to explain the data sets:")
print("Kraiser rule: ", sum(s2 > 1))
print("Scree plot: ", sum(s2 > scree_plot_treshold))
print("3/4 of the total variance rule: ", sum(variance > 0.75) + 1)

df = pd.DataFrame(pca_space)
df.insert(0, "name", sample_names, True)
df.insert(len(df.columns), "variance", variance, True) 
df.to_csv("output.txt", index=False, header=False)
print("-"*20, "generating output file", "-"*20)
print("-output.txt")