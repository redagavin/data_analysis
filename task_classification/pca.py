from load import load_data
from sklearn.decomposition import PCA
import pickle as pkl

if __name__ == '__main__':
    data = load_data("..\\Tmaze_Data\\252-1375\\2018-01-07_15-14-54\\04_tmaze1")
    pca = PCA(n_components=40)
    reduced = []
    for arr in data:
        reduced.append(pca.fit_transform(arr))

    with open("pca.pkl", "wb") as f:
        pkl.dump(reduced, f)