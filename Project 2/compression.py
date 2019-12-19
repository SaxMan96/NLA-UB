from matplotlib import pyplot as plt
import numpy as np

def clip(img):
    img[img < 0] = 0
    img[img > 1] = 1
    return img


def compress(img, k):
    U, S, V = np.linalg.svd(img, full_matrices=False)
    img_compressed = U[:, :k] @ (np.diag(S[:k]) @ V[:k, :])
    error = np.linalg.norm(img - img_compressed, ord=2) / np.linalg.norm(img, ord=2)
    return img_compressed, error


def save_images(img, k, name):
    out_img = np.zeros(img.shape)
    total_error = 0
    for i in range(4):
        tmp_img, error = compress(img[:, :, i], k)
        out_img[:, :, i] = tmp_img
        total_error += error
    total_error = str(round(total_error / out_img.shape[2], 4))

    clipped_img = clip(out_img)
    print("k:", k, "error:", total_error)
    plt.imsave("compressed/" + name + "_" + total_error + ".png", clipped_img)
#     plt.axis("off")
#     plt.imshow(clipped_img)
#     plt.show()

for name in ["cat", "sagrada", "gal"]:
    img = plt.imread("images/" + name + ".png")
    print("Image:", name, ", Size:", img.shape[0], "x", img.shape[1])
    for k in [1, 4, 10, 30, 100]:
        save_images(img, k, name)