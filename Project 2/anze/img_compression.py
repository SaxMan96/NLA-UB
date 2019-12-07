from matplotlib import pyplot as plt
import numpy as np
import math


def compress(img, k=None):
    um, s, vmt = np.linalg.svd(img, full_matrices=False)
    if k is None:
        diff = s[:-1] - s[1:]
        idx = np.argmax(diff)
        k = idx + 1
        print(k)
    img_c = np.dot(um[:, :k], np.dot(np.diag(s[:k]), vmt[:k, :]))
    rate_of_compression = (k * (1 + img.shape[0] + img.shape[1])
                           / img.shape[0] * img.shape[1])
    rel_err = (np.linalg.norm(img - img_c, ord='fro')
               / np.linalg.norm(img, ord='fro'))
    return img_c, rel_err


def to_boundaries(img):
    img[img < 0] = 0
    img[img > 1] = 1
    return img


def save_more_images(image, k_arr, image_name):
    for k in k_arr:
        print("k:", k)
        if len(image.shape) == 2:
            c_img, error = compress(image, k)
        else:
            c_img = np.zeros(image.shape)
            error = 0
            for i in range(image.shape[2]):
                tmp_img, err = compress(image[:, :, i], k)
                c_img[:, :, i] = tmp_img
                error += err
            error /= image.shape[2]

        c_img = to_boundaries(c_img)
        if len(image.shape) == 2:
            plt.imsave('results/{}-{:.5f}.png'.format(image_name, error), c_img, cmap='gray')
        else:
            plt.imsave('results/{}-{:.5f}.png'.format(image_name, error), c_img)


image_names = ['lena', 'lena_gray', 'fingerprint']
num = 10

for name in image_names:
    print('current image:', name)
    image = plt.imread('pictures/' + name + '.png')
    st = min(image.shape[:-1])
    save_more_images(image, [math.floor((i+1)*st/num) for i in range(num)], name)
