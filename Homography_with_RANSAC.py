# TRABALHO 2 DE VISÃO COMPUTACIONAL
# Nome: Tonim
import numpy as np
import matplotlib.pyplot as plt
import math
import cv2 as cv


########################################################################################################################
# Função para normalizar pontos
# Entrada: points (pontos da imagem a serem normalizados)
# Saída: norm_points (pontos normalizados)
#        T (matriz de normalização)
def normalize_points(points):
    # Calculate centroid
    points = points.reshape(-1, 2)  

    centroid = np.mean(points, axis=0) # Change axis to 0

    # Calculate the average distance of the points having the centroid as origin
    avg_distance = np.mean(np.linalg.norm(points - centroid, axis=1))

    # Define the scale to have the average distance as sqrt(2)
    scale = np.sqrt(2) / avg_distance

    # Define the translation to have the centroid as origin
    T = np.eye(3)
    T[0, 0] = scale
    T[1, 1] = scale
    T[0, 2] = -scale * centroid[0]
    T[1, 2] = -scale * centroid[1]

    points_hom = np.concatenate([points, np.ones((points.shape[0], 1))], axis=1)

    norm_pts = T @ points_hom.T

    return T, norm_pts

# Função para montar a matriz A do sistema de equações do DLT
# Entrada: pts1, pts2 (pontos "pts1" da primeira imagem e pontos "pts2" da segunda imagem que atendem a pts2=H.pts1)
# Saída: A (matriz com as duas ou três linhas resultantes da relação pts2 x H.pts1 = 0)
def compute_A(pts1, pts2):
    
    npoints = pts1.shape[1] #number of columns

    if pts2.shape[1] != npoints:
        raise ValueError("Number of points don't match.")



    # Essential Metrix Estimation
    # A_nx9.Es_9x1 = 0
    # Initialize matrix A

    A = np.zeros((npoints,9))

    # Stack the equation for each pair of matchings
    for i in range (npoints):

        # Create each line of matrix A
        # A[i] = [ptsl[0,i]*ptsr[0,i], ptsl[0,i]*ptsr[1,i], ptsl[0,i]*ptsr[2,i],
        #         ptsl[1,i]*ptsr[0,i], ptsl[1,i]*ptsr[1,i], ptsl[1,i]*ptsr[2,i],
        #         ptsl[2,i]*ptsr[0,i], ptsl[2,i]*ptsr[1,i], ptsl[2,i]*ptsr[2,i] ]
        # You can use Kronecker produt -> np.kron(a,b)

        A[i] = np.kron(pts1[:,i],pts2[:,i])


    return A

# Função do DLT Normalizado
# Entrada: pts1, pts2 (pontos "pts1" da primeira imagem e pontos "pts2" da segunda imagem que atendem a pts2=H.pts1)
# Saída: H (matriz de homografia estimada)
def compute_normalized_dlt(pts1, pts2):

  # Add homogeneous coordinates
  pts1 = np.hstack((pts1,np.ones((pts1.shape[0],1)))).T
  pts2 = np.hstack((pts2,np.ones((pts2.shape[0],1)))).T


  # Compute matrix A
  A = np.zeros((3*pts1.shape[1],9))

  for i in range(pts1.shape[1]):
      # Changed pts1[:, i] to pts1[:3, i] to select only the first 3 elements
      A[3*i, 3:6] = -pts2[2,i]*pts1[:3, i]
      A[3*i, 6:9] = pts2[1,i]*pts1[:3,i]

      A[3*i+1, 0:3] = -pts2[2,i]*pts1[:3,i]
      A[3*i+1, 6:9] = pts2[0,i]*pts1[:3,i]

      A[3*i+2, 0:3] = -pts2[1,i]*pts1[:3,i]
      A[3*i+2, 3:6] = pts2[0,i]*pts1[:3, i]
  # Perform SVD(A) = U.S.Vt to estimate the homography

  U, S, Vt = np.linalg.svd(A)

  H_matrix = Vt[-1,:].reshape(3,3)

  return H_matrix

def my_homography(pts1,pts2):

    # Normalize points
    T1, norm_pts1 = normalize_points(pts1)
    T2, norm_pts2 = normalize_points(pts2)

    # Perform DLT and obtain normalized matrix
    H = compute_normalized_dlt(norm_pts1.T, norm_pts2.T)

    # Denormalize the homography matrix
    H = np.dot(np.dot(np.linalg.inv(T2), H), T1)

    return H

# Função do RANSAC
# Entradas:
# pts1: pontos da primeira imagem
# pts2: pontos da segunda imagem 
# dis_threshold: limiar de distância a ser usado no RANSAC
# N: número máximo de iterações (pode ser definido dentro da função e deve ser atualizado 
#    dinamicamente de acordo com o número de inliers/outliers)
# Ninl: limiar de inliers desejado (pode ser ignorado ou não - fica como decisão de vocês)
# Saídas:
# H: homografia estimada
# pts1_in, pts2_in: conjunto de inliers dos pontos da primeira e segunda imagens


def RANSAC(pts1, pts2, dis_threshold, N, Ninl):
    
    max_inliers = 0
    best_H = None
    best_inliers = []
    n_points = pts1.shape[0]
    
    # Garantir que pts1 e pts2 tenham a forma correta para concatenação
    pts1 = pts1.reshape(-1, 2)  # Garantir que os pontos estejam no formato correto
    pts2 = pts2.reshape(-1, 2)  # Garantir que os pontos estejam no formato correto

    # Transformar os pontos para coordenadas homogêneas (adicionando a coluna de 1s)
    pts1_hom = np.hstack((pts1, np.ones((pts1.shape[0], 1))))  # Coordenadas homogêneas

    # Processo Iterativo
    for _ in range(N):
        # Sorteio aleatório de 4 amostras
        sample_indices = np.random.choice(n_points, 4, replace=False)
        sample_pts1 = pts1[sample_indices]
        sample_pts2 = pts2[sample_indices]

        H = my_homography(sample_pts1, sample_pts2)

        projected_pts = (H @ pts1_hom.T).T # Testar os pontos restantes com a homografia estimada
        projected_pts /= projected_pts[:, 2].reshape(-1, 1)  # Normalizar para coordenadas 2D

        distances = np.linalg.norm(projected_pts[:, :2] - pts2, axis=1)
        inliers = distances < dis_threshold
        num_inliers = np.sum(inliers)

        if num_inliers > max_inliers:  # Atualiza se o número de inliers for maior que o máximo encontrado até agora
            max_inliers = num_inliers
            best_H = H
            best_inliers = inliers

            if Ninl and num_inliers >= Ninl:
                break

    # Estima a homografia final usando todos os inliers
    pts1_in = pts1[best_inliers]
    pts2_in = pts2[best_inliers]
    final_H = my_homography(pts1_in, pts2_in)

    return final_H, pts1_in, pts2_in

def load_images(case_number):
    if case_number == 1:
        img1 = cv.imread('270_left.jpg', 0)
        img2 = cv.imread('270_right.jpg', 0)
        N, dis_threshold, Ninl = 300, 90, 500
        print('Par de imagem 1 carregado')
    elif case_number == 2:
        img1 = cv.imread('San_Francisco_2349.jpg', 0)
        img2 = cv.imread('San_Francisco_2348.jpg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 2 carregado')
    elif case_number == 3:
        img1 = cv.imread('outdoors01.jpg', 0)
        img2 = cv.imread('outdoors02.jpg', 0)
        N, dis_threshold, Ninl = 1000, 100, 700
        print('Par de imagem 3 carregado')
    elif case_number == 4:
        img1 = cv.imread('comicsStarWars02.jpg', 0)
        img2 = cv.imread('comicsStarWars01.jpg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 4 carregado')
    elif case_number == 5:
        img1 = cv.imread('colorLeft.png', 0)
        img2 = cv.imread('colorRight.png', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 5 carregado')
    elif case_number == 6:
        img1 = cv.imread('monitor1.jpeg', 0)
        img2 = cv.imread('monitor2.jpeg', 0)
        N, dis_threshold, Ninl = 1000, 50, 700
        print('Par de imagem 6 carregado')
    elif case_number == 7:
        img1 = cv.imread('labvisio1.jpeg', 0)
        img2 = cv.imread('labvisio4.jpeg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 7 carregado')
    elif case_number == 8:
        img1 = cv.imread('labvisio2.jpeg', 0)
        img2 = cv.imread('labvisio4.jpeg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 8 carregado')
    elif case_number == 9:
        img1 = cv.imread('labvisio3.jpeg', 0)
        img2 = cv.imread('labvisio4.jpeg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 9 carregado')
    elif case_number == 10:
        img1 = cv.imread('img_left.jpg', 0)
        img2 = cv.imread('img_right.jpg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 10 carregado')
    elif case_number == 11:
        img1 = cv.imread('001.jpg', 0)
        img2 = cv.imread('002.jpg', 0)
        N, dis_threshold, Ninl = 100, 50, 500
        print('Par de imagem 11 carregado')
    else:
        raise ValueError("Número inválido! Escolha um número entre 1 e 11.")
    
    return img1, img2, dis_threshold, N, Ninl

########################################################################################################################
# Exemplo de Teste da função de homografia usando o SIFT

MIN_MATCH_COUNT = 10
#Parte do codigo para visualizar caso a caso, abaixo eu fiz um loop para percorrer todas as imagens
'''case = 3
img1, img2 = load_images(case)

sift = cv.SIFT_create()

kp1, des1 = sift.detectAndCompute(img1, None)
kp2, des2 = sift.detectAndCompute(img2, None)


# FLANN
FLANN_INDEX_KDTREE = 1
index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
search_params = dict(checks=50)
flann = cv.FlannBasedMatcher(index_params, search_params)
matches = flann.knnMatch(des1, des2, k=2)

good = []
for m, n in matches:
    if m.distance < 0.75 * n.distance:
        good.append(m)

if len(good) > MIN_MATCH_COUNT:
    src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1, 1, 2)
    dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1, 1, 2)
    
    #################################################
    M, inliers_src, inliers_dst = RANSAC(src_pts, dst_pts, dis_threshold=50, N=100, Ninl=500)
    #################################################

    img4 = cv.warpPerspective(img1, M, (img2.shape[1], img2.shape[0])) 

else:
    print("Not enough matches are found - {}/{}".format(len(good), MIN_MATCH_COUNT))
    matchesMask = None

draw_params = dict(matchColor = (0,255,0), # draw matches in green color
                   singlePointColor = None,
                   flags = 2)
img3 = cv.drawMatches(img1, kp1, img2, kp2, good, None, **draw_params)

fig, axs = plt.subplots(2, 2, figsize=(30, 15))
fig.add_subplot(2, 2, 1)
plt.imshow(img3, 'gray')
fig.add_subplot(2, 2, 2)
plt.title('Primeira imagem')
plt.imshow(img1, 'gray')
fig.add_subplot(2, 2, 3)
plt.title('Segunda imagem')
plt.imshow(img2, 'gray')
fig.add_subplot(2, 2, 4)
plt.title('Primeira imagem após transformação')
plt.imshow(img4, 'gray')
plt.show()'''

########################################################################################################################
figs = []
for i in [3, 6]:#range(1,11):
    case = i
    img1, img2, dis_threshold, N, Ninl = load_images(case)

    sift = cv.SIFT_create()

    kp1, des1 = sift.detectAndCompute(img1, None)
    kp2, des2 = sift.detectAndCompute(img2, None)

    # FLANN
    FLANN_INDEX_KDTREE = 1
    index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
    search_params = dict(checks=50)
    flann = cv.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(des1, des2, k=2)

    good = []
    for m, n in matches:
        if m.distance < 0.75 * n.distance:
            good.append(m)

    if len(good) > MIN_MATCH_COUNT:
        src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1, 1, 2)
        dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1, 1, 2)
        
        #################################################
        M, inliers_src, inliers_dst = RANSAC(src_pts, dst_pts, dis_threshold=dis_threshold, N=N, Ninl=Ninl)
        #################################################

        img4 = cv.warpPerspective(img1, M, (img2.shape[1], img2.shape[0])) 

    else:
        print("Not enough matches are found - {}/{}".format(len(good), MIN_MATCH_COUNT))
        matchesMask = None

    draw_params = dict(matchColor = (0,255,0), # draw matches in green color
                    singlePointColor = None,
                    flags = 2)
    img3 = cv.drawMatches(img1, kp1, img2, kp2, good, None, **draw_params)

    fig = plt.figure(i, figsize=(5, 3))  # Criar janela separada
    figs.append(fig)
    fig, axs = plt.subplots(2, 2, figsize=(10, 5))
    fig.add_subplot(2, 2, 1)
    plt.title('Correlação utilizando SIFT')
    plt.imshow(img3, 'gray')
    fig.add_subplot(2, 2, 2)
    plt.title('Primeira imagem')
    plt.imshow(img1, 'gray')
    fig.add_subplot(2, 2, 3)
    plt.title('Segunda imagem')
    plt.imshow(img2, 'gray')
    fig.add_subplot(2, 2, 4)
    plt.title('Primeira imagem após transformação')
    plt.imshow(img4, 'gray')
    plt.pause(0.1)
plt.show()
