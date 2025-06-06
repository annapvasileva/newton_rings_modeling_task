import numpy as np
import matplotlib.pyplot as plt


EPSILON = 1e-9

VISIBLE_SPECTRUM_LOWER_LIMIT = 380.0        # нм
VISIBLE_SPECTRUM_UPPER_LIMIT = 750.0        # нм
MIN_RADIUS_OF_CURVATURE = 1.0               # м
MAX_DELTA_LAMBDA = 10.0                     # нм
MAX_DISTANCE_FROM_LENSE_CENTER = 5e-3       # м
NUM_OF_DOTS = 2000
NUM_OF_DOTS_GRID = 1000
NUM_OF_LAMBDAS = 100
I0 = 1                                      # W/m2
N_LAMBDA = 100                              # число отдельно взятых длин волн из спектра

x = np.linspace(-MAX_DISTANCE_FROM_LENSE_CENTER, MAX_DISTANCE_FROM_LENSE_CENTER, NUM_OF_DOTS_GRID)
y = np.linspace(-MAX_DISTANCE_FROM_LENSE_CENTER, MAX_DISTANCE_FROM_LENSE_CENTER, NUM_OF_DOTS_GRID)

def get_input():
    lense_radius = 0.0
    is_invalid = True
    while is_invalid:
        lense_radius = float(input("Введите радиус кривизны линзы (м): "))
        if lense_radius < MIN_RADIUS_OF_CURVATURE - EPSILON:
            print(f"Радиус должен быть не менее {MIN_RADIUS_OF_CURVATURE} м!")
        else:
            is_invalid = False
    
    lambda0 = 0.0
    is_invalid = True
    while is_invalid:
        lambda0 = float(input("Введите длину волны (нм): "))
        if (lambda0 < VISIBLE_SPECTRUM_LOWER_LIMIT - EPSILON or lambda0 > VISIBLE_SPECTRUM_UPPER_LIMIT + EPSILON):
            print(f"Длина волны должна быть в пределах видимого спектра ({VISIBLE_SPECTRUM_LOWER_LIMIT} - {VISIBLE_SPECTRUM_UPPER_LIMIT} нм)!")
        else:
            is_invalid = False
    
    is_invalid = True
    is_monochromatic = False
    while is_invalid:
        ch = input("Свет монохроматический? [y/n] ")
        if ch == 'y' or ch == 'д':
            is_monochromatic = True
            is_invalid = False
        elif ch =='n' or ch == 'н':
            is_monochromatic = False
            is_invalid = False
        else:
            print("Такого варианта ответа нет!")
    
    if is_monochromatic:
        return [lense_radius, lambda0, is_monochromatic]
    
    is_invalid = True
    delta_lambda = 0.0
    while is_invalid:
        delta_lambda = float(input("Введите ширину спектра (нм): "))
        if delta_lambda < 0.0:
            print("Ширина спектра должна быть неотрицательной!")
        elif delta_lambda > MAX_DELTA_LAMBDA + EPSILON:            
            print(f"Ширина спектра не должна превышать {MAX_DELTA_LAMBDA}, чтобы свет считался квазимонохроматическим!")
        elif lambda0 + delta_lambda / 2 > VISIBLE_SPECTRUM_UPPER_LIMIT + EPSILON or lambda0 - delta_lambda / 2 < VISIBLE_SPECTRUM_LOWER_LIMIT - EPSILON:
            print(f"Спектр должен быть в пределах видимого спектра ({VISIBLE_SPECTRUM_LOWER_LIMIT} - {VISIBLE_SPECTRUM_UPPER_LIMIT} нм)!")
        else:
            is_invalid = False
    
    return [lense_radius, lambda0, is_monochromatic, delta_lambda]

def build_intensity_graph_monochromatic(R : float, lambda0 : float):
    lambda0 *= 1e-9 # switching to meters
    r = np.linspace(-5e-3, 5e-3, NUM_OF_DOTS)
    d = R - np.sqrt(R * R - r * r)
    delta_phi = 4 * np.pi * d / lambda0 + np.pi
    I = 4 * I0 * np.cos(delta_phi / 2) ** 2
    
    plt.plot(r * 1e3, I)
    plt.xlabel("Радиус r (мм)")
    plt.ylabel("Интенсивность I(r)")
    plt.title("Кольца Ньютона (монохроматический свет)")
    plt.grid()
    plt.show()

def wavelength_to_rgb(wavelength):
    gamma = 0.8
    if 380 <= wavelength < 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif 440 <= wavelength < 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif 490 <= wavelength < 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif 510 <= wavelength < 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif 580 <= wavelength < 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif 645 <= wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R, G, B)

def build_visualization_monochromatic(R : float, lambda0 : float):
    X, Y = np.meshgrid(x, y)
    r = np.sqrt(X**2 + Y**2)
    d = R - np.sqrt(R * R - r * r)
    lambda0 *= 1e-9
    image_rgb = np.zeros((NUM_OF_DOTS_GRID, NUM_OF_DOTS_GRID, 3))

    rgb = wavelength_to_rgb(lambda0 * 1e9)
    phase = 4 * np.pi * d / lambda0 + np.pi
    I = 4 * I0 * np.cos(phase / 2)**2
    image_rgb[:, :, 0] += rgb[0] * I
    image_rgb[:, :, 1] += rgb[1] * I
    image_rgb[:, :, 2] += rgb[2] * I

    image_rgb /= image_rgb.max()
    plt.imshow(image_rgb, extent=[-5, 5, -5, 5])
    plt.title(f"Кольца Ньютона в монохроматическом свете\nλ = {lambda0 * 1e9} нм, R = {R} м")
    plt.xlabel("мм")
    plt.ylabel("мм")
    plt.show()

def build_intensity_graph_quasi_monochromatic(R : float, lambda0 : float, delta_lambda : float):
    lambda0 *= 1e-9
    delta_lambda *= 1e-9
    r = np.linspace(-MAX_DISTANCE_FROM_LENSE_CENTER, MAX_DISTANCE_FROM_LENSE_CENTER, NUM_OF_DOTS)
    d = R - np.sqrt(R * R - r * r)
    lambdas = np.linspace(lambda0 - delta_lambda/2, lambda0 + delta_lambda/2, N_LAMBDA)
    I_total = np.zeros_like(r)
    for l in lambdas:
        I_total += 4 * I0 * np.cos(2 * np.pi * d / l + np.pi/2)**2
    I_total /= N_LAMBDA

    plt.plot(r * 1e3, I_total)
    plt.xlabel("Радиус r (мм)")
    plt.ylabel("Интенсивность I(r)")
    plt.title("Кольца Ньютона (квазимонохроматический свет)")
    plt.grid()
    plt.show()

def build_visualization_quasi_monochromatic(R : float, lambda0 : float, delta_lambda : float):
    X, Y = np.meshgrid(x, y)
    r = np.sqrt(X**2 + Y**2)
    d = R - np.sqrt(R * R - r * r)
    lambda0 *= 1e-9
    delta_lambda *= 1e-9
    lambdas = np.linspace(lambda0-delta_lambda/2, lambda0+delta_lambda/2, NUM_OF_LAMBDAS)
    image_rgb = np.zeros((NUM_OF_DOTS_GRID, NUM_OF_DOTS_GRID, 3))

    for l in lambdas:
        rgb = wavelength_to_rgb(l * 1e9)
        phase = 4 * np.pi * d / l + np.pi
        I = 4 * I0 * np.cos(phase / 2)**2
        image_rgb[:, :, 0] += rgb[0] * I
        image_rgb[:, :, 1] += rgb[1] * I
        image_rgb[:, :, 2] += rgb[2] * I

    image_rgb /= image_rgb.max()
    plt.imshow(image_rgb, extent=[-5, 5, -5, 5])
    plt.title(f"Кольца Ньютона в квазимонохроматическом свете\nλ = {lambda0 * 1e9} нм, Δλ = {delta_lambda * 1e9} нм, R = {R} м")
    plt.xlabel("мм")
    plt.ylabel("мм")
    plt.show()

def monochromatic_light(R : float, lambda0 : float):
    build_intensity_graph_monochromatic(R, lambda0)
    build_visualization_monochromatic(R, lambda0)
    
def quasi_monochromatic_light(R : float, lambda0 : float, delta_lambda : float):
    build_intensity_graph_quasi_monochromatic(R, lambda0, delta_lambda)
    build_visualization_quasi_monochromatic(R, lambda0, delta_lambda)
    

def main():
    args = get_input()
    if args[2]:
        monochromatic_light(args[0], args[1])
    else:
        quasi_monochromatic_light(args[0], args[1], args[3])
    

if __name__ == '__main__':
    main()