import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, CheckButtons, TextBox


def update(val):

    global prev_t_max_a, prev_M_AN, prev_delta, prev_num_vueltas, prev_einstenian, prev_vy, prev_vx, prev_x0
    global line, y0

    if float(t_max_textbox.text) != prev_t_max_a:
        prev_t_max_a = float(t_max_textbox.text)
        change = True
    elif float(M_AN_textbox.text) != prev_M_AN:
        prev_M_AN = float(M_AN_textbox.text)
        change = True
    elif float(delta_textbox.text) != prev_delta:
        prev_delta = float(delta_textbox.text)
        change = True
    elif float(vuelta_textbox.text) != prev_num_vueltas:
        prev_num_vueltas = float(vuelta_textbox.text)
        change = True
    elif float(x0_textbox.text) != prev_x0:
        prev_x0 = float(x0_textbox.text)
        change = True
    elif checkbox_e.get_status()[0] != prev_einstenian:
        prev_einstenian = checkbox_e.get_status()[0]
        change = True
    elif vx_slider.val != prev_vx:
        prev_vx = vx_slider.val
        change = True
    elif vy_slider.val != prev_vy:
        prev_vy = vy_slider.val
        change = True
    else:
        change = False

    # Solo ejecuto si hay cambio
    if change:

        resultado = orbit_data(float(M_AN_textbox.text), float(x0_textbox.text), y0, vx_slider.val, vy_slider.val, float(vuelta_textbox.text),
                               float(t_max_textbox.text), float(delta_textbox.text), checkbox_e.get_status()[0], freq_escritura)

        updated_x = resultado[0]
        updated_y = resultado[1]
        updated_z = resultado[2]

        line.set_ydata(updated_y)
        line.set_xdata(updated_x)
        line.set_3d_properties(updated_z)

        margin = 1.1

        ax.set_xlim(min(0, np.min(line.get_xdata()), int(x0)) * margin, max(0, np.max(line.get_xdata()), int(x0)) * margin)
        ax.set_ylim(min(0, np.min(line.get_ydata()), int(y0)) * margin, max(0, np.max(line.get_ydata()), int(y0)) * margin)

        fig.canvas.draw_idle()

        if checkbox_e.get_status()[0]:
            print("Calculo einsteniano (update)")
        elif not checkbox_e.get_status()[0]:
            print("cálculo newtoniano (update)")


def orbit_data(masa_an, x_ini, y_ini, vx_ini, vy_ini, lap_num, max_t, dt, ein, fre_escri):

    max_t = max_t * 365 * 24 * 60 * 60      # Pasamos max_t a segundos

    # Inicializar vectores para almacenar
    lap = 0
    num_pasos = int(max_t / dt)
    last_value = num_pasos - 1

    masa_an = masa_an * M_sol

    pos_x = np.zeros(num_pasos)
    pos_y = np.zeros(num_pasos)
    pos_z = np.zeros(num_pasos)

    acc_x = np.zeros(num_pasos)
    acc_y = np.zeros(num_pasos)

    vel_x = np.zeros(num_pasos)
    vel_y = np.zeros(num_pasos)

    angulo = np.zeros(num_pasos)
    angulo_grad = np.zeros(num_pasos)

    vel_x[0] = vx_ini
    vel_y[0] = vy_ini
    pos_x[0] = x_ini
    pos_y[0] = y_ini

    tiempo = np.arange(0, max_t, dt)

    radio_an = 2 * G * masa_an / c ** 2 # Radio schwarzschild de agujero negro

    # Simular trayectoria planeta
    for i in range(num_pasos - 1):

        # Calcular el radio (distancia entre agujero negro y planeta)
        r = np.sqrt((pos_x[i] - xan0) ** 2 + (pos_y[i] - yan0) ** 2)

        # Cálculo de ángulo
        angulo[i] = math.atan2(pos_y[i], pos_x[i])
        if angulo[i] < 0:
            angulo[i] = 2*math.pi+angulo[i]
        angulo_grad[i] = math.degrees(angulo[i])
        # Cálculo aceleración
        if ein:
            # Angular momentum
            # From https://en.wikipedia.org/wiki/Angular_momentum. Here shoulb be per unit mass
            l_momentum = (vel_y[i] * pos_x[i] - vel_x[i] * pos_y[i])
            # From https://en.wikipedia.org/wiki/Two-body_problem_in_general_relativity
            # We found potential and d/dr to find acceleration
            acc_r1 = G * masa_an / r ** 2
            acc_r2 = - l_momentum ** 2 / r ** 3
            acc_r3 = 3 * (G * masa_an * l_momentum ** 2) / (c ** 2 * r ** 4)
            acc = acc_r1 + acc_r2 + acc_r3

            acc_x[i] = - acc * math.cos(angulo[i])
            acc_y[i] = - acc * math.sin(angulo[i])
        elif not ein:
            acc_x[i] = -((G * masa_an) / r ** 2) * math.cos(angulo[i])
            acc_y[i] = -((G * masa_an) / r ** 2) * math.sin(angulo[i])

        # Cálculo velocidad
        vel_x[i + 1] = vel_x[i] + acc_x[i] * dt
        vel_y[i + 1] = vel_y[i] + acc_y[i] * dt

        # Cálculo posición
        pos_x[i + 1] = pos_x[i] + vel_x[i] * dt
        pos_y[i + 1] = pos_y[i] + vel_y[i] * dt

        if not i % fre_escri:
            print(f"Iteración {i} de {num_pasos}. Ángulo actual {angulo_grad[i]:.2f}", "deg, ",
                  f"Posición {pos_x[i]}, {pos_y[i]}", end='\r')

        if angulo_grad[i] < angulo_grad[i - 1]:
            print(f"\nVuelta completada {lap}. Var.alfa_deg:", angulo_grad[i])
            lap = lap + 1
            if lap == lap_num:
                last_value = i
                break

    pos_x = pos_x[:last_value]
    pos_y = pos_y[:last_value]
    pos_z = pos_z[:last_value]

    acc_x = acc_x[:last_value]
    acc_y = acc_y[:last_value]

    vel_x = vel_x[:last_value]
    vel_y = vel_y[:last_value]

    angulo = angulo[:last_value]
    angulo_grad = angulo_grad[:last_value]

    tiempo = tiempo[:last_value]

    print("\nCálculo completado (orbit_data)")

    return pos_x, pos_y, pos_z, acc_x, acc_y, vel_x, vel_y, angulo, angulo_grad, tiempo


if __name__ == '__main__':

    # Zona de parámetros
    t_max_a = 10                  # Duración del cálculo en años
    delta = 1000            # Longitud del step en segundos
    x0 = 1e11               # Posicion inicial planeta en eje x
    y0 = 0                  # Posicion inicial planeta en eje y
    vx0 = 0                 # Velocidad inicial planeta en eje x
    vy0 = 1e5               # Velocidad inicial planeta en eje y
    xan0 = 0                # Posicion inicial agujero negro en eje x
    yan0 = 0                # Posicion inicial agujero negro en eje y
    num_vueltas = 10        # Número de vueltas que tendrá el calculo
    einstenian = False
    freq_escritura = 100    # Cada cuanto se escribe la iteración

    # Parámetros de física
    G = 6.67430e-11     # Constante gravitatoria universal [N*m^2/kg^2]
    c = 3e8             # Velocidad de la luz [m/s]
    M_sol = 1.989e30    # Masa del Sol [Kg]
    M_AN = 10

    # Variables que van a checkboxes
    prev_num_vueltas = num_vueltas
    prev_delta = delta
    prev_t_max_a = t_max_a
    prev_M_AN = M_AN
    prev_einstenian = einstenian
    prev_vy = vy0
    prev_vx = vx0
    prev_x0 = x0

    x, y, z, vx, vy, accx, accy, alfa, alfa_deg, t_vec = orbit_data(M_AN, x0, y0, vx0, vy0, num_vueltas, t_max_a, delta, einstenian, freq_escritura)

    # Gráfico

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(projection='3d')  # Creación del subplot en 3d
    ax.set_title('Gráfico de la orbita en 3D')
    line, = ax.plot(x, y, z, lw=2)
    ax.view_init(elev=15, azim=45)
    ax.scatter(x0, y0, 0, color='magenta', s=100)
    ax.scatter(xan0, yan0, 0, color='black', s=200)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')


    # Zona de sliders
    # Make a vertically oriented slider to control the amplitude
    vy_ax = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
    vy_slider = Slider(
        ax=vy_ax,
        label="Y speed [m/s]",
        valmin=-5e5,
        valmax=5e5,
        valinit=vy0,
        orientation="vertical",
        valstep=1e4
    )

    # Make a horizontal slider to control the frequency.
    vx_ax = fig.add_axes([0.2, 0.075, 0.65, 0.03])
    vx_slider = Slider(
        ax=vx_ax,
        label='X speed [m/s]',
        valmin=-1e5,
        valmax=1e5,
        valinit=vx0,
        orientation="horizontal",
        valstep=1e4
    )

    # Create a checkbox to add einstenian acceleration
    checkbox_ax = fig.add_axes([0.8, 0.9, 0.15, 0.04])  # Adjust position as needed
    checkbox_e = CheckButtons(checkbox_ax, ['Einstenian'], [einstenian])

    # Creación de Textboxes

    delta_ax_txt = fig.add_axes([0.8, 0.85, 0.15, 0.04])
    delta_textbox = TextBox(delta_ax_txt, "delta[s/sample]:", delta, textalignment="center")

    t_max_ax_txt = fig.add_axes([0.8, 0.80, 0.15, 0.04])
    t_max_textbox = TextBox(t_max_ax_txt, "t_max[years]:", t_max_a, textalignment="center")

    vuelta_ax_txt = fig.add_axes([0.8, 0.75, 0.15, 0.04])
    vuelta_textbox = TextBox(vuelta_ax_txt, "#laps[laps]:", num_vueltas, textalignment="center")

    M_AN_ax_txt = fig.add_axes([0.8, 0.70, 0.15, 0.04])
    M_AN_textbox = TextBox(M_AN_ax_txt, "Masa AN[UMA]:", M_AN, textalignment="center")

    x0_ax_txt = fig.add_axes([0.8, 0.65, 0.15, 0.04])
    x0_textbox = TextBox(x0_ax_txt, "X inicial[metros]:", x0, textalignment="center")

    # Llamadas a update
    # Sliders y checkbutton
    vx_slider.on_changed(update)
    vy_slider.on_changed(update)
    checkbox_e.on_clicked(update)
    # Textboxes
    delta_textbox.on_submit(update)
    t_max_textbox.on_submit(update)
    vuelta_textbox.on_submit(update)
    M_AN_textbox.on_submit(update)
    x0_textbox.on_submit(update)

    plt.tight_layout()

    plt.show()