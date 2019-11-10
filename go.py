import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def periodic(xtemp, ytemp, xflux, yflux, px, py, Lx, Ly):
    if xtemp < 0:
        xtemp = xtemp + Lx
        xflux = xflux - px

    if xtemp > Lx:
        xtemp = xtemp - Lx
        xflux = xflux + px

    if ytemp < 0:
        ytemp = ytemp + Ly
        yflux = yflux - py

    if ytemp > Ly:
        ytemp = ytemp - Ly
        yflux = yflux + py

    return xtemp, ytemp, xflux, yflux


def periodic_hole(xtemp, ytemp, xflux, yflux, px, py, Lx, Ly, i, local_N, mask):
    if xtemp <= 0:
        if ytemp >= Ly / 3 and ytemp <= 2 * Ly / 3:
            local_N[i] = 0
            mask[i] = False
            xtemp = -100
        else:
            xtemp = xtemp + Lx
            xflux = xflux - px[i]

    if xtemp > Lx:
        xtemp = xtemp - Lx
        xflux = xflux + px[i]

    if ytemp < 0:
        ytemp = ytemp + Ly
        yflux = yflux - py[i]

    if ytemp > Ly:
        ytemp = ytemp - Ly
        yflux = yflux + py[i]

    return xtemp, ytemp, xflux, yflux


def mixed_hole(xtemp, ytemp, xflux, yflux, px, py, Lx, Ly, i, local_N, mask):
    if xtemp <= 0:
        if ytemp >= Ly / 3 and ytemp <= 2 * Ly / 3:
            local_N[i] = 0
            mask[i] = False
            xtemp = -100
        else:
            px[i] = - px[i]
            xflux = xflux - px[i]

    if xtemp >= Lx:
        px[i] = - px[i]
        xflux = xflux + px[i]

    if ytemp < 0:
        ytemp = ytemp + Ly
        yflux = yflux - py[i]

    if ytemp > Ly:
        ytemp = ytemp - Ly
        yflux = yflux + py[i]

    return xtemp, ytemp, xflux, yflux


def strong(xtemp, ytemp, xflux, yflux, px, py, Lx, Ly, i, local_N, mask):
    if xtemp <= 0:
        if ytemp >= Ly / 3 and ytemp <= 2 * Ly / 3:
            local_N[i] = 0
            mask[i] = False
            xtemp = -100
        else:
            # xtemp = - xtemp
            px[i] = - px[i]
            xflux = xflux - px[i]

    if xtemp >= Lx:
        px[i] = - px[i]
        # xtemp = Lx - (xtemp - Lx)
        xflux = xflux + px[i]

    if ytemp <= 0:
        py[i] = - py[i]
        # ytemp = - ytemp
        yflux = yflux - py[i]

    if ytemp >= Ly:
        py[i] = - py[i]
        # ytemp = Ly - (ytemp - Ly)
        yflux = yflux + py[i]


    return xtemp, ytemp, xflux, yflux


def separation(dx, dy, Lx, Ly):
    if abs(dx) > (0.5 * Lx):
        dx = dx * (1.0 - Lx / abs(dx))
    if abs(dy) > (0.5 * Ly):
        dy = dy * (1.0 - Ly / abs(dy))
    return dx, dy


def separation_per_hole(dx, dy, Lx, Ly, x1, y1, x2, y2):

    # if abs(dx) > (0.5 * Lx):
    #     dx = dx * (1.0 - Lx / abs(dx))
    if abs(dy) > (0.5 * Ly):
        dy = dy * (1.0 - Ly / abs(dy))
    return dx, dy


def sep(dx, dy, Lx, Ly):
    if abs(dx) > Lx:
        print(dx)
        dx = Lx
    if abs(dy) > Ly:
        print(dy)
        dy = Ly
    return dx, dy

def f(r):
    ri = 1.0 / r
    # print(r)
    ri3 = ri * ri * ri
    ri6 = ri3 * ri3
    g = 24.0 * ri * ri6 * (2.0 * ri6 - 1.0)
    force = g * ri
    potential = 4.0 * ri6 * (ri6 - 1.0)
    return force, potential


def accel(x, y, ax, ay, N, Lx, Ly, pe, list1, list2, mask):
    for i in range(N):
        ax[i] = 0
        ay[i] = 0
    for i in range(N - 1):
        if mask[i]:
            for j in range(i + 1, N):
                if mask[j]:
                    dx = x[i] - x[j]
                    dy = y[i] - y[j]
                    # dx, dy = separation(dx, dy, Lx, Ly)
                    # dx, dy = separation_per_hole(dx, dy, Lx, Ly, x[i], y[i], x[j], y[j])
                    # dx, dy = sep(dx, dy, Lx, Ly)
                    r = np.sqrt(dx ** 2 + dy ** 2)
                    force, potential = f(r)
                    if force > 0:
                        # print("Расс \t, Сила \t, i \t, j \t")
                        # print(r, "\t", force, "\t", i, "\t", j, "\t")
                        if list1.__len__() == 0:
                            list1.append(i)
                            list2.append(j)
                        else:
                            if j != list2[-1]:
                                list1.append(i)
                                list2.append(j)
                    ax[i] = ax[i] + force * dx
                    ay[i] = ay[i] + force * dy

                    ax[j] = ax[j] - force * dx
                    ay[j] = ay[j] - force * dy

                    pe = pe + potential
    return pe


def Verlet(x, y, vx, vy, ax, ay, N, Lx, Ly, dt, virial, xflux, yflux, pe, ke, list1, list2, local_N, mask):
    for i in range(N):
        if mask[i]:
            xnew = x[i] + vx[i] * dt + 0.5 * ax[i] * dt ** 2
            ynew = y[i] + vy[i] * dt + 0.5 * ay[i] * dt ** 2
            # Частично меняем скорость, используя старое ускорение
            vx[i] = vx[i] + 0.5 * ax[i] * dt
            vy[i] = vy[i] + 0.5 * ay[i] * dt
            # Периодические краевые условия и вычисление потока
            # xnew, ynew, xflux, yflux = periodic(xnew, ynew, xflux, yflux, vx[i], vy[i], Lx, Ly)
            # Не забудь убрать сепарейшен
            xnew, ynew, xflux, yflux = strong(xnew, ynew, xflux, yflux, vx, vy, Lx, Ly, i, local_N, mask)
            # xnew, ynew, xflux, yflux = periodic_hole(xnew, ynew, xflux, yflux, vx, vy, Lx, Ly, i, local_N, mask)
            # xnew, ynew, xflux, yflux = mixed_hole(xnew, ynew, xflux, yflux, vx, vy, Lx, Ly, i, local_N, mask)
            x[i] = xnew
            y[i] = ynew
    pe = accel(x, y, ax, ay, N, Lx, Ly, pe, list1, list2, mask)
    for i in range(N):
        # Окончательно меняем скорость, используя новое ускорение
        if mask[i]:
            vx[i] = vx[i] + 0.5 * ax[i] * dt
            vy[i] = vy[i] + 0.5 * ay[i] * dt
            ke = ke + 0.5 * (vx[i] ** 2 + vy[i] ** 2)
            virial = virial + x[i] * ax[i] + y[i] * ay[i]
    return ke, virial, xflux, yflux, pe


def results(N, nave, Lx, Ly, dt, virial, ke, pe, xflux, yflux, time):
    if time == 0:

        time = time + dt * nave
        ke = ke / nave
        pe = pe / nave
        # Энергия, приходящаяся на частицу
        E = (pe + ke) / N
        # Температура
        T = ke / N
        ke = 0
        pe = 0
        # Приведенное давление из вычисления потока
        pflux = ((xflux / (2.0 * Lx)) + (yflux / (2.0 * Ly))) / (dt * nave)
        xflux = 0
        yflux = 0
        # Приведенное давление из вычисления анриала
        pvirial = (N * T) / (Lx * Ly) + 0.5 * virial / (nave * Lx * Ly)
        virial = 0
        # if E > 1000:
        print("Время \t, T \t, E \t, pflux \t, pvirial \t")
        print(time, T, "\t", E, "\t", pflux, "\t", pvirial, "\t", )


if __name__ == "__main__":
    print("Число частиц")
    # N = eval(input("N "))
    # N = 400
    N = 100
    #Координаты
    x = np.zeros(N)
    y = np.zeros(N)
    mask = np.zeros(N, dtype=np.bool)
    mask[:] = True

    #Скорости
    vx = np.zeros(N)
    vy = np.zeros(N)

    #Ускорения
    ax = np.zeros(N)
    ay = np.zeros(N)

    print("Размер ящика Lx")
    # Lx = eval(input("Lx "))
    # Lx = 50
    Lx = 20
    print("Размер ящика Ly")
    # Ly = eval(input("Ly "))
    # Ly = 50
    Ly = 20

    print("Шаг по времени")
    # dt = eval(input("dt "))
    # dt = 0.001
    dt = 0.002
    dt2 = dt * dt

    print("Число шагов по времени между усреднениями")
    # nave = eval(input("nave "))
    nave = 400
    print("Количество наборов усреднения")
    # nset = eval(input("nset "))

    nset = 10
    x_data = np.zeros((nset * nave + 1, N))
    y_data = np.zeros((nset * nave + 1, N))
    t = np.arange(N)


    # начинаем на треугольной решетке
    print("Число частиц в ряду")
    # nrow = eval(input("nrow "))
    nrow = 10
    print("Максимальная скорость")
    # vmax = eval(input("vmax "))
    # vmax = 50
    vmax = 30
    #Сколько рядов
    ncol = N // nrow
    # расстояние до близжайшего соседа в каждом направлении
    disty = Ly / nrow
    distx = Lx / ncol
    i = 0
    vxsum = 0
    vysum = 0

    #data = np.zeros((N, N))

    #По всем рядам
    for icol in range(ncol):
        #Внутри ряда
        for irow in range(nrow):
            #data[icol, irow] = 1
            y[i] = disty * (irow + 1 - 0.5)
            if irow % 2 == 0:
                x[i] = distx * (icol + 1 - 0.25)
            else:
                x[i] = distx * (icol + 1 - 0.75)
            # проверь функцию рандом в паскале
            vx[i] = np.random.rand() * vmax
            vy[i] = np.random.rand() * vmax
            vxsum = vxsum + vx[i]
            vysum = vysum + vy[i]
            i = i + 1

    x_data[0] = x.copy()
    y_data[0] = y.copy()

    vxsum = vxsum / N
    vysum = vysum / N
    for i in range(N):
        vx[i] = vx[i] - vxsum
        vy[i] = vy[i] - vysum
    pe = 0
    list1 = []
    list2 = []
    list1.clear()
    list2.clear()
    pe = accel(x, y, ax, ay, N, Lx, Ly, pe, list1, list2, mask)
    time = 0
    pe = 0
    ke = 0
    xflux = 0
    yflux = 0
    virial = 0

    sect = 0
    triple = 0

    global_counter = 0
    local_N = np.ones(N)
    f1 = open("demofile2.txt", "w")

    sum_ln_y = 0
    sum_t = 0
    sum_t2 = 0
    sum_ln_y_t = 0
    local_counter = 0
    t1 = []
    y1 = []
    last = -1
    #Количество наборов усреднения
    for iset in range(nset):
        #Число шагов по времени между усреднениями
        for iave in range(nave):
            ke, verial, xflux, yflux, pe = Verlet(x, y, vx, vy, ax, ay, N, Lx, Ly, dt, virial, xflux, yflux, pe, ke, list1, list2, local_N, mask)
            dct = {}
            dct.clear()
            local_triple = 0
            if list1.__len__() > 2:
                for i in list1:
                    if i in dct:
                        dct[i] += 1
                    else:
                        dct[i] = 1


                for i in sorted(dct):
                    if dct[i] > 1:
                        # print(i)
                        # print(list1)
                        # print(list2)
                        # print("###")
                        find_i1 = list1.index(i)
                        find_i2 = find_i1 + 1
                        # if list2[find_i1] in list1 and list2[find_i2] in list2:
                        if list2[find_i1] in list1:
                            if list2[find_i2] == list2[list1.index(list2[find_i1])]:
                                triple = triple + 1
                                local_triple = local_triple + 1

            sect = sect + list1.__len__() - local_triple

            global_counter += 1
            list1.clear()
            list2.clear()

            x_data[global_counter] = x.copy()
            y_data[global_counter] = y.copy()
            if global_counter % 50 == 0 and local_N.sum() != last:
                local_counter = local_counter + 1
                sum_ln_y = sum_ln_y + np.log(local_N.sum())
                sum_t = sum_t + global_counter * dt
                sum_t2 = sum_t2 + (global_counter * dt) ** 2
                sum_ln_y_t = sum_ln_y_t + np.log(local_N.sum()) * global_counter * dt

                t1.append(global_counter * dt)
                y1.append(local_N.sum())
                last = local_N.sum()
                f1.write("%f\t" % (global_counter * dt))
                f1.write("%s\n" % local_N.sum())
                # print(global_counter, local_N.sum())
            results(N, nave, Lx, Ly, dt, virial, ke, pe, xflux, yflux, time)

    print(sect)
    print(triple)
    print(mask.sum())
    log_A = (sum_t * sum_ln_y_t - sum_ln_y * sum_t2) / (sum_t ** 2 - t1.__len__() * sum_t2)
    B = (sum_ln_y - log_A * t1.__len__()) / (sum_t)
    f1.close()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)

    scatter = ax.scatter(x_data[0], y_data[0], c=t, zorder=1, s=25)  # Scatter plot
    time_text = ax.text(0.05, 0.95, '', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)


    def animate(i):
        if i % 2 == 0:
            scatter.set_offsets(np.c_[x_data[i, :], y_data[i, :]])
            # plt.vlines(0.1, ymin=Ly / 3, ymax=2 * Ly / 3)
        # xmin = x_data[i,:].min(); xmax = x_data[i,:].max()
        # ymin = y_data[i,:].min(); ymax = y_data[i,:].max()
        # ax.set_xlim(xmin - 0.1 * (xmax - xmin), xmax + 0.1 * (xmax - xmin))
        # ax.set_ylim(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin))
            time_text.set_text('time = %.1d' % i)



ani = animation.FuncAnimation(fig, animate, frames=len(x_data / 2),
                                  interval=1, blit=False)
plt.show()
#
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
my_func = []
for i in t1:
    my_func.append(np.exp(log_A) * np.exp(B * i))
plt.title(str(np.exp(log_A)) + "*e^" + str(B) + "t")
ax1.plot(t1, my_func)
plt.scatter(t1, y1, label= "stars", color= "green",
            marker= "*", s=30)
# plt.plot(t1[:], y1[:])
plt.show()
