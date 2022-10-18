#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 2022

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.wires
import matplotlib.pyplot as plt

# %% Function -plot-

def plot_cutted_contours2d(contour1, contour2, contours):
    n = int((len(contours) + 1) /2)
    if int((len(contours) + 1) /2) != ((len(contours) + 1) /2):
        n += 1

    count = 0
    fig, axs = plt.subplots(2, n)
    for i in range(0, 2):
        for j in range(0, n):
            if i == j and i == 0:
                contour1.plot(ax=axs[i][j])
                for prim in contour2.primitives:
                    prim.plot(ax=axs[i][j], width=2, color='r')
                axs[i][j].set_title("Contour2d 1 + Contour2d 2")

            else:
                contour1.plot(ax=axs[i][j])
                contour2.plot(ax=axs[i][j], color='r')
                for p in contours[count].primitives:
                    p.plot(ax=axs[i][j], width=2, color='b')
                    axs[i][j].set_title("Cutted Contour2d n° "+str(count+1))
                count += 1

# %% Contour 1

points = [
    vm.Point2D(0.20308817713481986, 0.04966773764193705),
    vm.Point2D(0.7969119765656952, 0.04966797879396336),
    vm.Point2D(0.8442697800221348, 0.04966805448106126),
    vm.Point2D(0.9031379981759905, 0.04545485716230839),
    vm.Point2D(0.9479619187253645, 0.10517762923653938),
    vm.Point2D(0.9545454713622077, 0.19659816573695035),
    vm.Point2D(0.9193352753725873, 0.27036785408838365),
    vm.Point2D(0.9193352733808331, 0.2891129034169921),
    vm.Point2D(0.9193352491395294, 0.7108874440320542),
    vm.Point2D(0.9193352445733572, 0.7296301344365137),
    vm.Point2D(0.9545453968718225, 0.8034021061594617),
    vm.Point2D(0.9479618025250582, 0.8948226395224322),
    vm.Point2D(0.9031378875027738, 0.9545453756361115),
    vm.Point2D(0.8442696587355882, 0.9503323074053969),
    vm.Point2D(0.7969118264212707, 0.950332267803071),
    vm.Point2D(0.203088054423346, 0.9503320443667395),
    vm.Point2D(0.15573023404235076, 0.9503319931770652),
    vm.Point2D(0.09686198756850553, 0.9545451189134176),
    vm.Point2D(0.05203811911069425, 0.8948223672618099),
    vm.Point2D(0.04545450905371345, 0.8034018054641981),
    vm.Point2D(0.08066468979315351, 0.7296322065226041),
    vm.Point2D(0.08066468965825305, 0.7108871704475791),
    vm.Point2D(0.08066474955300275, 0.289112539074918),
    vm.Point2D(0.08066475417207301, 0.27036750263798714),
    vm.Point2D(0.04545459207437027, 0.19659788603626316),
    vm.Point2D(0.05203818215846509, 0.10517733497317228),
    vm.Point2D(0.09686210943003962, 0.04545460320088811),
    vm.Point2D(0.15573034610061637, 0.04966777166777087)]

contour1 = volmdlr.wires.Contour2D.from_points(points)

# %% Contour 2

points = [
    vm.Point2D(0.2030881575366132, 0.04966771677601732),
    vm.Point2D(0.20308809125575447, 0.2891126765655333),
    vm.Point2D(0.08066474910267005, 0.2891125988482983),
    vm.Point2D(0.05674291581332371, 0.28911260559053886),
    vm.Point2D(0.04267127288273421, 0.311080039363985),
    vm.Point2D(0.04267123136672513, 0.6889196955746452),
    vm.Point2D(0.056742861492621116, 0.7108870816261534),
    vm.Point2D(0.08066469250960157, 0.7108871106086502),
    vm.Point2D(0.2030880435632797, 0.7108871156873701),
    vm.Point2D(0.20308804388447252, 0.9503320678576602),
    vm.Point2D(0.20308798646643267, 0.9876766506444887),
    vm.Point2D(0.21715969201684732, 1.0),
    vm.Point2D(0.7828401610275215, 1.0),
    vm.Point2D(0.7969118753894622, 0.9876767973117213),
    vm.Point2D(0.7969118376396574, 0.9503322697795179),
    vm.Point2D(0.796911890394745, 0.7108872963268008),
    vm.Point2D(0.9193352406490979, 0.7108873750750229),
    vm.Point2D(0.943257068363107, 0.7108873714262228),
    vm.Point2D(0.957328707411281, 0.6889199122213159),
    vm.Point2D(0.9573287795940568, 0.3110803289147694),
    vm.Point2D(0.9432571269950443, 0.28911295065634096),
    vm.Point2D(0.9193352813560651, 0.2891129208514467),
    vm.Point2D(0.7969119478569466, 0.289112868790721),
    vm.Point2D(0.7969119771240915, 0.049667935920999225),
    vm.Point2D(0.7969119981329132, 0.012323366834239619),
    vm.Point2D(0.7828402953854121, 0.0),
    vm.Point2D(0.21715982313378532, 0.0),
    vm.Point2D(0.2030881157414971, 0.012323183603131671)]


contour2 = volmdlr.wires.Contour2D.from_points(points)

# %% Cut_by_Wire

# %%% contour1.cut_by_wire(contour2)

results = contour1.cut_by_wire(contour2)
plot_cutted_contours2d(contour1, contour2, results)


# %%% contour2.cut_by_wire(contour1)

results = contour2.cut_by_wire(contour1)
plot_cutted_contours2d(contour2, contour1, results)
