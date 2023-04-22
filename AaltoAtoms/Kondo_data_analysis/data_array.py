# this file stores arrays with corral/spectrum dat file and analysis information
# these are used by scatter_params.py, which calls functions in analyze_data.py

# quick-and-dirty way to make an object
from collections import namedtuple

fields =  ['datfile', 'height_percentile','vertfile','marker1','marker2', "dataclipmin", "dataclipmax", "fit_order", "edge_cutoff", "chan", "type_fit"]
corralspectrum = namedtuple('corralspectrum', fields, defaults=(None,)*len(fields))
Co_Co_corrals = [
    corralspectrum("4-12 Co-Co/A220412.131659.dat", 99, "4-12 Co-Co/A220412.131810.VERT", -7, 30, -40, 40, None ), # 3.89 nm
    corralspectrum("03-28/A220328.182736.dat", 96, "03-28/A220328.183230.L0015.VERT", -20, 15, -30, 30, None, 0.3, "default"), # 2.9 nm
    corralspectrum("03-24/1/A220323.115455.dat", 97, "03-24/1/A220323.115740.VERT", -30, 13, -80, 40, 2, None, "default"),# 3.34 nm,
    corralspectrum("4-12 Co-Co/A220412.231331.dat", 99, "4-12 Co-Co/A220412.231512.VERT", 1, 10, -30, 30, None, None), #4.42 nm
    corralspectrum("03-26/A220325.102511.dat", 99, "03-26/A220325.102649.VERT", -5,20),# 4.5 nm
    corralspectrum("03-26/A220325.102511.dat", 99, "03-26/A220325.131310.VERT", 0,18, -20, 40, 3),# 4.51 nm
    corralspectrum("03-26/A220325.132236.dat", 99, "03-26/A220325.132407.VERT",  -5,14, -20, 40, None),# 4.5 nm
    corralspectrum("04-04/A220404.160833.dat", 99, "04-04/A220404.160935.VERT", -10, 15, -20, 45, 3), #4.478 nm
    corralspectrum("04-05 6nm Co/A220405.181958.dat", 99, "04-05 6nm Co/A220405.182713.VERT", -10, 20), # 6 nm
    corralspectrum("04-05 6nm Co/A220406.092143.dat", 99, "04-05 6nm Co/A220406.092807.L0010.VERT", -13, 20), # 6nm
    corralspectrum("4-12 Co-Co/A220412.115834.dat", 99, "4-12 Co-Co/A220412.120450.VERT", -10, 40, -60, 60, 3), # 5 nm radius
    corralspectrum("4-12 Co-Co/A220412.132510.dat", 99, "4-12 Co-Co/A220412.132822.L0017.VERT", -5, 15, -40, 40, None, None, None, "default"), # 3.8 nm
    corralspectrum("4-12 Co-Co/A220412.160314.dat", 99, "4-12 Co-Co/A220412.160413.VERT", 0, 50, 0, 50, 3), # 4.7 nm
    corralspectrum("4-12 Co-Co/A220412.163342.dat", 97, "4-12 Co-Co/A220412.163523.VERT", -30, 25, None, None, None, None, None, "default"), #3 nm
    corralspectrum("4-12 Co-Co/A220412.183001.dat", 99, "4-12 Co-Co/A220412.184138.L0025.VERT", -10, 20, None, None, None, 2.5), # 6 nm
    corralspectrum("4-12 Co-Co/A220412.223356.dat", 99, "4-12 Co-Co/A220412.223556.VERT", -16, 20, -30, 30, None, None, None), # 6.47
    corralspectrum("4-12 Co-Co/A220412.224607.dat", 99, "4-12 Co-Co/A220412.224743.VERT", -10, 15, -40, 45, None, None), # 7 nm
    corralspectrum("4-12 Co-Co/A220412.225804.dat", 99, "4-12 Co-Co/A220412.230002.VERT", -10, 20), # 5.14 nm
    corralspectrum("4-12 Co-Co/A220412.233728.dat", 99, "4-12 Co-Co/A220412.233926.VERT", -5, 20, -40, 40, None, None, None), # 4 nm
    corralspectrum("4-12 Co-Co/A220412.235220.dat", 98, "4-12 Co-Co/A220412.235427.VERT", -4, 20, -20, 20, None, None, None), # 3.5 nm
    corralspectrum("4-12 Co-Co/A220413.103724.dat", 95, "4-12 Co-Co/A220413.103854.VERT", -20, 20, -50, 60, None, None, None), # 2.9 nm
    corralspectrum("4-12 Co-Co/A220413.105651.dat", 96, "4-12 Co-Co/A220413.105932.VERT", -13, 31), # 2.65 nm
    corralspectrum("4-12 Co-Co/A220413.113527.dat", 99, "4-12 Co-Co/A220413.113654.VERT", -15, 15, -45, 60, 3), # 6.4 nm
    corralspectrum("03-29/A220329.113009.dat", 99, "03-29/A220329.113106.VERT", -20, 18, -40, 30, None, None, None, "default"), # 7.51 nm
    corralspectrum("04-01/A220401.010342.dat", 99, "04-01/A220401.011026.L0125.VERT",-2, 40, -5, 40, 3), # 8 nm
    corralspectrum("04-15 Co-Co/A220415.180247.dat", 99, "04-15 Co-Co/A220415.180404.VERT", 0, 12, None, None, None, None, 0), # 8.3 nm
    corralspectrum("04-14 Co Co/A220414.194555.dat", 99, "04-14 Co Co/A220414.195507.VERT", -2, 25, None, None, None, None, 0), # 4.19 nm
    corralspectrum("04-14 Co Co/A220414.201501.dat", 98, "04-14 Co Co/A220414.201655.VERT", -2, 11, None, None, 3, None, 0), # 3.7 nm
    corralspectrum("04-14 Co Co/A220414.202552.dat", 97, "04-14 Co Co/A220414.202911.VERT", -2, 60,  -20, 30, None, None, 0), # 3.38 nm
    corralspectrum("04-14 Co Co/A220414.204203.dat", 99, "04-14 Co Co/A220414.204346.VERT", -2, 30, -10, 70, None, None, 0, "default"), # 3.9 nm
    corralspectrum("04-14 Co Co/A220414.205921.dat", 99, "04-14 Co Co/A220414.210247.VERT", 0, 13, -60, 60, 3, None, 0), # 4.4 nm
    corralspectrum("04-14 Co Co/A220414.212955.dat", 99, "04-14 Co Co/A220414.213310.VERT", -8, 25, None, None, None, None, 0) , # 5.35 nm
    corralspectrum("04-17 Co Co/position dependence experiment/A220417.213810.dat", 99, "04-17 Co Co/position dependence experiment/A220417.214221.VERT", -25, 18, -50, 30, 3, None, 0)  # 6.9 nm
]

Co_Ag_corrals = [
    corralspectrum("04-11 Ag Co/A220411.215004.dat", 99, "04-11 Ag Co/A220411.215940.L0015.VERT", -5, 10, -30, 20, None, None, None, "default"), # 4.6 nm
    corralspectrum("04-11 Ag Co/A220412.010237.dat", 97, "04-11 Ag Co/A220412.010418.VERT", -2, 40, -10, 60, None ), # 3.63 nm
    corralspectrum("04-11 Ag Co/A220411.133438.dat", 99, "04-11 Ag Co/A220411.134351.L0013.VERT", -15, 25, -15, 25, 3, None, None), # 11.9 nm
    corralspectrum("04-11 Ag Co/A220411.141643.dat", 99, "04-11 Ag Co/A220411.141923.L0017.VERT", -7, 15, -80, 80, 3), # 10.07 nm
    corralspectrum("04-11 Ag Co/A220411.145437.dat", 99, "04-11 Ag Co/A220411.145852.VERT", -2, 15), # 9.16 nm
    corralspectrum("04-11 Ag Co/A220411.153007.dat", 99, "04-11 Ag Co/A220411.153513.VERT", -2, 12, -15, 28, 3), # 8.16 nm
    corralspectrum("04-11 Ag Co/A220411.153007.dat", 99, "04-11 Ag Co/A220411.154106.L0017.VERT", -8, 12, -20, 20, 3), # 8.1 nm
    corralspectrum("04-11 Ag Co/A220411.161126.dat", 99, "04-11 Ag Co/A220411.161806.L0017.VERT", -2, 12, None, None, None), # 7.6 nm
    # corralspectrum("04-11 Ag Co/A220411.165133.dat", 99, "04-11 Ag Co/A220411.165336.VERT", -20, 15, -40, 60, None, None, None, "default"), # 7.11 nm
    corralspectrum("04-11 Ag Co/A220411.165133.dat", 99, "04-11 Ag Co/A220411.165336.VERT", -0, 15, -10, 30, None, None, None, "default"), # 7.11 nm
    corralspectrum("04-11 Ag Co/A220411.183528.dat", 99, "04-11 Ag Co/A220411.183838.VERT", -13, 30), # 6.09 nm
    corralspectrum("04-11 Ag Co/A220411.173719.dat", 99, "04-11 Ag Co/A220411.174011.VERT", -12, 15, -25, 50, 3), # 6.82 nm
    corralspectrum("04-11 Ag Co/A220411.193017.dat", 99, "04-11 Ag Co/A220411.193232.VERT", -7, 25), # 5.8 nm
    corralspectrum("04-11 Ag Co/A220411.200858.dat", 99, "04-11 Ag Co/A220411.201104.VERT", -7, 25), # 5.55 nm
    corralspectrum("04-11 Ag Co/A220411.204552.dat", 99, "04-11 Ag Co/A220411.204741.VERT", -8, 17), # 5.14 nm
    corralspectrum("04-11 Ag Co/A220411.214417.dat", 99, "04-11 Ag Co/A220411.214626.VERT", -4, 15, -10, 22, 3), # 4.61 nm
    corralspectrum("04-11 Ag Co/A220411.222442.dat", 99, "04-11 Ag Co/A220411.222845.L0016.VERT", -4, 15, -80, 50, 3), # 4.08 nm
    corralspectrum("04-11 Ag Co/A220411.233446.dat", 99, "04-11 Ag Co/A220411.233625.VERT", -18, 19, -20, 22, 3), # 5.04 nm
    #corralspectrum("04-11 Ag Co/A220411.233446.dat", 99, "04-11 Ag Co/A220411.233844.L0015.VERT", -10, 20, -40, 40, None), # 4.19 nm

    corralspectrum("04-11 Ag Co/A220411.233446.dat", 99, "04-11 Ag Co/A220411.233844.L0015.VERT", -12, 30, -30, 30, 3), # 4.19 nm
    corralspectrum("04-06 6nm Ag walls/A220407.155505.dat", 99, "04-06 6nm Ag walls/A220407.160008.L0071.VERT", -5, 20), # 5.95 nm
    corralspectrum(r"Ag 2021-08-13 2p5 nm radius/2p5nm radius pm 20mV line spectrum/Createc2_210813.161840.dat", 97, r"Ag 2021-08-13 2p5 nm radius/2p5nm radius pm 20mV line spectrum/Createc2_210813.165557.L0030.VERT", -15, 20), # 2.44 nm
    corralspectrum(r"Ag 2021-08-13 2p5 nm radius/pm 100 mV 2p5 nm radius line spectrum/Createc2_210813.172359.dat", 97, r"Ag 2021-08-13 2p5 nm radius/pm 100 mV 2p5 nm radius line spectrum/Createc2_210813.173235.L0030.VERT", -20, 20), # 2.46 nm
    corralspectrum(r"Ag 2021-08-13 2p5 nm radius/300mV to -200mV line/Createc2_210813.231030.dat", 97, r"Ag 2021-08-13 2p5 nm radius/300mV to -200mV line/Createc2_210813.231403.L0030.VERT", -20, 20), # 2.46 nm
    corralspectrum(r"Ag 2021-08-13 3p8 nm radius/Createc2_210813.102220.dat", 99, r"Ag 2021-08-13 3p8 nm radius/3p8nm pm20mV line/Createc2_210813.105240.L0030.VERT", 0, 18), # 3.65 nm
    corralspectrum(r"Ag 2021-08-12 4p5 nm radius/4p5 nm line spectrum pm100mV/Createc2_210812.152302.dat", 99,  r"Ag 2021-08-12 4p5 nm radius/4p5 nm line spectrum pm100mV/Createc2_210812.154131.L0036.VERT", -15, 15),
    corralspectrum(r"Ag 2021-08-12 4p5 nm radius/pm 20mV line/Createc2_210812.170231.dat", 99, r"Ag 2021-08-12 4p5 nm radius/pm 20mV line/Createc2_210812.171526.L0036.VERT", -15, 15), #4.65 nm
    corralspectrum(r"Ag 2021-08-11/3p8 nm radius line spectra pm100mV/Createc2_210811.112143.dat",99, r"Ag 2021-08-11/3p8 nm radius line spectra pm100mV/Createc2_210811.113827.L0030.VERT", -2, 25), # 4.65 nm
    corralspectrum(r"Ag 2021-08-10 2p5 nm radius/2p5 nm radius pm20mV line spectrum/Createc2_210810.090437.dat",98, r"Ag 2021-08-10 2p5 nm radius/2p5 nm radius pm20mV line spectrum/Createc2_210810.091431.L0075.VERT", -20, 20) # 3.81 nm
]

H2Pc_corrals_center = [
    corralspectrum(r"2022-05-25 second deposition after noise issue/9p5 nm corral/A220527.000351.dat",99,  r"2022-05-25 second deposition after noise issue/9p5 nm corral/A220527.001034.L0050.VERT" , -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-25 second deposition after noise issue/9 nm corral/A220527.131257.dat", 99, r"2022-05-25 second deposition after noise issue/9 nm corral/A220527.132244.L0025.VERT", -20, 20, -20, 20, None, None, None, None ),
    corralspectrum(r"2022-05-25 second deposition after noise issue/8p57 nm corral/A220527.173842.dat",99, r"2022-05-25 second deposition after noise issue/8p57 nm corral/A220527.174243.L0024.VERT", -20, 20, -20, 20, None, None, None, None ),
    corralspectrum(r"2022-05-25 second deposition after noise issue/8p08 nm corral/A220527.185403.dat",99,  r"2022-05-25 second deposition after noise issue/8p08 nm corral/A220527.185837.L0024.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-25 second deposition after noise issue/7p5 nm corral/A220527.212711.dat",99,  r"2022-05-25 second deposition after noise issue/7p5 nm corral/A220527.214240.L0025.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/8 nm/A220524.132718.dat", 95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/8 nm/A220524.133236.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/7.3 nm corral/A220523.182042.dat", 94,r"2022-05-16 H2Pc (2 min 300 C -130 C)/7.3 nm corral/A220523.182756.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/7 nm corral/A220523.234720.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/7 nm corral/A220523.235603.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6p07 nm Monday corral/A220523.145912.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6p07 nm Monday corral/A220523.150635.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.7 nm corral/A220524.015654.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.7 nm corral/A220524.020013.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.3nm corral/A220524.094612.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.3nm corral/A220524.095007.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.1nm radius corral/A220524.105814.dat",99, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.1nm radius corral/A220524.110432.L0020.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/5p5 nm corral/A220522.235140.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/5p5 nm corral/A220522.235730.L0030.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/5nm corral/A220520.173130.dat",99, r"2022-05-16 H2Pc (2 min 300 C -130 C)/5nm corral/A220520.175130.L0040.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/4p5 nm corral/A220520.142331.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/4p5 nm corral/A220520.142629.L0030.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/4 nm corral/A220518.151648.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/4 nm corral/A220518.153129.L0030.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/3p5nm corral/A220518.181143.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/3p5nm corral/A220518.182241.L0025.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/3nm radius/A220519.211852.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/3nm radius/A220519.213027.L0016.VERT", -20, 20, -20, 20, None, None, None, None),
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/3.93 nm radius corral/A220520.094659.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/3.93 nm radius corral/A220520.095212.L0015.VERT", -20, 20, -20, 20, None, None, None, None),
]

H2Pc_corrals_edge = [
    corralspectrum(r"2022-05-25 second deposition after noise issue/9p5 nm corral/A220527.000351.dat", 99,  r"2022-05-25 second deposition after noise issue/9p5 nm corral/A220527.001034.L0045.VERT" , -20, 20, -20, 20, None, None, None, None), #9.48 nm
    corralspectrum(r"2022-05-25 second deposition after noise issue/9 nm corral/A220527.131257.dat", 99, r"2022-05-25 second deposition after noise issue/9 nm corral/A220527.132244.L0019.VERT", -20, 20, -20, 20, None, None, None, None ), #9.03 nm
    corralspectrum(r"2022-05-25 second deposition after noise issue/8p57 nm corral/A220527.173842.dat", 99, r"2022-05-25 second deposition after noise issue/8p57 nm corral/A220527.174243.L0020.VERT", -20, 20, -20, 20, None, None, None, None ), #8.55 nm
    corralspectrum(r"2022-05-25 second deposition after noise issue/8p08 nm corral/A220527.185403.dat", 99,  r"2022-05-25 second deposition after noise issue/8p08 nm corral/A220527.185837.L0019.VERT", -20, 20, -20, 20, None, None, None, None), #8.05 nm
    corralspectrum(r"2022-05-25 second deposition after noise issue/7p5 nm corral/A220527.212711.dat", 99,  r"2022-05-25 second deposition after noise issue/7p5 nm corral/A220527.214240.L0019.VERT", -20, 20, -20, 20, None, None, None, None), #7.5 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/8 nm/A220524.132718.dat", 95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/8 nm/A220524.133236.L0023.VERT", -20, 20, -20, 20, None, None, None, None), #8 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/7.3 nm corral/A220523.182042.dat", 94,r"2022-05-16 H2Pc (2 min 300 C -130 C)/7.3 nm corral/A220523.182756.L0015.VERT", -20, 20, -20, 20, None, None, None, None), #7.4 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/7 nm corral/A220523.234720.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/7 nm corral/A220523.235603.L0017.VERT", -20, 20, -20, 20, None, None, None, None), # 6.96 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6p07 nm Monday corral/A220523.145912.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6p07 nm Monday corral/A220523.150635.L0016.VERT", -20, 20, -20, 20, None, None, None, None), #6.05 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.7 nm corral/A220524.015654.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.7 nm corral/A220524.084536.L0014.VERT", -20, 20, -20, 20, None, None, None, None), # 6.67 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.3nm corral/A220524.094612.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.3nm corral/A220524.095007.L0016.VERT", -20, 20, -20, 20, None, None, None, None), #6.28 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.1nm radius corral/A220524.105814.dat",99, r"2022-05-16 H2Pc (2 min 300 C -130 C)/6.1nm radius corral/A220524.110432.L0023.VERT", -20, 20, -20, 20, None, None, None, None), # 6.10 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/5p5 nm corral/A220522.235140.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/5p5 nm corral/A220522.235730.L0026.VERT", -20, 20, -20, 20, None, None, None, None), #5.55 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/5nm corral/A220520.173130.dat",99, r"2022-05-16 H2Pc (2 min 300 C -130 C)/5nm corral/A220520.175130.L0047.VERT", -20, 20, -20, 20, None, None, None, None), # 5.10 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/4p5 nm corral/A220520.142331.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/4p5 nm corral/A220520.142629.L0026.VERT", -20, 20, -20, 20, None, None, None, None), # 4.53 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/4 nm corral/A220518.151648.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/4 nm corral/A220518.153129.L0026.VERT", -20, 20, -20, 20, None, None, None, None), #4.10 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/3p5nm corral/A220518.181143.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/3p5nm corral/A220518.182241.L0020.VERT", -20, 20, -20, 20, None, None, None, None), # 3.56 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/3nm radius/A220519.211852.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/3nm radius/A220519.213027.L0020.VERT", -20, 20, -20, 20, None, None, None, None), #3.62 nm
    corralspectrum(r"2022-05-16 H2Pc (2 min 300 C -130 C)/3.93 nm radius corral/A220520.094659.dat",95, r"2022-05-16 H2Pc (2 min 300 C -130 C)/3.93 nm radius corral\A220520.095212.L0009.VERT", -20, 20, -20, 20, None, None, None, None), #3.93 nm
]
