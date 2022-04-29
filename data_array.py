# this file just stores the arrays with the corral/spectrum dat file and analysis information
# these arrays are used by the file scatter_params.py, which calls the functions in analyze_data.py


from collections import namedtuple

fields =  ['datfile', 'height_percentile','vertfile','marker1','marker2', "dataclipmin", "dataclipmax", "fit_order", "edge_cutoff", "chan", "type_fit"]
corralspectrum = namedtuple('corralspectrum', fields, defaults=(None,)*len(fields))
# fit parameters for "wtimes"
Co_Co_corrals = [
    corralspectrum("03-28/A220328.182736.dat", 96, "03-28/A220328.183230.L0015.VERT", -14, 40, -60, 60, 3, 0.3, "default"), # 2.9 nm
    corralspectrum("03-24/1/A220323.115455.dat", 97, "03-24/1/A220323.115740.VERT", -30, 14, -80, 80, 2, None, "default"),# 3.34 nm,
    corralspectrum("4-12 Co-Co/A220412.231331.dat", 99, "4-12 Co-Co/A220412.231512.VERT", 1, 10, None, None,3, None), #4.42 nm
    corralspectrum("03-26/A220325.102511.dat", 99, "03-26/A220325.102649.VERT", -5,20),# 4.5 nm
    corralspectrum("04-04/A220404.160833.dat", 99, "04-04/A220404.160935.VERT", -10, 15, -20, 45, 3), #4.478 nm
    corralspectrum("04-05 6nm Co/A220405.181958.dat", 99, "04-05 6nm Co/A220405.182713.VERT", -10, 20), # 6 nm
    corralspectrum("04-05 6nm Co/A220406.092143.dat", 99, "04-05 6nm Co/A220406.092807.L0010.VERT", -13, 20), # 6nm
    corralspectrum("4-12 Co-Co/A220412.115834.dat", 99, "4-12 Co-Co/A220412.120450.VERT", -10, 40, -60, 60, 3), # 5 nm radius
    corralspectrum("4-12 Co-Co/A220412.132510.dat", 99, "4-12 Co-Co/A220412.132822.L0017.VERT", -5, 15, -40, 40, 3, None, None, "wtimes"), # 3.8 nm
    corralspectrum("4-12 Co-Co/A220412.154921.dat", 99, "4-12 Co-Co/A220412.160413.VERT", 0, 50, 0, 50, 3), #4.7
    corralspectrum("4-12 Co-Co/A220412.163342.dat", 97, "4-12 Co-Co/A220412.163523.VERT", -30, 18, None, None, 3, None, None, "default"), #3 nm
    corralspectrum("4-12 Co-Co/A220412.183001.dat", 99, "4-12 Co-Co/A220412.184138.L0025.VERT", -10, 20, None, None, None, 2.5), # 6 nm
    corralspectrum("4-12 Co-Co/A220412.223356.dat", 99, "4-12 Co-Co/A220412.223556.VERT", -16, 20, -70, 60, 3, None, None), # 6.47
    corralspectrum("4-12 Co-Co/A220412.224607.dat", 99, "4-12 Co-Co/A220412.224743.VERT", -12, 35, -20, 45, 3, None), # 7 nm
    corralspectrum("4-12 Co-Co/A220412.225804.dat", 99, "4-12 Co-Co/A220412.230002.VERT", -10, 20), # 5.14 nm
    corralspectrum("4-12 Co-Co/A220412.233728.dat", 99, "4-12 Co-Co/A220412.233926.VERT", 0, 13, -40, 60, 3, None, None), # 4 nm
    corralspectrum("4-12 Co-Co/A220412.235220.dat", 98, "4-12 Co-Co/A220412.235427.VERT", -20, 20, -20, 20, 3, None, None), # 3.5 nm
    corralspectrum("4-12 Co-Co/A220413.103724.dat", 95, "4-12 Co-Co/A220413.103854.VERT", -20, 40, -50, 60, 3, None, None), # 2.9 nm
    corralspectrum("4-12 Co-Co/A220413.105651.dat", 96, "4-12 Co-Co/A220413.105932.VERT", -13, 31), # 2.65 nm
    corralspectrum("4-12 Co-Co/A220413.113527.dat", 99, "4-12 Co-Co/A220413.113654.VERT", -15, 15, -45, 60, 3), #6.4
    corralspectrum("03-29/A220329.113009.dat", 99, "03-29/A220329.113106.VERT", -15, 15, None, None, 3, None, None), # 7.42 nm
    corralspectrum("04-01/A220401.010342.dat", 99, "04-01/A220401.011026.L0125.VERT",0, 50, None, None, 3), # 8 nm
    corralspectrum("04-15 Co-Co\A220415.180247.dat", 99, "04-15 Co-Co\A220415.180404.VERT", 0, 12, None, None, None, None, 0), # 8.3 nm
    corralspectrum("04-14 Co Co\A220414.194555.dat", 99, "04-14 Co Co\A220414.195507.VERT", -2, 25, None, None, 3, None, 0), # 4.19 nm
    corralspectrum("04-14 Co Co\A220414.201501.dat", 98, "04-14 Co Co\A220414.201655.VERT", -2, 11, None, None, 3, None, 0), # 3.7 nm
    corralspectrum("04-14 Co Co\A220414.202552.dat", 97, "04-14 Co Co\A220414.202911.VERT", -2, 60,  -50, 60, 3, None, 0), # 3.38 nm
    corralspectrum("04-14 Co Co\A220414.204203.dat", 99, "04-14 Co Co\A220414.204346.VERT", -2, 12, -10, 70, None, None, 0, "wtimes"), # 3.9 nm
    corralspectrum("04-14 Co Co\A220414.205921.dat", 99, "04-14 Co Co\A220414.210247.VERT", 0, 13, -60, 60, 3, None, 0), # 4.4 nm
    corralspectrum("04-14 Co Co\A220414.212955.dat", 99, "04-14 Co Co\A220414.213310.VERT", -8, 25, None, None, None, None, 0) , # 5.35 nm
    corralspectrum("04-17 Co Co\position dependence experiment\A220417.213810.dat", 99, "04-17 Co Co\position dependence experiment\A220417.214221.VERT", -18, 40, None, None, 3, None, 0)  # 6.9 nm
]

Co_Ag_corrals = [
    corralspectrum("04-11 Ag Co\A220412.010237.dat", 97, "04-11 Ag Co\A220412.010418.VERT", -5, 40, -10, 60, None ),
    corralspectrum("04-11 Ag Co\A220411.141643.dat", 99, "04-11 Ag Co\A220411.141923.L0017.VERT", -8, 17, -30, 64, 3),
    corralspectrum("04-11 Ag Co\A220411.145437.dat", 99, "04-11 Ag Co\A220411.145852.VERT", -2, 15),
    corralspectrum("04-11 Ag Co\A220411.153007.dat", 99, "04-11 Ag Co\A220411.153513.VERT", -9, 15, -10, 20, 3),
    corralspectrum("04-11 Ag Co\A220411.153007.dat", 99, "04-11 Ag Co\A220411.154106.L0017.VERT", -8, 12, -10, 21, 3),
    corralspectrum("04-11 Ag Co\A220411.161126.dat", 99, "04-11 Ag Co\A220411.161806.L0017.VERT", 0, 11, -5, 20, 3),
    corralspectrum("04-11 Ag Co\A220411.165133.dat", 99, "04-11 Ag Co\A220411.165336.VERT", -13, 50, -25, 60, 3),
    corralspectrum("04-11 Ag Co\A220411.183528.dat", 99, "04-11 Ag Co\A220411.183838.VERT", -13, 30),
    corralspectrum("04-11 Ag Co\A220411.173719.dat", 99, "04-11 Ag Co\A220411.174011.VERT", -12, 15, -25, 50, 3),
    corralspectrum("04-11 Ag Co\A220411.193017.dat", 99, "04-11 Ag Co\A220411.193232.VERT", -7, 25),
    corralspectrum("04-11 Ag Co\A220411.200858.dat", 99, "04-11 Ag Co\A220411.201104.VERT", -7, 25),
    corralspectrum("04-11 Ag Co\A220411.204552.dat", 99, "04-11 Ag Co\A220411.204741.VERT", -8, 17),
    corralspectrum("04-11 Ag Co\A220411.214417.dat", 99, "04-11 Ag Co\A220411.214626.VERT", -4, 15, -10, 22, 3),
    corralspectrum("04-11 Ag Co\A220411.222442.dat", 99, "04-11 Ag Co\A220411.222845.L0016.VERT", -20, 15, -20, 22, 3),
    corralspectrum("04-11 Ag Co\A220411.233446.dat", 99, "04-11 Ag Co\A220411.233625.VERT", -15, 19, -20, 22, 3),
    corralspectrum("04-11 Ag Co\A220411.233446.dat", 99, "04-11 Ag Co\A220411.233844.L0015.VERT", -15, 14, -20, 22, 3),
    corralspectrum("04-06 6nm Ag walls\A220407.155505.dat", 99, "04-06 6nm Ag walls\A220407.160008.L0071.VERT", -5, 20)
]
