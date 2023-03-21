import seaborn as sns

#color = sns.color_palette("hls", 9).as_hex()

color = sns.color_palette("viridis", 4).as_hex()

#print(color)

hls           = ['#db5f57', '#dbb757', '#a7db57', '#57db5f', '#57dbb7', '#57a7db', '#5f57db', '#b757db', '#db57a7']
colorblind    = ['#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#ece133', '#56b4e9']
spectral      = ['#d43d4f', '#f46d43', '#fdad60', '#fee08b', '#ffffbe', '#e6f598', '#aadca4', '#66c2a5', '#3387bc']
seagreen_dark = ['#232724', '#24332a', '#254031', '#274c37', '#28593e', '#2a6644', '#2b724a', '#2d7f51', '#2e8b57']
rocket_r      = ['#4c1d4b', '#a11a5b', '#e83f3f', '#f69c73']
viridis       = ['#414487', '#2a788e', '#22a884', '#7ad151']

colors = [['#a11a5b', '#029e73', '#de8f05', '#0173b2'],
 		             ['#66c2a5', '#dbb757', '#0173b2'],
		  ['#87429b', '#66c2a5', '#d55e00', '#56b4e9']]

color_rad = colors[0][::-1]

print(color_rad)

['#0173b2', '#de8f05', '#029e73', '#a11a5b']

'''
1 1x1
2 2x1
3 2x2
4 3x1
5 3x1.5
6 3x2.5
7 3x3
8 3x5
9 4x1
'''