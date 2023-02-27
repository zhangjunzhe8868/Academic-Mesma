# select the most useful combo based on the spectral lib file and initial 1000 combo

import numpy as np
data1 = np.loadtxt('result\best_combo\models_subset_1000_combined_library_isric_aster' +
                   '_thoralf_814x2151.txt', dtype='float')
data2 = np.loadtxt('result\best_combo\mesma_stats.txt', dtype='float')
temp = np.zeros((len(data2), 3))

for i in range(len(data2)):
    print(data1[int(data2[i, 0]), :])
    temp[i, :] = (data1[int(data2[i, 0]), :])

np.savetxt('result\models_subset_1000_combined_library_isric_aster' +
                   '_thoralf_814x2151_sub_'+'CA.txt', temp, delimiter='    ',
           fmt='%d')
