from os import path, makedirs, listdir, getcwd
import numpy as np
from crystal_math import GTitoGT1_cubic_corres
from crystal_math import Qu
from crystal_math import OIMExport, Polefig
from crystal_math import cubicsymmetry
from crystal_math import Euler as Eu
from crystal_math import SpecialTrans as Sp
from crystal_math import AngleCalc as Ac
import math
import time
import json
import matplotlib.pyplot as plt

def plot_variant_selection(directory, schmid_factors_dic, quaveaus):
    slip_system_2_NW = {
        '1'   : '3',
        '-1'  : '2',
        '2'   : '1',
        '-2'  : '3',
        '3'   : '2',
        '-3'  : '1',
        '4'   : '5',
        '-4'  : '6',
        '5'   : '6',
        '-5'  : '4',
        '6'   : '4',
        '-6'  : '5',
        '7'   : '8',
        '-7'  : '9',
        '8'   : '9',
        '-8'  : '7',
        '9'   : '7',
        '-9'  : '8',
        '10'  : '12',
        '-10' : '11',
        '11'  : '10',
        '-11' : '12',
        '12'  : '11',
        '-12' : '10'
    }

    sf = [abs(float(schmid_factors_dic[str(x+1)])) for x in range(12)]

    sorted_sf = sorted(sf)[::-1]
    
    selected_variant = {}

    for y in range(12):
        for x in range(12):
            if sorted_sf[y] == abs(float(schmid_factors_dic[str(x+1)])):
                if float(schmid_factors_dic[str(x+1)]) >= 0:
                    if str(y+1) == '1':
                        selected_variant['Max N_W_variant'] = slip_system_2_NW[str(x+1)]
                    elif str(y+1) == '2':
                        selected_variant['2nd N_W_variant'] = slip_system_2_NW[str(x+1)]
                    elif str(y+1) == '3':
                        selected_variant['3rd N_W_variant'] = slip_system_2_NW[str(x+1)]
                    else:
                        selected_variant[str(y+1) + 'th ' + 'N_W_variant'] = slip_system_2_NW[str(x+1)]
                
                elif float(schmid_factors_dic[str(x+1)]) < 0:
                    if str(y+1) == '1':
                        selected_variant['Max N_W_variant'] = slip_system_2_NW[str(-(x+1))]
                    elif str(y+1) == '2':
                        selected_variant['2nd N_W_variant'] = slip_system_2_NW[str(-(x+1))]
                    elif str(y+1) == '3':
                        selected_variant['3rd N_W_variant'] = slip_system_2_NW[str(-(x+1))]
                    else:
                        selected_variant[str(y+1) + 'th ' + 'N_W_variant'] = slip_system_2_NW[str(-(x+1))]
    
    print('Deformation induced variant selection:')
    
    for k, v in selected_variant.items():
        print(k,'=',v)
    
    selected_variant = list(selected_variant.values())
    selected_variant = [int(selected_variant[x]) for x in range(12)]
    
    OIMExport.nw_variant_selecion_or_plot(quaveaus, directory, selected_variant)


def plot_schmid_factors(direc, filename, schmid_factors_dic):
    sf = [schmid_factors_dic[str(x+1)] for x in range(12)]
    n_groups = 12
    index = np.arange(n_groups)
    bar_width = 0.6
    opacity = 0.4
    

    plt.title('Schmid Factors on different slip systems')
    plt.xlabel('slip system')
    plt.ylabel('Schmid Factor')
    plt.xticks(index, ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'))
    plt.bar(left=index, height=tuple(sf), width=bar_width, align="center", color="r", label=str(filename), alpha=opacity)
    plt.axhline(color='k', linewidth=0.5)
    plt.legend()
    
    plt.savefig(direc + filename[:-4] + '_schmid_factor.png')
    

def know_schmid_factors(direc, filename, quaveaus):
 
    slip_systems_planes = {1 : [1.0, 1.0, 1.0],\
		                   2 : [1.0, 1.0, 1.0],\
		                   3 : [1.0, 1.0, 1.0],\
		                   4 : [1.0, 1.0, 1.0],\
		                   5 : [1.0, 1.0, 1.0],\
		                   6 : [1.0, 1.0, 1.0],\
		                   7 : [-1.0, 1.0, 1.0],\
		                   8 : [-1.0, 1.0, 1.0],\
		                   9 : [-1.0, 1.0, 1.0],\
		                   10 : [-1.0, 1.0, 1.0],\
		                   11 : [-1.0, 1.0, 1.0],\
		                   12 : [-1.0, 1.0, 1.0],\
		                   13 : [1.0, -1.0, 1.0],\
		                   14 : [1.0, -1.0, 1.0],\
		                   15 : [1.0, -1.0, 1.0],\
		                   16 : [1.0, -1.0, 1.0],\
		                   17 : [1.0, -1.0, 1.0],\
		                   18 : [1.0, -1.0, 1.0],\
		                   19 : [-1.0, -1.0, 1.0],\
		                   20 : [-1.0, -1.0, 1.0],\
		                   21 : [-1.0, -1.0, 1.0],\
		                   22 : [-1.0, -1.0, 1.0],\
		                   23 : [-1.0, -1.0, 1.0],\
		                   24 : [-1.0, -1.0, 1.0]}

    slip_systems_directions = {1 : [0.0, -1.0, 1.0],\
		                   2 : [0.0, 1.0, -1.0],\
		                   3 : [1.0, 0.0, -1.0],\
		                   4 : [-1.0, 0.0, 1.0],\
		                   5 : [-1.0, 1.0, 0.0],\
		                   6 : [1.0, -1.0, 0.0],\
		                   7 : [0.0, 1.0, -1.0],\
		                   8 : [0.0, -1.0, 1.0],\
		                   9 : [1.0, 0.0, 1.0],\
		                   10 : [-1.0, 0.0, -1.0],\
		                   11 : [-1.0, -1.0, 0.0],\
		                   12 : [1.0, 1.0, 0.0],\
		                   13 : [0.0, -1.0, -1.0],\
		                   14 : [0.0, 1.0, 1.0],\
		                   15 : [-1.0, 0.0, 1.0],\
		                   16 : [1.0, 0.0, -1.0],\
		                   17 : [1.0, 1.0, 0.0],\
		                   18 : [-1.0, -1.0, 0.0],\
		                   19 : [0.0, 1.0, 1.0],\
		                   20 : [0.0, -1.0, -1.0],\
		                   21 : [-1.0, 0.0, -1.0],\
		                   22 : [1.0, 0.0, 1.0],\
		                   23 : [1.0, -1.0, 0.0],\
		                   24 : [-1.0, 1.0, 0.0]}

    def_vector_on_OIM = [1.0, 0.0, 0.0]
    quaveausi = Qu.qu_inv(quaveaus)
    def_vector_center = Qu.qu_vec_rotation(quaveausi, def_vector_on_OIM)
    print(def_vector_center)
    
   
    d = def_vector_center
    ssp = slip_systems_planes
    ssd = slip_systems_directions

    schmid_factors_list = []
    schmid_factors_dic = {}

    for x in range(24):
        schmid_factors = {}
        cos_theta1 = Ac.cosVector(d, ssp[x+1])
        cos_theta2 = Ac.cosVector(d, ssd[x+1])

        schmid_factor = cos_theta1 * cos_theta2
        schmid_factors['{0}'.format(str(x+1))] = str(schmid_factor)
        schmid_factors_list.append(schmid_factors)
    
    for schmid_factor in schmid_factors_list:
        print(schmid_factor)

    with open(direc + filename[:-4] + '_schmid_factor.txt', 'w', encoding='utf-8') as filew:
        json.dump(schmid_factors_list, filew, indent=2, sort_keys=True, ensure_ascii=False)
    
    for x in range(12):
        cos_theta1 = Ac.cosVector(d, ssp[2*x+1])
        cos_theta2 = Ac.cosVector(d, ssd[2*x+1])

        schmid_factor = cos_theta1 * cos_theta2
        schmid_factors_dic['{0}'.format(str(x+1))] = str(schmid_factor)
    
    print('Schmid Factors on different slip systems:')
    print(schmid_factors_dic)

    return schmid_factors_list, schmid_factors_dic

def output_variant_distribution(direc, filename, variants_counts):
    filevar = open (direc + filename[:-4] + '_var_distribution.txt', 'w')
    
    variant_sum = np.sum(variants_counts)
    variants_counts = Sp.variant24_convert(variants_counts)
    

    for index, var_num in enumerate(variants_counts):
        strprint = '{0:2d}\t {1:6.4f} \n'\
                        .format(index + 1, var_num / variant_sum)
        filevar.write(strprint)

def plot_OIM_OR(direc, quave, quavefer, quavefer_original, quave_angle_std, quaveaus_angle_std, quaveaus_angle110_std, quaveaus_angle111_std):
    
    filew = OIMExport.oim_or_plot (quavefer_original, quave, direc)
    
    # 計算平均方位等之差異
    var_num, misang = Qu.qu_fervar_num_from_gt(quavefer)
    aus111angle = Qu.aus111fer110_deviation(quavefer, var_num)
    aus110angle = Qu.aus110fer111_deviation(quavefer, var_num)

	# 計算標準差，注意degree of freedom = 1，代表N-1
    quavestd = np.std(quave_angle_std, ddof=1)
    quaveferstd = np.std(quaveaus_angle_std, ddof=1)
    quaveferstd110 = np.std(quaveaus_angle110_std, ddof=1)
    quaveferstd111 = np.std(quaveaus_angle111_std, ddof=1)

    filew.write('\n\nOR')
    filew.write('{0}'.format(quavefer_original))

    filew.write('\ndeviation of austenite mean grain: {0:4.2f}'.format(quavestd))
    filew.write('\ndeviation of ferrite mean grain:   {0:4.2f}'.format(quaveferstd))

    filew.write('\nDeviated From Aus111: {0:4.4f}'.format(aus111angle))
    filew.write('\t deviation:  {0:4.2f}'.format(quaveferstd111))

    filew.write('\nDeviated From Aus110: {0:4.4f}'.format(aus110angle))
    filew.write('\t deviation:  {0:4.2f}'.format(quaveferstd110))

    

    filew.close()

def draw_deviation_map(direc, filename, phase, quave, quave_angle_std):
    
    file_oimw = open(direc + filename[:-4] + '_qu.txt', 'r+')
    file_oimr= open(filename, 'r')
    file_oimw3 = open(direc + filename[:-4] + '_deviation_map.ang', 'w')
    quavei= Qu.qu_inv(quave)
    quaveaus_angle_std = []
    

    file_oimw.seek(0) # 重新打檔案

    for qustr, prnt in zip(file_oimw.readlines(), file_oimr.readlines()):
        if qustr.startswith('#'):
            # 寫入開頭檔案
            file_oimw3.write(qustr)
            continue

        qustr = qustr.replace('[','')
        qustr = qustr.replace(']','')
        qustr = qustr.replace(',','')
        qustr = qustr.strip()
        qulist = [float(qustri) for qustri in qustr.split()]

        if qulist[4] == phase['Austenite']:

            quaus = Qu.qu_mult(qulist[0:4], quavei)
            # 乘上inversed quave會把每一筆austenite轉回100 010 001標準位置

            aus111angle = Qu.aus111aus111_deviation (quaus) * 100 # 100 for output
            aus110angle = Qu.aus110aus110_deviation (quaus) * 100 # 100 for output

            devang = Qu.qu_angle(qulist[0:4], quave) * 100  # 100 for output
            
            line = oim_data_output(qulist[0:4], prnt, devang)
            file_oimw3.write(line)

            if aus111angle > 700 or aus110angle > 700:
                continue

            quave_angle_std.append(devang/100)

            continue

        elif qulist[4] == phase['Ferrite']:

            qufer = Qu.qu_mult(qulist[0:4], quavei)
            qufer = Qu.qu_std(qufer)

            var_num, misang = Qu.qu_fervar_num_from_gt(qufer)
            aus111angle = Qu.aus111fer110_deviation(qufer, var_num)
            aus110angle = Qu.aus110fer111_deviation(qufer, var_num)

            misang = misang * 100
            # IQ-111   SEM-110
            

            # 我乘了qulist，應該是要先作會標準！
            qufer = Qu.qu_mult(qufer, GTitoGT1_cubic_corres[var_num])
            qufer = Qu.qu_std(qufer)

            # print the adjucted grain
            quferprnt = Qu.qu_mult(qufer, Qu.qu_inv(GTitoGT1_cubic_corres[var_num]))
            quferprnt = Qu.qu_mult(quferprnt, quave)


            line = oim_data_output(quferprnt, prnt,\
                        SEM_signal = misang)
            file_oimw3.write(line)


            # 要來計算平均肥粒鐵方位：> 7度的都不算進去了：
            if aus111angle > 700 or aus110angle > 700:
                continue

            quaveaus_angle_std.append(misang/100)

            continue

    return quaveaus_angle_std, quave_angle_std


    
def draw_dev110_111_map(direc, filename, phase, quave):
    
    file_oimw = open(direc + filename[:-4] + '_qu.txt', 'r+')
    file_oimr= open(filename, 'r')
    file_oimw2 = open(direc + filename[:-4] + '_dev110_111_map.ang', 'w')
    quavei= Qu.qu_inv(quave)
    quavefer = [0.0,0.0,0.0,0.0]
    ferrite_grains = 0
    ferrite_grains_analyzed = 0
    variants_counts = [0] * 24
    quaveaus_angle110_std = []
    quaveaus_angle111_std = []
    
    file_oimw.seek(0) # 重新打檔案


    for qustr, prnt in zip(file_oimw.readlines(), file_oimr.readlines()):
        if qustr.startswith('#'):
            # 寫入開頭檔案
            file_oimw2.write(qustr)
            continue

        qustr = qustr.replace('[','')
        qustr = qustr.replace(']','')
        qustr = qustr.replace(',','')
        qustr = qustr.strip()
        qulist = [float(qustri) for qustri in qustr.split()]

        if qulist[4] == phase['Austenite']:
    
            quaus = Qu.qu_mult(qulist[0:4], quavei)
            # 乘上inversed quave會把每一筆austenite轉回100 010 001標準位置

            # 第一個file
            aus111angle = Qu.aus111aus111_deviation (quaus) * 100 # 100 for output
            aus110angle = Qu.aus110aus110_deviation (quaus) * 100 # 100 for output

            line = oim_data_output(qulist[0:4], prnt,\
                        SEM_signal = aus110angle, \
                        IQ = aus111angle)
            file_oimw2.write(line)

        elif qulist[4] == phase['Ferrite']:
            ferrite_grains += 1

            qufer = Qu.qu_mult(qulist[0:4], quavei)
            qufer = Qu.qu_std(qufer)

            var_num, misang = Qu.qu_fervar_num_from_gt(qufer)
            aus111angle = Qu.aus111fer110_deviation(qufer, var_num)
            aus110angle = Qu.aus110fer111_deviation(qufer, var_num)

            # IQ-111   SEM-110
            aus111angle = aus111angle * 100
            aus110angle = aus110angle * 100
            fit     = float(var_num / 10)
            variants_counts [var_num - 1] += 1

            # 我乘了qulist，應該是要先作會標準！
            qufer = Qu.qu_mult(qufer, GTitoGT1_cubic_corres[var_num])
            qufer = Qu.qu_std(qufer)

            # print the adjucted grain
            quferprnt = Qu.qu_mult(qufer, Qu.qu_inv(GTitoGT1_cubic_corres[var_num]))
            quferprnt = Qu.qu_mult(quferprnt, quave)

            line = oim_data_output(quferprnt, prnt,\
                        SEM_signal = aus110angle, \
                        IQ = aus111angle, \
                        fit = fit)
            file_oimw2.write(line)

            # 要來計算平均肥粒鐵方位：> 7度的都不算進去了：
            if aus111angle > 700 or aus110angle > 700:
                continue

            quavefer = Qu.qu_add(quavefer,qufer)

            quaveaus_angle110_std.append(aus110angle/100)
            quaveaus_angle111_std.append(aus111angle/100)

            ferrite_grains_analyzed += 1

            continue
    
    quavefer = quavefer / np.linalg.norm(quavefer)

    print('Ferrite average orientation:', quavefer)
    print ('Input Ferrite grains:', ferrite_grains)
    print ('Analyzed Ferrite grains:', ferrite_grains_analyzed)
    

    return \
    variants_counts, \
    ferrite_grains, \
    ferrite_grains_analyzed, \
    quavefer, quaveaus_angle110_std, \
    quaveaus_angle111_std 

def calc_aus_average(direc, filename, refgrain, phase):
    file_oimw= open(direc + filename[:-4] + '_qu.txt', 'r+')
    file_oimw.seek(0)

    qu_sum = [0.0, 0.0, 0.0, 0.0]
    qusum = [0.0,0.0,0.0,0.0]

    quave = refgrain

    meandev = 15
    itr = 3
    iter_dev = meandev / itr
    quave_angle_std = []
    while (meandev > iter_dev):

        for qustr in file_oimw.readlines():
            if qustr.startswith('#'):
                continue

            qustr = qustr.replace('[','')
            qustr = qustr.replace(']','')
            qustr = qustr.replace(',','')
            qustr = qustr.strip()
            qulist = [float(qustri) for qustri in qustr.split()]

            if qulist[4] == phase['Austenite']:
                devang = Qu.qu_angle(qulist[0:4], quave)
                if devang < (iter_dev) :
                    qusum = Qu.qu_add(qulist[0:4], qusum)
                    if itr == 1:
                        quave_angle_std.append(devang)
        quave = qusum / np.linalg.norm(qusum)
        itr -= 1
        iter_dev = meandev / itr
        # 優化austenite grain的平均

    return quave , quave_angle_std  

def create_qu_phase(direc, filename, filename_aus, phase):
    
    path_exist = path.isfile(direc + filename[:-4] + '_qu.txt')
    
    if not path_exist: # 如果檔案不存在，就代表還沒轉換，建立一個含有qu以及phase資訊的ang檔案
        print('Create a folder containing {0}'.format(filename))
        makedirs(direc)
        file_oimr = open(filename, 'r')
        file_oimw = open(direc + filename[:-4] + '_qu.txt', 'w')
        '''
        第一步，作qu standardization當中的autneite 整理：
        q0 for austenite 必須一次到位
        q0 for ferrite 則需要轉換
        找到作為平均值的參考方位 refgrain
        '''
        refgrain = austenite_mean_grain(filename_aus)
        print('Reference austenite (in quaternion): ', refgrain)
        file_oimw_mean_grain= open(direc + filename[:-4] + '_qu_mean.txt', 'w')
        file_oimw_mean_grain.write('{0}'.format(refgrain))
        file_oimw_mean_grain.close

        ausgrain_in_oim = 0

        for line in file_oimr.readlines():
            if line.startswith('#'):
                file_oimw.write(line)
                continue

            line = line.strip()
            eulerangles = line.split()

            if int(eulerangles[7]) == -1:
                qu = [0,0,0,0]
                #如果是-1的，代表這個grain不重要不必進行運算。
                #但是還是必須寫入qu

            elif int(eulerangles[7]) == phase['Ferrite']:
                euler = [float (phis) for phis in eulerangles[0:3]]
                qu = Eu.euler2qu(euler)
                # ferrite僅作轉換即可，因為接下來才能正式作轉換


            else:
                euler = [float (phis) for phis in eulerangles[0:3]]
                qu = Eu.euler2qu(euler)
                qu = Qu.qu_std(qu, refgrain)
                devan = Qu.qu_angle(qu, refgrain)
                # 將其與ref grain作標準化
                # if devan > 15:
                #     eulerangles[7] = -1
                #     qu = [0, 0, 0, 0]
                #     continue

                ausgrain_in_oim += 1

            file_oimw.write('{0}  {1:d}  \n'.format(qu, int(eulerangles[7])))
    
    else: # 如果檔案存在，直接進行計算。
        file_oimw_mean_grain= open(direc + filename[:-4] + '_qu_mean.txt', 'r')
        line = str(file_oimw_mean_grain.readline())
        line = line.replace('[','')
        line = line.replace(']','')
        refgrains = line.split()
        refgrain = [float(ref) for ref in refgrains]
        file_oimw_mean_grain.close
        print('Reference austenite:', refgrain)
    
    return refgrain
 

def know_oim_phase(filename):
    file_oimr= open(filename, 'r')
    idn = 0
    phaseid = {}
    for line in file_oimr.readlines():
        if line.startswith('#') == False:
            break
        if idn > 0:
            phase = line.split()
            phaseid[phase[2]] = idn
            idn = 0
        if line.startswith('#') and ('Phase' in line):
            number = line.split()
            idn = int(number[2])
        else:
            continue
    return (phaseid)

def austenite_mean_grain(filename_aus):
    itr = 1
    file_oimr = open(filename_aus,'r')

    for line in file_oimr.readlines():

        if line.startswith('#'):
            continue

        line = line.strip()
        eulerangles = line.split()

        euler = [float (phis) for phis in eulerangles[0:3]]
        qu = Eu.euler2qu(euler)

        if itr == 1:
            itr += 1
            qusum = Qu.qu_std(qu)
            qu1 = Qu.qu_std(qu)
            continue

        qu = Qu.qu_std(qu, qu1)
        qusum = Qu.qu_add(qusum, qu)

    quave = qusum / np.linalg.norm(qusum)
    return quave



def oim_data_output (qu, oimline = 'without input', SEM_signal = 0, IQ = 0, fit = 0):
    '''
    將檔案輸出成OIM的樣子，已使得OIM可以讀取。
    oimline與qu之指稱之grain應相同。
    SEM_signal: 代表沃斯田鐵的dev. angle

    對肥粒鐵:
    SEM_signal: 代表與110的角度偏差*100
    IQ: 代表與111的角度偏差*100
    CI: 代表Variant Number / 24

    要*100是因為OIM的SEM與IQ只讀取整數
    要/10，此外fit只讀取到小數第一位
    '''
    global x_coordinate_of_oim
    # 要在其他function之前，否則會出現error

    # except在第一次發生，之後都發生try
    try:
        if oimline == 'without input':
            oimline = '2.57290   0.31547   3.04730  ' +\
                            str(x_coordinate_of_oim) + \
                                    '.00000 0.00000 100 1 0 1 1.355'
            x_coordinate_of_oim = x_coordinate_of_oim + 1
    except:
        x_coordinate_of_oim = 0
        if oimline == 'without input':
            oimline = '2.57290   0.31547   3.04730  ' +\
                            str(x_coordinate_of_oim) + \
                                    '.00000 0.00000 100 1 0 1 1.355'
            x_coordinate_of_oim = x_coordinate_of_oim + 1

    oimline_str = oimline.split()
    euler = Qu.qu2euler(qu)

    oimline_str[0] = '{0:7.5f}'.format(euler[0])
    oimline_str[1] = '{0:7.5f}'.format(euler[1])
    oimline_str[2] = '{0:7.5f}'.format(euler[2])

    if SEM_signal != 0:
        oimline_str[8] = '{0:4.0f}'.format(SEM_signal)

    if IQ != 0:
        oimline_str[5] = '{0:4.0f}'.format(IQ)

    if fit != 0:
        oimline_str[9] = '{0:3.1f}'.format(fit)

    line = ('  ').join(oimline_str) +'\n'

    return line


def main(filename, filename_aus):
    
    cur_path = getcwd()
    direc = cur_path + "\\" + filename[:-4] + '\\'
    print('File path:', cur_path)
    print ('Start time:', time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) )
    starttime = time.time()
    print('File name:', filename)
    
    # 建立phase_id的dictionary
    phase = know_oim_phase(filename) 

    # 建立一個txt含有qu + phase的一個txt檔案來運用，順便回傳austenite方位給程式參考
    refgrain = create_qu_phase(direc, filename, filename_aus, phase)
    
    # 把refgrain餵入，回傳檔案中的austenite平均方位以及角度的偏差值
    quave, quave_angle_std = calc_aus_average(direc, filename, refgrain, phase)

    '''
    輸出檔案，順便回傳:
    總共的ferrite資料點數量、分析的ferrite資料點數量、ferrite方位資料點的平均(qu)
    每個austenite資料點與平均的110及111偏差角(list)、每個variant的數量統計
    '''
    variants_counts, \
    ferrite_grains,\
    ferrite_grains_analyzed,\
    quavefer, \
    quaveaus_angle110_std, \
    quaveaus_angle111_std = draw_dev110_111_map(direc, filename, phase, quave)
    
    '''
    畫出角度差的map，回傳:
    每一筆austenite的偏差角(list)
    每一筆ferrite的偏差角(list)
    '''
    quaveaus_angle_std, \
    quave_angle_std = draw_deviation_map(direc, filename, phase, quave, quave_angle_std)


    # 畫出: 轉成 100 010 001方位的 pole figure
    quavefer_original = Qu.qu_mult(quavefer, quave)
    Polefig.origin_polefig (quavefer_original, quave, direc)
    
    # 以平均方位來看的 pole figure (原始的)
    plot_OIM_OR(direc, quave, quavefer, quavefer_original, quave_angle_std, \
    quaveaus_angle_std, quaveaus_angle110_std, quaveaus_angle111_std)
    
    # 把dominanting的variant轉成V1，並統計24個variant的grain數
    output_variant_distribution(direc, filename, variants_counts)

    # 畫出平均austenite方位的pole figure
    OIMExport.aus_oim_or_plot(quave, direc)

    # 以平均austenite方位為準來模擬NW的方位關係
    OIMExport.nw_oim_or_plot(quave, direc)

    # 以平均austenite方位為準來模擬KS的方位關係
    OIMExport.ks_oim_or_plot(quave, direc)
    
    # 算出austenite晶體下的schmid factors，並回傳list以及dic
    schmid_factors_list, schmid_factors_dic = know_schmid_factors(direc, filename, quave)
    
    # 畫出schmid factor的分布圖
    plot_schmid_factors(direc, filename, schmid_factors_dic)
    
    # 畫出因為deformation造成的variant selection的pole figure
    plot_variant_selection(direc, schmid_factors_dic, quave)

    endtime = time.time()
    print ('it costs ',  endtime - starttime  , ' s ')
    print('Done!')


if __name__ == '__main__' :
    
    filedir = (path.dirname(path.abspath(__file__)))
    angs = [f for f in listdir(filedir) if f.endswith('.ang')]
    ang_aus = [f for f in listdir(filedir) if f.endswith('aus.txt')]

    for ang in angs:
        if ang == 'head.ang':
            continue
        if ang in ang_aus:
            continue
        main(ang, ang[:-4] + '_aus.txt')
   
