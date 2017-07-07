from os import path, makedirs, listdir
import numpy as np
from crystal_math import Qu
from crystal_math import OIMExport, Polefig
from crystal_math import cubicsymmetry
from crystal_math import Euler as Eu
from crystal_math import SpecialTrans as Sp
import math
import time

sqrt2i = 1 / math.sqrt(2)
sqrt3i = 1 / math.sqrt(3)
GTitoGT1_cubic_corres = {
		1 : [1,0,0,0],\
		2 : [  0,-sqrt2i,   0, sqrt2i],\
		3 : [0.5, 0.5, 0.5, 0.5],\
		4 : [  0,-sqrt2i, sqrt2i,   0],\
		5 : [0.5,-0.5,-0.5,-0.5],\
		6 : [  0,   0,-sqrt2i, sqrt2i],\
		7 : [  0, sqrt2i,   0, sqrt2i],\
		8 : [0, 0, 1, 0],\
		9 : [sqrt2i,-sqrt2i,   0,   0],\
		10: [0.5, 0.5, 0.5,-0.5],\
		11: [sqrt2i,   0,   0, sqrt2i],\
		12: [0.5, 0.5,-0.5,-0.5],\
		13: [sqrt2i,   0,   0,-sqrt2i],\
		14: [0.5,-0.5,-0.5, 0.5],\
		15: [sqrt2i,   0, sqrt2i,   0],\
		16: [0, 1, 0, 0],\
		17: [  0,   0, sqrt2i, sqrt2i],\
		18:	[0.5, 0.5,-0.5, 0.5],\
		19: [sqrt2i, sqrt2i,   0,   0],\
		20: [0.5,-0.5, 0.5, 0.5],\
		21: [  0, sqrt2i, sqrt2i,   0],\
		22: [0.5,-0.5, 0.5,-0.5],\
		23: [sqrt2i,   0,-sqrt2i,   0],\
		24: [0, 0, 0, 1]	}


# Qu這個class已經被import，可直接使用:
# 如 Qu.qu_inv

def know_oim_phase(filename):
    '''
    輸出為一個dictionary，
    phaseid ['Austenite'] = 1 (integer)
    '''
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

    oimline_str[0] = '   {0:7.5f}'.format(euler[0])
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


def oim_eulers_to_mean_grain(filename, filename_aus):
    '''
    以用沃斯田鐵平均出發，作q0校正
    '''
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) )
    starttime = time.time()

    direc = r'C:\Users\\Andrew Tung\Desktop\Python Crystal\\' + filename[:-4] + '\\'
    phase = know_oim_phase(filename)
    print(filename)

    pathexist = path.isfile(direc + filename[:-4] + '_qu.txt')
    qusum = [0.0,0.0,0.0,0.0]
    if not pathexist :
        '''
        如果檔案存在時，代表已進行轉換
        '''
        makedirs(direc)
        file_oimr= open(filename, 'r')
        file_oimw= open(direc + filename[:-4] + '_qu.txt', 'w')

        # 第一步，作qu standardization當中的autneite 整理：
        # q0 for austenite 必須一次到位
        # q0 for ferrite 則需要轉換。
        ferrite_grains = 0   # 用以計算沃斯田鐵的顆粒數
        ferrite_qus = []
            # 用以計算沃斯田鐵的暫訂平均
        try_time = 0
        tempave = 'not found'

        '''
        已找到作為平均值的參考方位 refgrain
        '''
        refgrain = austenite_mean_grain(filename_aus)
        print (refgrain)
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

    else:
        file_oimw_mean_grain= open(direc + filename[:-4] + '_qu_mean.txt', 'r')
        line = str(file_oimw_mean_grain.readline())
        line = line.replace('[','')
        line = line.replace(']','')
        refgrains = line.split()
        refgrain = [float(ref) for ref in refgrains]
        file_oimw_mean_grain.close
        print (refgrain)


    # 如果檔案存在，直接進行計算。
    file_oimw= open(direc + filename[:-4] + '_qu.txt', 'r+')
    file_oimw.seek(0)

    qusum = [0.0,0.0,0.0,0.0]

    quaveaus_angle110_std = []
    quaveaus_angle111_std = []
    quaveaus_angle_std = []

    quave = refgrain
    # iteration 1
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



    quavei= Qu.qu_inv(quave)
    quavefer = [0.0,0.0,0.0,0.0]
    ferrite_grains = 0
    ferrite_grains_analyzed = 0
    variants_counts = [0] * 24

    # 輸出檔案，順便計算該KS的variant number! 首先重新打開文件，這次要讀非-1的
    file_oimw.seek(0)
    file_oimw2 = open(direc + filename[:-4] + '_dev110_111_map.ang', 'w')
    file_oimw3 = open(direc + filename[:-4] + '_deviation_map.ang', 'w')
    file_oimr= open(filename, 'r')

    # 用zip function將此二list 作榜定
    for qustr, prnt in zip(file_oimw.readlines(), file_oimr.readlines()):
        if qustr.startswith('#'):
            # 寫入開頭檔案
            file_oimw2.write(qustr)
            file_oimw3.write(qustr)
            continue

        qustr = qustr.replace('[','')
        qustr = qustr.replace(']','')
        qustr = qustr.replace(',','')
        qustr = qustr.strip()
        qulist = [float(qustri) for qustri in qustr.split()]

        if qulist[4] == phase['Austenite']:

            quaus = Qu.qu_mult(qulist[0:4], quavei)

            # 第一個file
            aus111angle = Qu.aus111aus111_deviation (quaus) * 100 # 100 for output
            aus110angle = Qu.aus110aus110_deviation (quaus) * 100 # 100 for output

            line = oim_data_output(qulist[0:4], prnt,\
                        SEM_signal = aus110angle, \
                        IQ = aus111angle)
            file_oimw2.write(line)


            # 第二個file
            devang = Qu.qu_angle(qulist[0:4], quave) * 100  # 100 for output


            line = oim_data_output(qulist[0:4], prnt, devang)
            file_oimw3.write(line)

            if aus111angle > 700 or aus110angle > 700:
                continue

            quave_angle_std.append(devang/100)

            continue

        elif qulist[4] == phase['Ferrite']:
            ferrite_grains += 1

            qufer = Qu.qu_mult(qulist[0:4], quavei)
            qufer = Qu.qu_std(qufer)

            var_num, misang = Qu.qu_fervar_num_from_gt(qufer)
            aus111angle = Qu.aus111fer110_deviation(qufer, var_num)
            aus110angle = Qu.aus110fer111_deviation(qufer, var_num)

            misang = misang * 100
            # IQ-111   SEM-110
            aus111angle = aus111angle * 100
            aus110angle = aus110angle * 100
            fit     = float(var_num / 10)
            variants_counts [var_num - 1] += 1

                #  我乘了qulist，應該是要先作會標準！
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

            line = oim_data_output(quferprnt, prnt,\
                        SEM_signal = misang)
            file_oimw3.write(line)


            # 要來計算平均肥粒鐵方位：> 7度的都不算進去了：
            if aus111angle > 700 or aus110angle > 700:
                continue

            quavefer = Qu.qu_add(quavefer,qufer)

            quaveaus_angle110_std.append(aus110angle/100)
            quaveaus_angle111_std.append(aus111angle/100)
            quaveaus_angle_std.append(misang/100)

            ferrite_grains_analyzed += 1

            continue

        line = oim_data_output(qulist[0:4], prnt)
        file_oimw2.write(line)
        file_oimw3.write(line)

    # 計算平均方位等之差異
    quavefer = quavefer / np.linalg.norm(quavefer)
    var_num, misang = Qu.qu_fervar_num_from_gt(quavefer)
    aus111angle = Qu.aus111fer110_deviation(quavefer, var_num)
    aus110angle = Qu.aus110fer111_deviation(quavefer, var_num)

	# 計算標準差，注意degree of freedom = 1，代表N-1
    quavestd = np.std(quave_angle_std, ddof=1)
    quaveferstd = np.std(quaveaus_angle_std, ddof=1)
    quaveferstd110 = np.std(quaveaus_angle110_std, ddof=1)
    quaveferstd111 = np.std(quaveaus_angle111_std, ddof=1)

    print ('Input Ferrite grains:', ferrite_grains)
    print ('Analyzed grains     :', ferrite_grains_analyzed)

    quavefer = Qu.qu_mult(quavefer, quave)

    Polefig.origin_polefig (quavefer, quave, direc)
    filew = OIMExport.oim_or_plot (quavefer, quave, direc)
    filew.write('\n\nOR')
    filew.write('{0}'.format(quavefer))

    filew.write('\ndeviation of austenite mean grain: {0:4.2f}'.format(quavestd))
    filew.write('\ndeviation of ferrite mean grain:   {0:4.2f}'.format(quaveferstd))

    filew.write('\nDeviated From Aus111: {0:4.4f}'.format(aus111angle))
    filew.write('\t deviation:  {0:4.2f}'.format(quaveferstd111))

    filew.write('\nDeviated From Aus110: {0:4.4f}'.format(aus110angle))
    filew.write('\t deviation:  {0:4.2f}'.format(quaveferstd110))

    variant_sum = np.sum(variants_counts)

    filew.close()
    filevar = open (direc + filename[:-4] + '_var_distribution.txt', 'w')

    variants_counts = Sp.variant24_convert(variants_counts)

    for index, var_num in enumerate(variants_counts):
        strprint = '{0:2d}\t {1:6.4f} \n'\
                        .format(index + 1, var_num / variant_sum)
        filevar.write(strprint)

    endtime = time.time()
    print ('it costs ',  endtime - starttime  , ' s ')


if __name__ == '__main__' :

    filedir = (path.dirname(path.abspath(__file__)))
    angs = [f for f in listdir(filedir) if f.endswith('.ang')]
    ang_aus = [f for f in listdir(filedir) if f.endswith('aus.txt')]

    for ang in angs:
        if ang == 'head.ang':
            continue
        if ang in ang_aus:
            continue
        oim_eulers_to_mean_grain(ang, ang[:-4] + '_aus.txt')
