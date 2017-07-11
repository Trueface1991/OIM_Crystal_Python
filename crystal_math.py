'''
this file contains the basic mathematical methods:
'''
import numpy as np
from numpy import dot
import math as math
from shutil import copyfile

global cubicsymmetry
sqrt2i = 1 / math.sqrt(2)
sqrt3i = 1 / math.sqrt(3)
cubicsymmetry = [[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0],
                 [0, 0, 0, 1],
                 [0.5, 0.5, 0.5, 0.5],
                 [0.5,-0.5,-0.5,-0.5],
                 [0.5, 0.5,-0.5, 0.5],
                 [0.5,-0.5, 0.5,-0.5],
                 [0.5,-0.5, 0.5, 0.5],
                 [0.5, 0.5,-0.5,-0.5],
                 [0.5,-0.5,-0.5, 0.5],
                 [0.5, 0.5, 0.5,-0.5],
                 [sqrt2i, sqrt2i,   0,   0],
                 [sqrt2i,   0, sqrt2i,   0],
                 [sqrt2i,   0,   0, sqrt2i],
                 [sqrt2i,-sqrt2i,   0,   0],
                 [sqrt2i,   0,-sqrt2i,   0],
                 [sqrt2i,   0,   0,-sqrt2i],
                 [  0, sqrt2i, sqrt2i,   0],
                 [  0,-sqrt2i, sqrt2i,   0],
                 [  0,   0, sqrt2i, sqrt2i],
                 [  0,   0,-sqrt2i, sqrt2i],
                 [  0, sqrt2i,   0, sqrt2i],
                 [  0,-sqrt2i,   0, sqrt2i]]


'''gt was adjusted according to KS'''
global gtdict
# GTquDict_KSorder: GT in KS order
gtdict = {\
		1 : [0.9264174573, -0.0710698455, 0.0500952943, -0.3663198513],\
		2 : [0.9264174573, 0.3663198513, -0.0500952943, 0.0710698455],\
		3 : [0.9264174573, 0.0500952943, -0.3663198513, -0.0710698455],\
		4 : [0.9264174573, -0.0500952943, 0.0710698455, 0.3663198513],\
		5 : [0.9264174573, -0.3663198513, -0.0710698455, 0.0500952943],\
		6 : [0.9264174573, 0.0710698455, 0.3663198513, -0.0500952943],\
		7 : [0.9264174573, -0.3663198513, -0.0500952943, -0.0710698455],\
		8 : [0.9264174573, 0.0710698455, 0.0500952943, 0.3663198513],\
		9 : [0.9264174573, -0.0710698455, 0.3663198513, 0.0500952943],\
		10 : [0.9264174573, 0.3663198513, -0.0710698455, -0.0500952943],\
		11 : [0.9264174573, 0.0500952943, 0.0710698455, -0.3663198513],\
		12 : [0.9264174573, -0.0500952943, -0.3663198513, 0.0710698455],\
		13 : [0.9264174573, -0.0500952943, -0.0710698455, -0.3663198513],\
		14 : [0.9264174573, 0.0500952943, 0.3663198513, 0.0710698455],\
		15 : [0.9264174573, 0.3663198513, 0.0500952943, -0.0710698455],\
		16 : [0.9264174573, -0.0710698455, -0.0500952943, 0.3663198513],\
		17 : [0.9264174573, 0.0710698455, -0.3663198513, 0.0500952943],\
		18 : [0.9264174573, -0.3663198513, 0.0710698455, -0.0500952943],\
		19 : [0.9264174573, -0.0710698455, -0.3663198513, -0.0500952943],\
		20 : [0.9264174573, 0.3663198513, 0.0710698455, 0.0500952943],\
		21 : [0.9264174573, 0.0500952943, -0.0710698455, 0.3663198513],\
		22 : [0.9264174573, -0.0500952943, 0.3663198513, -0.0710698455],\
		23 : [0.9264174573, -0.3663198513, 0.0500952943, 0.0710698455],\
		24 : [0.9264174573, 0.0710698455, -0.0500952943, -0.3663198513],\
}

global ksdict
ksdict = {\
		1 : [0.9309036628, -0.0648782514, 0.0648782514, -0.3535533844],\
		2 : [0.9309036597, 0.3535533904, -0.0648782571, 0.0648782578],\
		3 : [0.9309036597, 0.0648782571, -0.3535533904, -0.0648782578],\
		4 : [0.9309036628, -0.0648782514, 0.0648782514, 0.3535533844],\
		5 : [0.9309036588, -0.3535533927, -0.0648782580, 0.0648782567],\
		6 : [0.9309036588, 0.0648782580, 0.3535533927, -0.0648782567],\
		7 : [0.9309036597, -0.3535533904, -0.0648782571, -0.0648782578],\
		8 : [0.9309036628, 0.0648782514, 0.0648782514, 0.3535533844],\
		9 : [0.9309036588, -0.0648782580, 0.3535533927, 0.0648782567],\
		10 : [0.9309036588, 0.3535533927, -0.0648782580, -0.0648782567],\
		11 : [0.9309036628, 0.0648782514, 0.0648782514, -0.3535533844],\
		12 : [0.9309036597, -0.0648782571, -0.3535533904, 0.0648782578],\
		13 : [0.9309036628, -0.0648782514, -0.0648782514, -0.3535533844],\
		14 : [0.9309036597, 0.0648782571, 0.3535533904, 0.0648782578],\
		15 : [0.9309036597, 0.3535533904, 0.0648782571, -0.0648782578],\
		16 : [0.9309036628, -0.0648782514, -0.0648782514, 0.3535533844],\
		17 : [0.9309036588, 0.0648782580, -0.3535533927, 0.0648782567],\
		18 : [0.9309036588, -0.3535533927, 0.0648782580, -0.0648782567],\
		19 : [0.9309036588, -0.0648782580, -0.3535533927, -0.0648782567],\
		20 : [0.9309036588, 0.3535533927, 0.0648782580, 0.0648782567],\
		21 : [0.9309036628, 0.0648782514, -0.0648782514, 0.3535533844],\
		22 : [0.9309036597, -0.0648782571, 0.3535533904, -0.0648782578],\
		23 : [0.9309036597, -0.3535533904, 0.0648782571, 0.0648782578],\
		24 : [0.9309036628, 0.0648782514, -0.0648782514, -0.3535533844],\
}

global nwdict
nwdict = {\
		1 : [0.9205472246, -0.3813031426, -0.0783977010, 0.0324733912],\
		2 : [0.9205472246, 0.0783977010, 0.3813031426, -0.0324733912],\
		3 : [0.9205472240, -0.0783976975, 0.0324733896, -0.3813031450],\
		4 : [0.9205472246, -0.3813031426, 0.0783977010, -0.0324733912],\
		5 : [0.9205472246, 0.0783977010, -0.3813031426, 0.0324733912],\
		6 : [0.9205472240, -0.0783976975, -0.0324733896, 0.3813031450],\
		7 : [0.9205472246, 0.3813031426, -0.0783977010, -0.0324733912],\
		8 : [0.9205472246, -0.0783977010, 0.3813031426, 0.0324733912],\
		9 : [0.9205472240, 0.0783976975, 0.0324733896, 0.3813031450],\
		10 : [0.9205472254, 0.3813031416, 0.0783976989, 0.0324733873],\
		11 : [0.9205472254, -0.0783976989, -0.3813031416, -0.0324733873],\
		12 : [0.9205472239, 0.0324733902, -0.0783976990, 0.3813031449],\
}

global GTitoGT1_cubic_corres
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

''' for calculating aus110_deviation'''
global aus110dict
aus110dict = {          1 : [sqrt2i, 0, -sqrt2i] ,\
                        2 : [-sqrt2i, 0, sqrt2i] ,\
                        3 : [0, -sqrt2i, sqrt2i] ,\
                        4 : [0, sqrt2i, -sqrt2i] ,\
                        5 : [-sqrt2i, sqrt2i, 0] ,\
                        6 : [sqrt2i, -sqrt2i, 0] ,\
                        7 : [sqrt2i, 0, -sqrt2i] ,\
                        8 : [-sqrt2i, 0, sqrt2i] ,\
                        9 : [sqrt2i, sqrt2i, 0] ,\
                        10 : [sqrt2i, sqrt2i, 0] ,\
                        11 : [0, sqrt2i, sqrt2i] ,\
                        12 : [0, sqrt2i, sqrt2i] ,\
                        13 : [0, -sqrt2i, sqrt2i] ,\
                        14 : [0, sqrt2i, -sqrt2i] ,\
                        15 : [sqrt2i, 0, sqrt2i] ,\
                        16 : [sqrt2i, 0, sqrt2i] ,\
                        17 : [sqrt2i, sqrt2i, 0] ,\
                        18 : [sqrt2i, sqrt2i, 0] ,\
                        19 : [-sqrt2i, sqrt2i, 0] ,\
                        20 : [sqrt2i, -sqrt2i, 0] ,\
                        21 : [0, sqrt2i, sqrt2i] ,\
                        22 : [0, sqrt2i, sqrt2i] ,\
                        23 : [sqrt2i, 0, sqrt2i] ,\
                        24 : [sqrt2i, 0, sqrt2i]}
global fer111dict
fer111dict ={          1 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        2 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        3 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        4 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        5 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        6 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        7 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        8 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        9 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        10 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        11 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        12 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        13 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        14 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        15 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        16 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        17 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        18 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        19 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        20 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        21 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        22 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        23 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        24 : [sqrt3i, sqrt3i, sqrt3i]}

''' for calculating aus111_deviation'''
global aus111dict
aus111dict= {           1 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        2 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        3 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        4 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        5 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        6 : [sqrt3i, sqrt3i, sqrt3i] ,\
                        7 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        8 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        9 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        10 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        11 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        12 : [sqrt3i, -sqrt3i, sqrt3i] ,\
                        13 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        14 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        15 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        16 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        17 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        18 : [-sqrt3i, sqrt3i, sqrt3i] ,\
                        19 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        20 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        21 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        22 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        23 : [sqrt3i, sqrt3i, -sqrt3i] ,\
                        24 : [sqrt3i, sqrt3i, -sqrt3i]}
global fer110dict
fer110dict= {           1 : [0, sqrt2i, sqrt2i] ,\
                        2 : [sqrt2i, sqrt2i, 0] ,\
                        3 : [sqrt2i, sqrt2i, 0] ,\
                        4 : [sqrt2i, 0, sqrt2i] ,\
                        5 : [sqrt2i, 0, sqrt2i] ,\
                        6 : [0, sqrt2i, sqrt2i] ,\
                        7 : [sqrt2i, -sqrt2i, 0] ,\
                        8 : [0, -sqrt2i, sqrt2i] ,\
                        9 : [0, -sqrt2i, sqrt2i] ,\
                        10 : [sqrt2i, 0, sqrt2i] ,\
                        11 : [sqrt2i, 0, sqrt2i] ,\
                        12 : [sqrt2i, -sqrt2i, 0] ,\
                        13 : [-sqrt2i, 0, sqrt2i] ,\
                        14 : [-sqrt2i, sqrt2i, 0] ,\
                        15 : [-sqrt2i, sqrt2i, 0] ,\
                        16 : [0, sqrt2i, sqrt2i] ,\
                        17 : [0, sqrt2i, sqrt2i] ,\
                        18 : [-sqrt2i, 0, sqrt2i] ,\
                        19 : [0, sqrt2i, -sqrt2i] ,\
                        20 : [sqrt2i, 0, -sqrt2i] ,\
                        21 : [sqrt2i, 0, -sqrt2i] ,\
                        22 : [sqrt2i, sqrt2i, 0] ,\
                        23 : [sqrt2i, sqrt2i, 0] ,\
                        24 : [0, sqrt2i, -sqrt2i] }

variant_convert_table =np.array(\
[[  1,   2,   5,   4,   3,   6,   7,   8,  19,  14,  13,  20,  11,  10,   23,  16,  17,  22,   9,  12,  21,  18,  15,  24],\
 [  2,   1,   6,   3,   4,   5,   8,   7,  20,  13,  14,  19,  12,   9,   24,  15,  18,  21,  10,  11,  22,  17,  16,  23],\
 [  3,   6,   1,   2,   5,   4,   9,  12,  21,  18,  15,  24,   7,   8,   19,  14,  13,  20,  11,  10,  23,  16,  17,  22],\
 [  4,   5,   2,   1,   6,   3,  10,  11,  22,  17,  16,  23,   8,   7,   20,  13,  14,  19,  12,   9,  24,  15,  18,  21],\
 [  5,   4,   3,   6,   1,   2,  11,  10,  23,  16,  17,  22,   9,  12,   21,  18,  15,  24,   7,   8,  19,  14,  13,  20],\
 [  6,   3,   4,   5,   2,   1,  12,   9,  24,  15,  18,  21,  10,  11,   22,  17,  16,  23,   8,   7,  20,  13,  14,  19],\
 [  7,   8,  19,  14,  13,  20,   1,   2,   5,   4,   3,   6,  22,  17,   16,  23,  10,  11,  18,  21,  12,   9,  24,  15],\
 [  8,   7,  20,  13,  14,  19,   2,   1,   6,   3,   4,   5,  21,  18,   15,  24,   9,  12,  17,  22,  11,  10,  23,  16],\
 [  9,  12,  21,  18,  15,  24,   3,   6,   1,   2,   5,   4,  20,  13,   14,  19,   8,   7,  16,  23,  10,  11,  22,  17],\
 [ 10,  11,  22,  17,  16,  23,   4,   5,   2,   1,   6,   3,  19,  14,   13,  20,   7,   8,  15,  24,   9,  12,  21,  18],\
 [ 11,  10,  23,  16,  17,  22,   5,   4,   3,   6,   1,   2,  24,  15,   18,  21,  12,   9,  14,  19,   8,   7,  20,  13],\
 [ 12,   9,  24,  15,  18,  21,   6,   3,   4,   5,   2,   1,  23,  16,   17,  22,  11,  10,  13,  20,   7,   8,  19,  14],\
 [ 13,  20,   7,   8,  19,  14,  18,  21,  12,   9,  24,  15,   1,   2,    5,   4,   3,   6,  22,  17,  16,  23,  10,  11],\
 [ 14,  19,   8,   7,  20,  13,  17,  22,  11,  10,  23,  16,   2,   1,    6,   3,   4,   5,  21,  18,  15,  24,   9,  12],\
 [ 15,  24,   9,  12,  21,  18,  16,  23,  10,  11,  22,  17,   3,   6,    1,   2,   5,   4,  20,  13,  14,  19,   8,   7],\
 [ 16,  23,  10,  11,  22,  17,  15,  24,   9,  12,  21,  18,   4,   5,    2,   1,   6,   3,  19,  14,  13,  20,   7,   8],\
 [ 17,  22,  11,  10,  23,  16,  14,  19,   8,   7,  20,  13,   5,   4,    3,   6,   1,   2,  24,  15,  18,  21,  12,   9],\
 [ 18,  21,  12,   9,  24,  15,  13,  20,   7,   8,  19,  14,   6,   3,    4,   5,   2,   1,  23,  16,  17,  22,  11,  10],\
 [ 19,  14,  13,  20,   7,   8,  22,  17,  16,  23,  10,  11,  18,  21,   12,   9,  24,  15,   1,   2,   5,   4,   3,   6],\
 [ 20,  13,  14,  19,   8,   7,  21,  18,  15,  24,   9,  12,  17,  22,   11,  10,  23,  16,   2,   1,   6,   3,   4,   5],\
 [ 21,  18,  15,  24,   9,  12,  20,  13,  14,  19,   8,   7,  16,  23,   10,  11,  22,  17,   3,   6,   1,   2,   5,   4],\
 [ 22,  17,  16,  23,  10,  11,  19,  14,  13,  20,   7,   8,  15,  24,    9,  12,  21,  18,   4,   5,   2,   1,   6,   3],\
 [ 23,  16,  17,  22,  11,  10,  24,  15,  18,  21,  12,   9,  14,  19,    8,   7,  20,  13,   5,   4,   3,   6,   1,   2],\
 [ 24,  15,  18,  21,  12,   9,  23,  16,  17,  22,  11,  10,  13,  20,    7,   8,  19,  14,   6,   3,   4,   5,   2,   1]])

class Qu:
    @staticmethod
    def qu_inv(qu):
        '''
        quartarion 的inverse 相當於axis倒過來!
        '''
        qu = [-i for i in qu]
        qu [0] = - qu[0]
        return qu

    @staticmethod
    def qu_add(qu1, qu2):
        qusum = [qu1[i] + qu2[i] for i in range(4)]
        return qusum

    @staticmethod
    def qu_print(qu):
        qustr = ['{0:7.5f}'.format(qu[i]) for i in range(4)]
        qustr = (', ').join(qustr)
        return ('[' + qustr + ']')

    @staticmethod
    def qu_mult (qu1, qu2):
        '''
        P = -1
        notice the quanterion product

        The quaternion product of two vectors
        (qu1[1], qu1[2], qu1[3]) and (qu2[1], qu2[2], qu2[3]) is the
        product of q = qu1[1]i + qu1[2]j + zk and q' = x’i + y’j + z’k as quaternions.
        The quaternion product qq´ works out to be
        – (xx´ + yy´ + zz´) + (qu1[2]* qu2[3] –qu1[3]* qu2[2])i
                            + ( qu1[3]* qu2[1] –qu1[1]* qu2[3])j
                            + (qu1[1]* qu2[2] –qu1[2]*qu2[1])k
        '''
        # 如果用np.cross寫的話，輸出的結果為  ndarray!

        qu31 = qu1[0]*qu2[0] - dot(qu1[1:4], qu2[1:4])

        qu32 = - qu1[2]*qu2[3] + qu1[3]*qu2[2] \
                            + qu1[0]*qu2[1] + qu1[1]*qu2[0]

        qu33 = - qu1[3]*qu2[1] + qu1[1]*qu2[3] \
                            + qu1[0]*qu2[2] + qu1[2]*qu2[0]

        qu34 = - qu1[1]*qu2[2] + qu1[2]*qu2[1] \
                            + qu1[0]*qu2[3] + qu1[3]*qu2[0]

        if qu31 < 0:
            qu3 = [-qu31, -qu32, -qu33, -qu34]
        else:
            qu3 = [qu31, qu32, qu33, qu34]
        return qu3

    @staticmethod
    def qu_vec_rotation (qu, vector):
        '''
        P = -1
        Lp(r)
        'passive rotation'
        '''
        quarray = np.array(qu[1:])
        vectorarray = np.array(vector)

        Lpr = ( qu[0]**2 - np.sum(quarray**2) ) * vectorarray  \
                + 2 * dot(quarray, vectorarray) * quarray      \
                - 2 * np.cross (quarray, vectorarray) * qu[0]

        return list(Lpr)

    @classmethod
    def qu_angle (cls, qu1, qu2, smallest = 'no'):
        '''
        qu間夾角實際上要兩倍才會是晶體間旋轉的角度
        但為求方便因此直接*2，從而獲得晶體間旋轉的角度。
        若smallest是yes時，作對稱操作。

        '''
        if  'yes' in smallest or 'y' in smallest:
            qu24 = cls.qu_sym(qu1)
            qudotlist = list(map(lambda qu1i: np.dot(qu1i, qu2), qu24))
            dotpro = max(qudotlist)

        else:
            dotpro = np.dot(qu1,qu2)

        try:
            angle = 2 * math.degrees(math.acos(dotpro))
            return angle

        except ValueError:
            '''
            由於浮點數誤差的關係，可能出現dot product > 1之情形，因此
            對此二者進行排除！
            '''
            if dotpro > 1:
                return 0.0
            else:
                return 180.0

    @classmethod
    def vec_angle (cls, vec1, vec2):
        '''
        指的是一般向量之夾角。
        注意不能跟qu的角度混用，

        '''
        dotpro = np.dot(vec1,vec2)
        try:
            angle = math.degrees(math.acos(dotpro))
            return angle
        except ValueError:
            '''
            由於浮點數誤差的關係，可能出現dot product > 1之情形，因此
            對此二者進行排除！
            '''
            if dotpro > 1:
                return 0.0
            else:
                return 180.0

    @classmethod
    def qu_std (cls, qu1, qu2 = 'something'):
        '''
        將qu1轉成與qu2最接近之symmetry-related variant!
        若是qu2 == None的時候，則是找q0最大的variant
        '''
        qu1_24 = np.array(cls.qu_sym(qu1))
        if len(qu2) > 4:
        # 利用np.array 將 qu_sym得來的 map object轉換為array，接著再找首行q0的最大值
            qu1_ind = np.argmax(qu1_24, axis = 0)
            qu1 = qu1_24[qu1_ind[0],:]
            return qu1
        else:
            result = list(map(lambda qu24: dot(qu2, qu24) , qu1_24))
            result_ind = np.argmax(result)
            qu1 = qu1_24[result_ind,:]
            return qu1

    @staticmethod
    def qu2euler(qu):
        '''
        the input qu is an array
        using the atan2 function
        P = -1
        then export the euler angles in radians
        '''
        qu03 = qu[0]*qu[0] + qu[3]*qu[3]
        qu12 = qu[1]*qu[1] + qu[2]*qu[2]
        chi = math.sqrt(qu03*qu12)
        if chi == 0 and qu12 == 0:
            phi1 = math.atan2 (2*qu[0]*qu[3], qu[0]*qu[0] - qu[3]*qu[3])
            Phi  = 0
            phi2 = 0
        elif chi == 0 and qu03 == 0:
            phi1 = math.atan2 (2*qu[1]*qu[2], qu[1]*qu[1] - qu[2]*qu[2])
            Phi  = math.pi
            phi2 = 0
        else:
            phi1 = math.atan2 ((qu[1]*qu[3] + qu[0]*qu[2])/chi,   \
                               (qu[0]*qu[1] - qu[2]*qu[3])/chi)
            Phi  = math.atan2 ( 2*chi, qu03 - qu12)
            phi2 = math.atan2 ((-qu[0]*qu[2] + qu[1]*qu[3])/chi,   \
                               ( qu[2]*qu[3] + qu[0]*qu[1])/chi)

        if phi1 < 0: phi1 = phi1 + 2*math.pi
        if phi2 < 0: phi2 = phi2 + 2*math.pi

        # eliminate the absurdly tiny numbers
        if abs(phi1) < 10**-5: phil = 0
        if abs(Phi)  < 10**-5: Phi  = 0
        if abs(phi2) < 10**-5: phi2 = 0

        return [phi1, Phi, phi2]

    @classmethod
    def qu_sym (cls, qu1):
        ''' returns a list of 24 symmetry-related variants'''
        result = map(lambda cubic24: cls.qu_mult(cubic24, qu1) , cubicsymmetry)
        return list(result)



    # Belows are advanced calculations used for deriving ORs

    @classmethod
    def qu_aus_to_gt (cls, qu1):
        '''returns a dictionary, consists of 24 variants from austenite orientation
        qu1'''
        my_dict = {i: cls.qu_mult(gt, qu1) for (i, gt) in gtdict.items()}
        return my_dict

    @classmethod
    def qu_fervar_num_from_gt(cls, qu_or):
        '''
        from a calculated orientation relationship,
        determine the ferrite variant
        where qu_or == qumultiply_maxq0(ferritequ, austenitequ^-1)
        Q0 MUST BE MAX
        '''

        gtdotqu  = [dot(gt, qu_or) for (i, gt) in gtdict.items()]
        var_num  = gtdotqu.index(max(gtdotqu)) + 1      # 用GT計算variant num
        misangle = cls.qu_angle(ksdict[var_num], qu_or) # 用KS回推deviation angle

        # misangle is the deviation from KS
        return var_num, misangle

    @classmethod
    def qu_fervar_num_from_nw(cls, qu_or):
        '''
        from a calculated orientation relationship,
        determine the ferrite variant
        where qu_or == qumultiply_maxq0(ferritequ, austenitequ^-1)
        Q0 MUST BE MAX
        '''

        nwdotqu  = [dot(nw, qu_or) for (i, nw) in nwdict.items()]
        var_num  = nwdotqu.index(max(nwdotqu)) + 1      # 用GT計算variant num
        misangle = cls.qu_angle(nwdict[var_num], qu_or) # 用KS回推deviation angle

        # misangle is the deviation from KS
        return var_num, misangle


    @classmethod
    def qu_austvar_num_from_ks(cls, qu_or):
        '''
        from a calculated orientation relationship,
        determine the ferrite variant
        where qu_or == qumultiply_maxq0(ferritequ, austenitequ^-1)
        Q0 MUST BE MAX
        '''

        ksdotqu  = [dot(ks, qu_or) for (i, ks) in ksdict.items()]
        var_num  = ksdotqu.index(max(ksdotqu)) + 1      # 用GT計算variant num
        misangle = cls.qu_angle(ksdict[var_num], qu_or) # 用KS回推deviation angle

        # misangle is the deviation from KS
        return var_num, misangle


    @classmethod
    def aus111fer110_deviation (cls, qu, var_num):
        '''
        to calculate the deviation of austenite 111 from ferrite 110,
        one must know that
        for KS variant 1 - 6:
        (0 1 1)ferrite  ||   ( 1  1  1)austenite
        for KS variant 7 - 12:
        (0 1 1)ferrite  ||   ( 1 -1  1)austenite
        for KS variant 13 -18:
        (0 1 1)ferrite  ||   (-1  1  1)austenite
        for KS variant 13 -18:
        (0 1 1)ferrite  ||   ( 1  1 -1)austenite

        as a result,
        calculations must be conducted using different austenite vector
        然而會跟認知相反：
        所輸入的 [sqrt3i,sqrt3i,sqrt3i] 是沃斯田鐵方向，從得到肥粒鐵方向，再計算偏差
        '''
        fer110 = fer110dict[var_num]
        aus111 = aus111dict[var_num]

        fer110_from_aus = cls.qu_vec_rotation(qu, aus111)
        devangle = cls.vec_angle(fer110_from_aus, fer110)

        return devangle

    @classmethod
    def aus110fer111_deviation (cls, qu, var_num):
        '''
        to calculate the deviation of austenite 110 from ferrite 111,
        one must know that
        odd ferrite variant:
        ( -1 -1 1 )ferr
        even variant
        ( -1 1 -1 )ferr

        every two ausenite is different:
        for  1, 2    ( -1  0  1)aus
        for  3, 4    (  0  1 -1)aus
        for  5, 6    (  1 -1  0)aus
        for  7, 8    (  1  0 -1)aus   負的1,2
        for  9,10    ( -1 -1  0)aus
        for 11,12    (  0  1  1)aus
        for 13,14    (  0 -1  1)aus
        for 15,16    ( -1  0 -1)aus
        for 17,18    (  1  1  0)aus
        for 19,20    ( -1  1  0)aus
        for 21,22    (  0 -1 -1)aus
        for 23,23    (  1  0  1)aus

        as a result,
        calculations must be conducted using different austenite vector
        然而會跟認知相反：
        所輸入的 [sqrt3i,sqrt3i,sqrt3i] 是沃斯田鐵方向，從得到肥粒鐵方向，再計算偏差
        '''

        fer111 = fer111dict[var_num]
        aus110 = aus110dict[var_num]
        fer111_from_aus = cls.qu_vec_rotation(qu, aus110)
        devangle = cls.vec_angle(fer111_from_aus, fer111)

        return devangle

    @classmethod
    def aus111aus111_deviation (cls, qu):
        '''
        由於轉軸未知，因此找出111當中，最大的deviation!
        aus111dict中的，1, 7, 13, 19分別為四個austenite variant
        '''
        angles = list(map(lambda i: cls.vec_angle(cls.qu_vec_rotation\
                        (qu, aus111dict[i]), aus111dict[i]), [1, 7, 13, 19]))
        return max(angles)

    @classmethod
    def aus110aus110_deviation (cls, qu):
        '''
        由於轉軸未知，因此找出110當中，最大的deviation!
        aus111dict中的，9, 15, 21, 6, 7, 14分別為四個austenite variant
        '''
        angles = list(map(lambda i: cls.vec_angle(cls.qu_vec_rotation\
                  (qu, aus110dict[i]), aus110dict[i]), [9, 15, 21, 6, 7, 14]))
        return max(angles)


    # @classmethod
    # def qu_or_deviation_calc(cls, qufer, quaus):
    #     '''
    #     calculate the orientation relationship between ferritequ and austenitequ,
    #     and export the oim file for plot
    #
    #     if 'point' is as input, the exported file is included with point number.
    #
    #     the function returns the deviation angle from ideal KS
    #
    #     '''
    #         # (αJγ)=[α;v]*[γ;v]
    #         #      = qferrite * qaustenite ^ -1
    #         # q^-1 = q* ; [1,1,1,1] = [1,-1,-1,-1]
    #
    #     # qaus^-1
    #     quaus = cls.qu_inv(quaus)
    #     qu_or = cls.qu_mult(qufer, quaus)
    #     qu_or = cls.qu_std(qu_or)
    #
    #     var_num, misangle = cls.qu_fervar_num_from_gt(qu)
    #     # determine the ferrite variant:
    #
    #
    #     angles = np.empty([24,2])
    #     def orqu24(qu):
    #         qu24 = np.empty([24,4])
    #         for i in range(24):
    #             qu24[i, :] = qumultiply_maxq0(qu, cubicsymmetry[i+1])
    #             eu1, eu2, eu3 = qu2euler(qu24[i, :])
    #
    #             # austenite 111可由 [1/sqrt(2),0,1/sqrt(2)]計算而成
    #             # austenite 110可由 [1/sqrt(3),1/sqrt(3),-1/sqrt(3)]
    #             aus111 = np.matmul([0,sqrt2i,sqrt2i], euler2rotm(eu1,eu2,eu3))
    #             aus110 = np.matmul([-sqrt3i,-sqrt3i,sqrt3i],\
    #                                 euler2rotm(eu1,eu2,eu3))
    #
    #             # 計算角度；內積之後取角度
    #             aus111angle = vector_angle(aus111, [sqrt3i,sqrt3i,sqrt3i])
    #             aus110angle = vector_angle(aus110, [-sqrt2i,0,sqrt2i])
    #             angles [i, :] = [aus111angle , aus110angle]
    #
    #
    #
    #     orqu24(qu)
    #     '''
    #     the definition of function be before the function is used
    #     '''
    #     anglemin = np.amin(angles,axis=0)
    #     return anglemin, [var_num, misangle], qu
    #


    @staticmethod
    def qu2axisangle(qu, unit = 'deg'):
        '''
        In default, the export of angles are in degrees
        '''
        angle = 2 * math.acos(qu[0])

        axis = [0, 0, 0]
        axis[0] = qu[1] / math.sin(angle/2)
        axis[1] = qu[2] / math.sin(angle/2)
        axis[2] = qu[3] / math.sin(angle/2)

        if unit == 'deg' or 'd' in unit:
            angle = math.degrees(angle)

        return axis, angle

class Euler:
    @staticmethod
    def euler2qu(eulerangles):
        '''
        following the description of
        'Consistent representations of and conversions between 3D rotations'
        with P = -1 to be in accord with passive transformation

        See the appendix for formulas
        '''
        [phi1, Phi, phi2]= eulerangles
        sigma = (phi1 + phi2)/2
        delta = (phi1 - phi2)/2
        quar0 = math.cos(Phi/2)*math.cos(sigma)
        quar1 = math.sin(Phi/2)*math.cos(delta)
        quar2 = math.sin(Phi/2)*math.sin(delta)
        quar3 = math.cos(Phi/2)*math.sin(sigma)

        if quar0 < 0:
            qu = ([-quar0, -quar1, -quar2, -quar3])
        else:
            qu = ([quar0, quar1, quar2, quar3])
        return qu

    @staticmethod
    def euler2rotm(eulerangles):
        '''
        convet the euler angles to rotation matrix
        the input eulerangles is a 1 x 3 list
        the input euler angles must be in radians.
        '''

        # the import is an str; converstion of str to float
        [phi1, Phi, phi2]= eulerangles

        rotm11 =   math.cos(phi1)*math.cos(phi2) \
                 - math.sin(phi1)*math.sin(phi2)*math.cos(Phi)
        rotm12 =   math.sin(phi1)*math.cos(phi2) \
                 + math.cos(phi1)*math.sin(phi2)*math.cos(Phi)
        rotm13 =   math.sin(phi2)*math.sin(Phi)

        rotm21 = - math.cos(phi1)*math.sin(phi2)  \
                 - math.sin(phi1)*math.cos(phi2)*math.cos(Phi)
        rotm22 = - math.sin(phi1)*math.sin(phi2)  \
                 + math.cos(phi1)*math.cos(phi2)*math.cos(Phi)
        rotm23 =   math.cos(phi2)*math.sin(Phi)

        rotm31 =   math.sin(phi1)*math.sin(Phi)
        rotm32 = - math.cos(phi1)*math.sin(Phi)
        rotm33 =   math.cos(Phi)

        # to use array function, one has to note that the double [[
        # to set up the second dimension
        rotmat = np.array([[rotm11, rotm12, rotm13], [rotm21, rotm22, rotm23],  \
                          [rotm31, rotm32, rotm33]])

        return rotmat

    # @staticmethod
    # def euler2ax (eulerangles):
    #     '''
    #     角度是負數的時候不知道怎麼辦?
    #     '''
    #     [phi1, Phi, phi2]= eulerangles
    #
    #     if phi1 == 0 and Phi == 0 and phi2 == 0:
    #     # in case that there is an exception of 0, 0, 0 euler angle
    #         axis  = [0, 0, 1]
    #         angle = 0
    #         return axis, angle
    #
    #     sigma = (phi1 + phi2)/2
    #     delta = (phi1 - phi2)/2
    #     t     = math.tan(Phi/2)
    #     tau   = math.sqrt(t**2 + math.sin(sigma) * math.sin(sigma))
    #     alpha = 2 * math.atan (tau / math.cos(sigma))
    #
    #     if alpha <= math.pi:
    #         axis = [t*math.cos(delta)/tau, t*math.sin(delta)/tau, math.sin(sigma)/tau]
    #         angle= alpha
    #     else:
    #         axis = [-t*math.cos(delta)/tau, -t*math.sin(delta)/tau, -math.sin(sigma)/tau]
    #         angle= 2*math.pi - alpha
    #
    #     angle = math.degrees(angle)
    #
    #     return axis, angle

class RotMat:
    @staticmethod
    def rotm2euler(rotmat):
        '''
        one must note that in numpy, the index start also from 0;
        [[ 00  01  02 ],
         [ 10  11  12 ],
         [ 20  21  22 ]]
        '''
        if rotmat[0,0] == 1 and rotmat[1,1] == 1 and rotmat[2,2] == 1:
            phi1 = 0
            Phi  = 0
            phi2 = 0
            print ('euler angle is 0, 0, 0')

        elif abs(rotmat[2,2]) == 1:
            phi1 = math.atan2(rotmat[0,1], rotmat[0,0])
            Phi  = math.pi / 2 * (1 - rotmat[2,2])
            phi2 = 0

        else:

            zeta = 1 / math.sqrt(1 - rotmat[2,2]**2)
            phi1 = math.atan2 (rotmat[2,0]*zeta, -rotmat[2,1]*zeta)
            Phi  = math.acos  (rotmat[2,2])
            phi2 = math.atan2 (rotmat[0,2]*zeta,  rotmat[1,2]*zeta)

            # handling the atan2 property;
            # atan2 has the value domain of -pi to pi, but not useful in pi to 0.
        if phi1 < 0: phi1 = phi1 + 2*math.pi
        if phi2 < 0: phi2 = phi2 + 2*math.pi

        return [phi1, Phi, phi2]

    @staticmethod
    def rotm2ax (rotmat):
        angle = math.acos((np.trace(rotmat) - 1)/2)
        axis = [0, 0, 0]

        if angle == 0:
            axis = [0, 0 , 1]

        else:
            '''
            refer to the geometry of crystals by Bhadeshia;
            the formula is alreading being taken with P = -1.
            note the index
            '''
            axis [0] = (rotmat[1,2] - rotmat [2,1]) / 2 / math.sin(angle)
            axis [1] = (rotmat[2,0] - rotmat [0,2]) / 2 / math.sin(angle)
            axis [2] = (rotmat[0,1] - rotmat [1,0]) / 2 / math.sin(angle)

        return axis, angle

class AxisAngle:
    @staticmethod
    def axisangle2qu(axis, angle, unit = 'deg'):
        '''
        in default, the angles are in degrees
        '''
        if unit == 'deg' or 'd' in unit:
            angle = math.radians(angle)

        axis = axis / np.linalg.norm(axis)
        qu = [0.0,0.0,0.0,0.0]
        qu[0] = math.cos(angle/2)
        qu[1] = math.sin(angle/2) * axis[0]
        qu[2] = math.sin(angle/2) * axis[1]
        qu[3] = math.sin(angle/2) * axis[2]
        return qu

class OIMExport:
    @staticmethod
    def oim_or_plot (quavefer, quaveaus, directory):
        quaveausi = Qu.qu_inv(quaveaus)
        qu_or = Qu.qu_mult(quavefer, quaveausi)
        direc = directory
        direc.split('\\')

        copyfile('head.ang', directory + direc[-1] + '\\oim_or_centered.ang')
        filew = open(directory + direc[-1] + '\\oim_or_centered.ang', 'a')
        filew.write('\n')

        copyfile('head.ang', directory + direc[-1] + 'oim_or.ang')
        filew2 = open(directory + direc[-1] + 'oim_or.ang', 'a')
        filew2.write('\n')



        for x in range(24):
            oimline = '2.57290   0.31547   3.04730  ' +\
                                    str(x) + '.00000 0.00000 100 1 0 1 1.355'

            oimline_str = oimline.split()
            oimline_str2= oimline.split()

            qu = Qu.qu_mult(qu_or, cubicsymmetry[x])
            quaveaus24 = Qu.qu_mult(cubicsymmetry[x], quaveaus)


            euler = Qu.qu2euler(qu)
            oimline_str[0] = '   {0:7.5f}'.format(euler[0])
            oimline_str[1] = '{0:7.5f}'.format(euler[1])
            oimline_str[2] = '{0:7.5f}'.format(euler[2])
            line = ('  ').join(oimline_str) +'\n'
            filew.write(line)
 

            euler2= Qu.qu2euler(Qu.qu_mult(qu_or, quaveaus24))
            oimline_str2[0] = '   {0:7.5f}'.format(euler2[0])
            oimline_str2[1] = '{0:7.5f}'.format(euler2[1])
            oimline_str2[2] = '{0:7.5f}'.format(euler2[2])
            line2 = ('  ').join(oimline_str2) +'\n'
            filew2.write(line2)

        return filew

    @staticmethod
    def nw_oim_or_plot(quaveaus, directory):
    
        direc = directory
        direc.split('\\')

        copyfile('head.ang', directory + direc[-1] + 'nw_oim_or.ang')
        filew = open(directory + direc[-1] + 'nw_oim_or.ang', 'a')
        filew.write('\n')

        for x in range(12):
            oimline = '2.57290   0.31547   3.04730  ' +\
                                    str(x) + '.00000 0.00000 100 1 0 1 1.355'

            oimline_str = oimline.split()

            qu = Qu.qu_mult(nwdict[x+1], quaveaus)

            euler = Qu.qu2euler(qu)
            oimline_str[0] = '   {0:7.5f}'.format(euler[0])
            oimline_str[1] = '{0:7.5f}'.format(euler[1])
            oimline_str[2] = '{0:7.5f}'.format(euler[2])
            line = ('  ').join(oimline_str) +'\n'
            filew.write(line)

        filew.close

    @staticmethod
    def ks_oim_or_plot(quaveaus, directory):
    
        direc = directory
        direc.split('\\')

        copyfile('head.ang', directory + direc[-1] + 'ks_oim_or.ang')
        filew = open(directory + direc[-1] + 'ks_oim_or.ang', 'a')
        filew.write('\n')

        for x in range(24):
            oimline = '2.57290   0.31547   3.04730  ' +\
                                    str(x) + '.00000 0.00000 100 1 0 1 1.355'

            oimline_str = oimline.split()

            qu = Qu.qu_mult(ksdict[x+1], quaveaus)

            euler = Qu.qu2euler(qu)
            oimline_str[0] = '   {0:7.5f}'.format(euler[0])
            oimline_str[1] = '{0:7.5f}'.format(euler[1])
            oimline_str[2] = '{0:7.5f}'.format(euler[2])
            line = ('  ').join(oimline_str) +'\n'
            filew.write(line)

        filew.close
    
    @staticmethod
    def aus_oim_or_plot(quaveaus, directory):
        direc = directory
        direc.split('\\')

        copyfile('head.ang', directory + direc[-1] + 'aus_oim_or.ang')
        with open(directory + direc[-1] + 'aus_oim_or.ang', 'a') as filew:
            filew.write('\n')

            oimline = '2.57290   0.31547   3.04730  ' +\
                                    str(1) + '.00000 0.00000 100 1 0 1 1.355'

            oimline_str = oimline.split()

            euler = Qu.qu2euler(quaveaus)
            oimline_str[0] = '   {0:7.5f}'.format(euler[0])
            oimline_str[1] = '{0:7.5f}'.format(euler[1])
            oimline_str[2] = '{0:7.5f}'.format(euler[2])
            line = ('  ').join(oimline_str) +'\n'
            filew.write(line)
    
    
    @staticmethod
    def nw_variant_selecion_or_plot(quaveaus, directory, selected_variant):
    
        direc = directory
        direc.split('\\')

        copyfile('head.ang', directory + direc[-1] + 'nw_variant_selecion_or_plot.ang')
        filew = open(directory + direc[-1] + 'nw_variant_selecion_or_plot.ang', 'a')
        filew.write('\n')

        
        selected_variant = sorted(selected_variant)

        for x in selected_variant:
            oimline = '2.57290   0.31547   3.04730  ' +\
                                    str(x) + '.00000 0.00000 100 1 0 1 1.355'

            oimline_str = oimline.split()

            qu = Qu.qu_mult(nwdict[x], quaveaus)

            euler = Qu.qu2euler(qu)
            oimline_str[0] = '   {0:7.5f}'.format(euler[0])
            oimline_str[1] = '{0:7.5f}'.format(euler[1])
            oimline_str[2] = '{0:7.5f}'.format(euler[2])
            line = ('  ').join(oimline_str) +'\n'
            filew.write(line)

        filew.close


class Polefig:
    @staticmethod
    def stereo_conversion(line_vector):
        '''
        conversion of a line vector into r, theta for stereographic plots
        '''
        line_vector = line_vector / np.linalg.norm(line_vector)
        [x, y, z] = line_vector

        if z >= 0:
            R2=(1-z)/(1+z);
            S2=x*(1+R2)/2;
            T2=y*(1+R2)/2;
        else:
            R2= (1+z)/(1-z);
            S2=-x*(1+R2)/2;
            T2=-y*(1+R2)/2;

        r = math.sqrt( S2**2 + T2**2)

        if r == 0:
            theta = 0
        elif (T2 >= 0):
            theta = math.acos( S2 / r) - math.pi / 2
            # 主要原因是立體投影圖，往下為正x，往右為 +y，因此差了90度

        else:
            theta = math.acos( -S2 / r) + math.pi / 2

        theta = math.degrees(theta)

        return r, theta

    @classmethod
    def origin_polefig (cls, quavefer, quaveaus, directory):
        quaveaus = Qu.qu_inv(quaveaus)
        qu_or = Qu.qu_mult(quavefer, quaveaus)
        direc = directory
        direc.split('\\')
        # direc的最後一個就是檔案名稱
        filew = open(directory + direc[-1] + '\\origin_or_100.txt', 'w')
        filew110 = open(directory + direc[-1] + '\\origin_or_110.txt', 'w')
        filew111 = open(directory + direc[-1] + '\\origin_or_111.txt', 'w')

        for i in range(24):
            euler = Qu.qu2euler(Qu.qu_mult(qu_or, cubicsymmetry[i]))
            rotmat= Euler.euler2rotm(euler)

            for y in range(3):
                r, theta = cls.stereo_conversion(rotmat[y, :])
                filew.write('{0:4.4f}\t{1:4.4f}\n'.format(theta, r))

            for aus110 in [2,6,8,12,10,20]:
                fer111 = np.matmul(aus110dict[aus110], rotmat)
                r, theta = cls.stereo_conversion(fer111)
                filew111.write('{0:4.4f}\t{1:4.4f}\n'.format(theta, r))


            for aus111 in [1,7,13,19]:
                fer110 = np.matmul(aus111dict[aus111], rotmat)
                r, theta = cls.stereo_conversion(fer110)
                filew110.write('{0:4.4f}\t{1:4.4f}\n'.format(theta, r))

class SpecialTrans:
    @staticmethod
    def variant24_convert(variants_counts):
        '''
        最多的轉成V1，有24組
        '''
        varmax = np.argmax(variants_counts)
        print ('Number of max variant:', varmax+1)
        variant_temp = [0]*24
        print (variants_counts)
    

        if varmax > 0:
            print(list(variant_convert_table[0,:]))
            column = list(variant_convert_table[0,:]).index(varmax+1)
            print ('Number of column in convert table:', column)
            print('Transform max variant to variant_1:')
            for i in range(24):
                variant_temp[i] = variants_counts\
                            [variant_convert_table[i,column] - 1]

            variants_counts = variant_temp
        
        print(variants_counts)
        print('Transformation finished!')

        return variants_counts

class AngleCalc:
    @staticmethod
    def cosVector(x,y):
        if(len(x)!=len(y)):
            print('error input,x and y is not in the same space')
            return
        result1=0.0
        result2=0.0
        result3=0.0
        for i in range(len(x)):
            result1+=x[i]*y[i]   #sum(X*Y)
            result2+=x[i]**2     #sum(X*X)
            result3+=y[i]**2     #sum(Y*Y)
        #print(result1)
        #print(result2)
        #print(result3)
        return float(result1/((result2*result3)**0.5))


if __name__ == '__main__':
    Eu = Euler()
    Qu = Qu()
    convert = np.empty([24,24])
    minnum = np.zeros(24)

    maxnum = np.empty([24,24])
    for (i, gtcor) in GTitoGT1_cubic_corres.items():
        for (j, gt) in gtdict.items():
            for (k, gt2) in gtdict.items():
                minnum[k-1] = Qu.qu_angle(gt2, Qu.qu_mult(gt, gtcor), 'yes')
            convert[j-1, i-1] = (np.argmin(minnum)+1)

    sortin = np.argsort(convert[0,:])
    print(sortin)
