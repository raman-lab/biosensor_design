#!/usr/bin/env python

import plotly.plotly as py
from plotly.graph_objs import *
from plotly.tools import FigureFactory as FF

data = [['', 'Emma', 'Isabella', 'Ava', 'Olivia', 'Sophia', 'row-sum'],
        ['Emma', 16, 3, 28, 0, 18, 65],
        ['Isabella', 18, 0, 12, 5, 29, 64],
        ['Ava', 9, 11, 17, 27, 0, 64],
        ['Olivia', 19, 0, 31, 11, 12, 73],
        ['Sophia', 23, 17, 10, 0, 34, 84]]

table = FF.create_table(data, index=True)
py.iplot(table, filename='Data-Table')

import numpy as np

matrix=np.array([[16,  3, 28,  0, 18],
                 [18,  0, 12,  5, 29],
                 [ 9, 11, 17, 27,  0],
                 [19,  0, 31, 11, 12],
                 [23, 17, 10,  0, 34]], dtype=int)

def check_data(data_matrix):
    L, M=data_matrix.shape
    if L!=M:
        raise ValueError('Data array must have (n,n) shape')
    return L

L=check_data(matrix)

PI=np.pi

def moduloAB(x, a, b): #maps a real number onto the unit circle identified with
                       #the interval [a,b), b-a=2*PI
        if a>=b:
            raise ValueError('Incorrect interval ends')
        y=(x-a)%(b-a)
        return y+b if y<0 else y+a

def test_2PI(x):
    return 0<= x <2*PI

row_sum=[np.sum(matrix[k,:]) for k in range(L)]

#set the gap between two consecutive ideograms
gap=2*PI*0.005
ideogram_length=2*PI*np.asarray(row_sum)/sum(row_sum)-gap*np.ones(L)

def get_ideogram_ends(ideogram_len, gap):
    ideo_ends=[]
    left=0
    for k in range(len(ideogram_len)):
        right=left+ideogram_len[k]
        ideo_ends.append([left, right])
        left=right+gap
    return ideo_ends

ideo_ends=get_ideogram_ends(ideogram_length, gap)
ideo_ends

def make_ideogram_arc(R, phi, a=50):
    # R is the circle radius
    # phi is the list of ends angle coordinates of an arc
    # a is a parameter that controls the number of points to be evaluated on an arc
    if not test_2PI(phi[0]) or not test_2PI(phi[1]):
        phi=[moduloAB(t, 0, 2*PI) for t in phi]
    length=(phi[1]-phi[0])% 2*PI
    nr=5 if length<=PI/4 else int(a*length/PI)

    if phi[0] < phi[1]:
        theta=np.linspace(phi[0], phi[1], nr)
    else:
        phi=[moduloAB(t, -PI, PI) for t in phi]
        theta=np.linspace(phi[0], phi[1], nr)
    return R*np.exp(1j*theta)

z=make_ideogram_arc(1.3, [11*PI/6, PI/17])
print z

labels=['Emma', 'Isabella', 'Ava', 'Olivia', 'Sophia']
ideo_colors=['rgba(244, 109, 67, 0.75)',
             'rgba(253, 174, 97, 0.75)',
             'rgba(254, 224, 139, 0.75)',
             'rgba(217, 239, 139, 0.75)',
             'rgba(166, 217, 106, 0.75)']

def map_data(data_matrix, row_value, ideogram_length):
    mapped=np.zeros(data_matrix.shape)
    for j  in range(L):
        mapped[:, j]=ideogram_length*data_matrix[:,j]/row_value
    return mapped

mapped_data=map_data(matrix, row_sum, ideogram_length)
mapped_data


idx_sort=np.argsort(mapped_data, axis=1)
idx_sort

def make_ribbon_ends(mapped_data, ideo_ends,  idx_sort):
    L=mapped_data.shape[0]
    ribbon_boundary=np.zeros((L,L+1))
    for k in range(L):
        start=ideo_ends[k][0]
        ribbon_boundary[k][0]=start
        for j in range(1,L+1):
            J=idx_sort[k][j-1]
            ribbon_boundary[k][j]=start+mapped_data[k][J]
            start=ribbon_boundary[k][j]
    return [[(ribbon_boundary[k][j],ribbon_boundary[k][j+1] ) for j in range(L)] for k in range(L)]

ribbon_ends=make_ribbon_ends(mapped_data, ideo_ends,  idx_sort)
print 'ribbon ends starting from the ideogram[2]\n', ribbon_ends[2]

def control_pts(angle, radius):
    #angle is a  3-list containing angular coordinates of the control points b0, b1, b2
    #radius is the distance from b1 to the  origin O(0,0)

    if len(angle)!=3:
        raise InvalidInputError('angle must have len =3')
    b_cplx=np.array([np.exp(1j*angle[k]) for k in range(3)])
    b_cplx[1]=radius*b_cplx[1]
    return zip(b_cplx.real, b_cplx.imag)

def ctrl_rib_chords(l, r, radius):
    # this function returns a 2-list containing control poligons of the two quadratic Bezier
    #curves that are opposite sides in a ribbon
    #l (r) the list of angular variables of the ribbon arc ends defining
    #the ribbon starting (ending) arc
    # radius is a common parameter for both control polygons
    if len(l)!=2 or len(r)!=2:
        raise ValueError('the arc ends must be elements in a list of len 2')
    return [control_pts([l[j], (l[j]+r[j])/2, r[j]], radius) for j in range(2)]