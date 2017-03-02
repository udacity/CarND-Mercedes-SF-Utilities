
# coding: utf-8

# # Loading from the ekf input file

# In[42]:

import plotly.offline as py
from plotly.graph_objs import *
import pandas as pd
import math
py.init_notebook_mode()

my_cols=['c1','c2','c3','c4','c5']
with open('obj_pose-laser-radar-ekf-input.txt') as f:
    table = pd.read_table(f, sep='\t', header=None, names=my_cols, lineterminator='\n')
    
table_vis_in = pd.DataFrame();

def f_x(row):
    if(row[0]=='L'):
        return row[1];
    else:
        ro = row[1]
        theta = row[2]
        return ro*math.cos(theta)

def f_y(row):
    if(row[0]=='L'):
        return row[2];
    else:
        ro = row[1]
        theta = row[2]
        return ro*math.sin(theta)

table_vis_in['c1'] = table['c1']
table_vis_in['c2'] = table.apply(f_x, axis=1)
table_vis_in['c3'] = table.apply(f_y, axis=1)


# In[43]:

table_vis_in


# In[ ]:




# In[ ]:




# In[38]:

table_vis


# In[ ]:

x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
x_rev = x[::-1]

# Line 1
y1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
y1_upper = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
y1_lower = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y1_lower = y1_lower[::-1]

# Line 2
y2 = [5, 2.5, 5, 7.5, 5, 2.5, 7.5, 4.5, 5.5, 5]
y2_upper = [5.5, 3, 5.5, 8, 6, 3, 8, 5, 6, 5.5]
y2_lower = [4.5, 2, 4.4, 7, 4, 2, 7, 4, 5, 4.75]
y2_lower = y2_lower[::-1]

# Line 3
y3 = [10, 8, 6, 4, 2, 0, 2, 4, 2, 0]
y3_upper = [11, 9, 7, 5, 3, 1, 3, 5, 3, 1]
y3_lower = [9, 7, 5, 3, 1, -.5, 1, 3, 1, -1]
y3_lower = y3_lower[::-1]

trace1 = Scatter(
    x=x+x_rev,
    y=y1_upper+y1_lower,
    fill='tozerox',
    fillcolor='rgba(0,100,80,0.2)',
    line=Line(color='transparent'),
    showlegend=False,
    name='Fair',
)
trace2 = Scatter(
    x=x+x_rev,
    y=y2_upper+y2_lower,
    fill='tozerox',
    fillcolor='rgba(0,176,246,0.2)',
    line=Line(color='transparent'),
    name='Premium',
    showlegend=False,
)
trace3 = Scatter(
    x=x+x_rev,
    y=y3_upper+y3_lower,
    fill='tozerox',
    fillcolor='rgba(231,107,243,0.2)',
    line=Line(color='transparent'),
    showlegend=False,
    name='Fair',
)
trace4 = Scatter(
    x=x,
    y=y1,
    line=Line(color='rgb(0,100,80)'),
    mode='lines',
    name='Fair',
)
trace5 = Scatter(
    x=x,
    y=y2,
    line=Line(color='rgb(0,176,246)'),
    mode='lines',
    name='Premium',
)
trace6 = Scatter(
    x=x,
    y=y3,
    line=Line(color='rgb(231,107,243)'),
    mode='lines',
    name='Ideal',
)

data = Data([trace1, trace2, trace3, trace4, trace5, trace6])

layout = Layout(
    paper_bgcolor='rgb(255,255,255)',
    plot_bgcolor='rgb(229,229,229)',
    xaxis=XAxis(
        gridcolor='rgb(255,255,255)',
        range=[1,10],
        showgrid=True,
        showline=False,
        showticklabels=True,
        tickcolor='rgb(127,127,127)',
        ticks='outside',
        zeroline=False
    ),
    yaxis=YAxis(
        gridcolor='rgb(255,255,255)',
        showgrid=True,
        showline=False,
        showticklabels=True,
        tickcolor='rgb(127,127,127)',
        ticks='outside',
        zeroline=False
    ),
)
fig = Figure(data=data, layout=layout)
py.iplot(fig, filename= 'shaded_lines')

