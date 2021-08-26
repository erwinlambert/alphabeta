from parameters import *

parms = parameters()

for i in parms.keys():
    exec(i+"=parms['"+i+"']['basicval']")
#    print i, parms[i]['basicval']
