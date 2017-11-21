# My class travel carbon footprint calculator
# Version 5.0

# import required libraries
import numpy as np
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#######################################################
# functions                                           #
#######################################################

def footprint(emcar,emtrain,emplane):
    # open the file
    fin = open('distances2017.csv','r')

    # skip the first line
    fin.readline()

    # initialize the variable
    foot = []

    # loop over input lines
    for line in fin:
        # get distance, frequency, and transport mode (name can be skipped)
        _, dist, freq, mode = line.strip().split(',')
        # convert required values to floats 
        dist = float(dist)
        freq = float(freq)
        # computing total distance 
        totdist = dist * freq
       
        # establish the emission coefficients
        if mode == 'car':
            em = emcar
        elif mode == 'train':
            em = emtrain
        elif mode == 'airplane':
            em = emplane
       
        # compute the footprint for each student
        foottemp =  totdist * em
        # append to the list
        foot.append(foottemp)
         
    # calculating average footprint and median
    mean = np.mean(foot)
    median = np.median(foot)

    return {'mean': mean, 'median': median}

#######################################################
coef1=[10,2.5,200]
    
emissions1=footprint(coef1[0],coef1[1],coef1[2])

random_dic=[]
random_tup=[]
for i in range(1000):
   tuple_random_number=np.random.uniform(5,50),np.random.uniform(1,10),np.random.uniform(199,355)
   
   item_to_list=footprint(tuple_random_number[0],tuple_random_number[1],tuple_random_number[2])
   
   random_dic.append(item_to_list) #mean and median
   random_tup.append(tuple_random_number)

means=[]
for i in range(len(random_dic)):
    means.append(random_dic[i]['mean'])

n, bins, patches = plt.hist(means, 50, normed=1, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=1)
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'Histogram')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()
#test1


